/*

   This file is part of the regenie software package.

   Copyright (c) 2020-2024 Joelle Mbatchou, Andrey Ziyatdinov & Jonathan Marchini

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in all
   copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.

*/

#include <fstream>
#include <iostream>
#include <sstream>

#include "bgen/bgen.h"

#include "BgenReader.hpp"
#include "S3_Utils.hpp"

namespace fs = boost::filesystem;

// Sequential variant iteration. Implemented in bgen-limix's src/variant.c and
// exported, but not yet declared in its installed headers.
extern "C" {
  struct bgen_variant* bgen_variant_begin(struct bgen_file*, int* error);
  struct bgen_variant* bgen_variant_next(struct bgen_file*, int* error);
  void bgen_variant_destroy(struct bgen_variant const*);
}

std::string bgen_metafile_path(std::string const& bgen_file){
  return bgen_file + ".metafile";
}

namespace {

bool is_readable(std::string const& p){
  std::ifstream probe(p, std::ios::in | std::ios::binary);
  return probe.good();
}

bool allow_remote_metafile_build = false;
bool force_scan_only = false;
std::string explicit_metafile_path;

// A stable per-input key, so a cached metafile is reused across runs rather
// than rebuilt. Size is included so an input replaced in place is not matched.
std::string cache_key(std::string const& path, int64_t size){
  std::hash<std::string> h;
  std::ostringstream o;
  o << std::hex << h(path) << "_" << size;
  return o.str();
}

// Directory to build a metafile in when it cannot live beside the input.
fs::path metafile_cache_dir(){

  char const* env[] = {"REGENIE_CACHE_DIR", "XDG_CACHE_HOME", "TMPDIR"};
  for(size_t i = 0; i < 3; i++){
    char const* v = getenv(env[i]);
    if(v == nullptr || *v == '\0') continue;
    fs::path d(v);
    if(i == 1) d /= "regenie";            // XDG_CACHE_HOME is shared
    boost::system::error_code ec;
    fs::create_directories(d, ec);
    if(!ec && fs::is_directory(d)) return d;
  }

  boost::system::error_code ec;
  fs::path d = fs::temp_directory_path(ec) / "regenie";
  if(!ec){
    fs::create_directories(d, ec);
    if(!ec) return d;
  }
  return fs::current_path();
}

bool dir_is_writable(fs::path const& dir){
  boost::system::error_code ec;
  fs::path probe = dir / ".regenie_write_test";
  { std::ofstream f(probe.string()); if(!f.good()) return false; }
  fs::remove(probe, ec);
  return true;
}

} // anonymous namespace

BgenParser::~BgenParser(){
  if(mfile != nullptr) bgen_metafile_close(mfile);
  if(bfile != nullptr) bgen_file_close(bfile);
}

void BgenParser::set_allow_remote_metafile_build(bool allow){
  allow_remote_metafile_build = allow;
}

void BgenParser::set_metafile_path(std::string const& path){
  explicit_metafile_path = path;
}

void BgenParser::set_force_scan_only(bool force){
  force_scan_only = force;
}

// Resolution order: an explicit path, then one beside the input, then one in
// the cache, then build one. Building never writes beside a remote input, and
// falls back to the cache when the input's directory is read-only.
std::string BgenParser::resolve_metafile_path(bool& must_create) const {

  if(!explicit_metafile_path.empty()){
    must_create = !is_readable(explicit_metafile_path);
    return explicit_metafile_path;
  }

  std::string const beside = bgen_metafile_path(path);
  bool const remote = is_remote_path(path);

  if(!remote && is_readable(beside)){ must_create = false; return beside; }

  fs::path const cached =
    metafile_cache_dir() / (cache_key(path, file_size) + ".metafile");
  if(is_readable(cached.string())){ must_create = false; return cached.string(); }

  must_create = true;

  if(!remote){
    fs::path dir = fs::path(beside).parent_path();
    if(dir.empty()) dir = fs::current_path();
    if(dir_is_writable(dir)) return beside;
  }

  return cached.string();
}

void BgenParser::load_metafile(){

  if(force_scan_only && !is_remote_path(path)){ scan_variants(); return; }

  bool must_create = false;
  std::string const mpath = resolve_metafile_path(must_create);

  if(!must_create) mfile = bgen_metafile_open(mpath.c_str());

  if(mfile != nullptr){ read_metafile(mpath); return; }

  // No metafile, or one written by an incompatible version.
  if(!is_remote_path(path)){

    // Cheap, local, and pays for itself once a run repeats or an --extract
    // subset makes random access from the metafile worthwhile: build one
    // unless the caller asked not to.
    mfile = bgen_metafile_create(bfile, mpath.c_str(), 1, 0);
    if(mfile != nullptr){ read_metafile(mpath); return; }

    scan_variants();
    return;
  }

  // Enumerating a remote file costs either one request per variant or the
  // whole object in egress, depending on read-ahead. Neither is something to
  // spend on the user's behalf.
  if(!allow_remote_metafile_build)
    throw "no metafile found for remote file : " + path + "\n"
      "       Building one requires walking every variant, which costs either\n"
      "       one request per variant or a full download, and can take hours\n"
      "       for a large file. Generate\n"
      "       " + bgen_metafile_path(path) + "\n"
      "       alongside the data instead -- doing that in the same region as\n"
      "       the bucket avoids the egress and the round trips entirely.\n"
      "       Pass --allow-remote-metafile-build to accept the cost and build\n"
      "       it now; the result is cached and paid for only once.";

  std::cerr << "WARNING: walking every variant of " << path
            << " to build a metafile; this may take hours for a large file.\n";

  mfile = bgen_metafile_create(bfile, mpath.c_str(), 1, 0);
  if(mfile == nullptr)
    throw "cannot create bgen metafile : " + mpath;

  read_metafile(mpath);
}

// Walks the file, collecting the same metadata a metafile holds. Genotype
// blocks are seeked over, not read.
void BgenParser::scan_variants(){

  variants.clear();
  variants.reserve(bgen_file_nvariants(bfile));

  int err = 0;
  for(bgen_variant* v = bgen_variant_begin(bfile, &err); v != nullptr;
      v = bgen_variant_next(bfile, &err)){

    variants.emplace_back();
    variant_meta& m = variants.back();

    m.genotype_offset = v->genotype_offset;
    m.position = v->position;
    m.chromosome.assign(v->chrom->data, v->chrom->length);
    m.rsid.assign(v->rsid->data, v->rsid->length);

    m.alleles.resize(v->nalleles);
    for(uint16_t ia = 0; ia < v->nalleles; ia++)
      m.alleles[ia].assign(v->allele_ids[ia]->data, v->allele_ids[ia]->length);

    bgen_variant_destroy(v);

    if(err) throw "error reading variant " + std::to_string(variants.size())
      + " of " + path;
  }

  if(err) throw "error walking variants of " + path;

  n_variants = static_cast<uint32_t>(variants.size());
}

void BgenParser::read_metafile(std::string const& mpath){

  n_variants = bgen_metafile_nvariants(mfile);
  variants.resize(n_variants);

  uint32_t const nparts = bgen_metafile_npartitions(mfile);
  size_t ivar = 0;

  for(uint32_t ipart = 0; ipart < nparts; ipart++){

    bgen_partition const* part = bgen_metafile_read_partition(mfile, ipart);
    if(part == nullptr)
      throw "cannot read partition " + std::to_string(ipart) + " of " + mpath;

    uint32_t const nvars = bgen_partition_nvariants(part);
    for(uint32_t iv = 0; iv < nvars; iv++, ivar++){

      bgen_variant const* v = bgen_partition_get_variant(part, iv);
      variant_meta& m = variants[ivar];

      m.genotype_offset = v->genotype_offset;
      m.position = v->position;
      m.chromosome.assign(v->chrom->data, v->chrom->length);
      m.rsid.assign(v->rsid->data, v->rsid->length);

      m.alleles.resize(v->nalleles);
      for(uint16_t ia = 0; ia < v->nalleles; ia++)
        m.alleles[ia].assign(v->allele_ids[ia]->data, v->allele_ids[ia]->length);
    }

    bgen_partition_destroy(part);
  }
}

void BgenParser::open(std::string const& filename){

  path = filename;

  // Only used to key the index cache; remote inputs simply get 0.
  if(!is_remote_path(filename)){
    boost::system::error_code ec;
    uintmax_t const sz = fs::file_size(filename, ec);
    if(!ec) file_size = static_cast<int64_t>(sz);
  }

  bfile = bgen_file_open(filename.c_str());
  if(bfile == nullptr) throw "cannot open bgen file : " + filename;

  n_samples = bgen_file_nsamples(bfile);
  layout = bgen_file_layout(bfile);
  compression = bgen_file_compression(bfile);
  have_sample_ids = bgen_file_contain_samples(bfile);

  load_metafile();

  if(have_sample_ids){
    bgen_samples* s = bgen_file_read_samples(bfile);
    if(s != nullptr){
      sample_ids.resize(n_samples);
      for(uint32_t i = 0; i < n_samples; i++){
        bgen_string const* id = bgen_samples_get(s, i);
        sample_ids[i].assign(id->data, id->length);
      }
      bgen_samples_destroy(s);
    }
  }

  if(sample_ids.empty()){ // anonymous samples get a placeholder, as before
    sample_ids.resize(n_samples);
    for(uint32_t i = 0; i < n_samples; i++)
      sample_ids[i] = "(unknown_sample_" + std::to_string(i + 1) + ")";
  }

  if(n_variants > 0){
    bgen_genotype* gt = bgen_file_open_genotype(bfile, variants[0].genotype_offset);
    if(gt == nullptr) throw "cannot read first genotype block of : " + filename;
    nbits = bgen_genotype_nbits(gt);
    bgen_genotype_close(gt);
  }

  cursor = 0;
  last_read = 0;
}

std::string BgenParser::summarise() const {

  std::ostringstream o;
  o << "   -summary : bgen file ("
    << (layout == 2 ? "v1.2 layout" : "v1.1 layout") << ", ";
  if(compression == 1) o << "zlib ";
  else if(compression == 2) o << "zstd ";
  o << (compression ? "compressed" : "uncompressed") << ")"
    << " with " << n_samples << " " << (have_sample_ids ? "named" : "anonymous")
    << " samples and " << n_variants << " variants";

  return o.str();
}

uint64_t BgenParser::get_position() const {
  if(cursor >= variants.size()) return 0;
  return variants[cursor].genotype_offset;
}

void BgenParser::build_offset_map(){
  offset_to_index.reserve(variants.size());
  for(size_t i = 0; i < variants.size(); i++)
    offset_to_index[ variants[i].genotype_offset ] = i;
}

void BgenParser::jumpto(uint64_t genotype_offset){

  if(offset_to_index.empty() && !variants.empty()) build_offset_map();

  auto it = offset_to_index.find(genotype_offset);
  if(it == offset_to_index.end())
    throw "bgen offset not present in index : " + std::to_string(genotype_offset);

  cursor = it->second;
}

bool BgenParser::read_variant(std::string* chromosome, uint32_t* position,
    std::string* rsid, std::vector<std::string>* alleles){

  if(cursor >= variants.size()) return false;

  variant_meta const& m = variants[cursor];
  *chromosome = m.chromosome;
  *position = m.position;
  *rsid = m.rsid;
  *alleles = m.alleles;

  last_read = cursor++;
  return true;
}

void BgenParser::read_probs(std::vector<std::vector<double>>* probs){

  if(last_read >= variants.size())
    throw std::string("no bgen variant has been read");

  bgen_genotype* gt = bgen_file_open_genotype(bfile, variants[last_read].genotype_offset);
  if(gt == nullptr)
    throw "cannot read genotypes for variant : " + variants[last_read].rsid;

  unsigned const ncombs = bgen_genotype_ncombs(gt);

  std::vector<double> flat(static_cast<size_t>(n_samples) * ncombs);
  if(bgen_genotype_read(gt, flat.data()) != 0){
    bgen_genotype_close(gt);
    throw "cannot decode genotypes for variant : " + variants[last_read].rsid;
  }

  probs->resize(n_samples);
  for(uint32_t i = 0; i < n_samples; i++){
    (*probs)[i].resize(ncombs);
    bool const missing = bgen_genotype_missing(gt, i);
    for(unsigned j = 0; j < ncombs; j++)
      (*probs)[i][j] = missing ? -1 : flat[static_cast<size_t>(i) * ncombs + j];
  }

  bgen_genotype_close(gt);
}
