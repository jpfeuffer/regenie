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
#include <sstream>

#include "bgen/bgen.h"

#include "BgenReader.hpp"

std::string bgen_metafile_path(std::string const& bgen_file){
  return bgen_file + ".metafile";
}

BgenParser::~BgenParser(){
  if(mfile != nullptr) bgen_metafile_close(mfile);
  if(bfile != nullptr) bgen_file_close(bfile);
}

void BgenParser::load_index(){

  std::string const mpath = bgen_metafile_path(path);

  std::ifstream probe(mpath, std::ios::in | std::ios::binary);
  if(probe.good()){
    probe.close();
    mfile = bgen_metafile_open(mpath.c_str());
  }

  if(mfile == nullptr){ // absent or unreadable (e.g. written by an older version)
    mfile = bgen_metafile_create(bfile, mpath.c_str(), 1, 0);
    if(mfile == nullptr)
      throw "cannot create bgen index file : " + mpath;
  }

  n_variants = bgen_metafile_nvariants(mfile);
  variants.resize(n_variants);

  uint32_t const nparts = bgen_metafile_npartitions(mfile);
  size_t ivar = 0;

  for(uint32_t ipart = 0; ipart < nparts; ipart++){

    bgen_partition const* part = bgen_metafile_read_partition(mfile, ipart);
    if(part == nullptr)
      throw "cannot read partition " + std::to_string(ipart) + " of " + bgen_metafile_path(path);

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

  bfile = bgen_file_open(filename.c_str());
  if(bfile == nullptr) throw "cannot open bgen file : " + filename;

  n_samples = bgen_file_nsamples(bfile);
  layout = bgen_file_layout(bfile);
  compression = bgen_file_compression(bfile);
  have_sample_ids = bgen_file_contain_samples(bfile);

  load_index();

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
