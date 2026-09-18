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

#ifndef BGEN_READER_H
#define BGEN_READER_H

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include <boost/filesystem.hpp>

struct bgen_file;
struct bgen_metafile;

// Reader for BGEN files, backed by the bgen-limix C API.
//
// Variant metadata and addresses come from bgen-limix's metafile index, which
// is created next to the BGEN file on first use. The older .bgi index is a
// SQLite database and is no longer read.
//
// Positions handed out by get_position() and accepted by jumpto() are
// genotype-block offsets (bgen_variant::genotype_offset), i.e. they point past
// the variant's identifying data straight at the genotype block. regenie's own
// fast-path decoder seeks to exactly that point.
class BgenParser {

  public:
    BgenParser() = default;
    ~BgenParser();

    BgenParser(const BgenParser&) = delete;
    BgenParser& operator=(const BgenParser&) = delete;

    void open(std::string const& filename);

    std::string summarise() const;

    bool get_layout() const { return layout == 2; }        // v1.2/v1.3
    bool get_compression() const { return compression == 1; } // zlib, else zstd

    // Bits per probability in the genotype blocks.
    uint32_t get_nbits() const { return nbits; }

    int number_of_samples() const { return n_samples; }
    int number_of_variants() const { return n_variants; }

    // Genotype-block offset of the variant the cursor is sitting on. Returns 0
    // once the cursor has run off the end.
    uint64_t get_position() const;

    void jumpto(uint64_t genotype_offset);

    template <typename Setter>
      void get_sample_ids(Setter setter) {
        for(size_t i = 0; i < sample_ids.size(); i++) setter( sample_ids[i] );
      }

    // Reads the variant at the cursor and advances it. Returns false at the end
    // of the file, matching the previous sequential-iteration contract.
    bool read_variant(std::string* chromosome, uint32_t* position,
        std::string* rsid, std::vector<std::string>* alleles);

    // Probabilities for the variant read by the last read_variant(), laid out
    // one vector per sample. Missing genotypes are encoded as -1.
    void read_probs(std::vector<std::vector<double>>* probs);

    void ignore_probs() {} // genotype blocks are only touched on demand

  private:
    // Variant metadata, copied out of the metafile partition so the partition
    // can be released rather than kept open for the life of the reader.
    struct variant_meta {
      uint64_t genotype_offset = 0;
      uint32_t position = 0;
      std::string chromosome, rsid;
      std::vector<std::string> alleles;
    };

    void load_index();
    void build_offset_map();

    // Where to read or build the index. Never returns a path beside a remote
    // input, and falls back to a cache directory when the input's directory is
    // not writable.
    std::string resolve_index_path(bool& must_create) const;

    bgen_file* bfile = nullptr;
    bgen_metafile* mfile = nullptr;

    std::string path;
    int64_t file_size = 0;
    uint32_t layout = 0, compression = 0, nbits = 0;
    uint32_t n_samples = 0, n_variants = 0;
    bool have_sample_ids = false;

    std::vector<std::string> sample_ids;
    std::vector<variant_meta> variants;

    // Only needed by the random-access paths, so built on first jumpto().
    std::unordered_map<uint64_t, size_t> offset_to_index;

    size_t cursor = 0;      // variant the next read_variant() will return
    size_t last_read = 0;   // variant the last read_variant() returned
};

// Path of the metafile index for a BGEN file.
std::string bgen_metafile_path(std::string const& bgen_file);

#endif
