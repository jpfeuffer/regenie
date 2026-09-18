// Micro-benchmark for the genotype read path: replicates readChunkFromBGEN's
// exact access pattern (seek to genotype offset, read two 4-byte sizes, read
// the compressed block) against std::ifstream and against bgen-limix's
// stream_handle, on the same file and the same offsets.
//
// The end-to-end step 2 run is compute-dominated (~110 ms/block of which the
// read is a small part), so it cannot resolve an I/O change; this can.

#include <bgen/bgen.h>
#include <bgen/s3stream.h>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <numeric>
#include <string>
#include <vector>

extern "C" {
struct bgen_variant* bgen_variant_begin(struct bgen_file*, int*);
struct bgen_variant* bgen_variant_next(struct bgen_file*, int*);
void bgen_variant_destroy(struct bgen_variant const*);
}

typedef unsigned char uchar;

static double now_s() {
  using namespace std::chrono;
  return duration<double>(steady_clock::now().time_since_epoch()).count();
}

// Mirrors readChunkFromBGEN() as it is today.
static uint64_t read_ifstream(std::istream* bfile,
                              std::vector<uint64_t> const& indices) {
  uint64_t total = 0;
  std::vector<uchar> geno_block;
  for (size_t i = 0; i < indices.size(); ++i) {
    uint32_t size1 = 0, size2 = 0;
    bfile->seekg(indices[i]);
    bfile->read(reinterpret_cast<char*>(&size1), 4);
    bfile->read(reinterpret_cast<char*>(&size2), 4);
    geno_block.resize(size1 - 4);
    bfile->read(reinterpret_cast<char*>(geno_block.data()), size1 - 4);
    total += size1;
  }
  return total;
}

// The same pattern through bgen-limix's stream_handle, which is what the
// migration switches to (fopen/fread/fseeko locally, ranged HTTP remotely).
static uint64_t read_handle(stream_handle* h,
                            std::vector<uint64_t> const& indices) {
  uint64_t total = 0;
  std::vector<uchar> geno_block;
  for (size_t i = 0; i < indices.size(); ++i) {
    uint32_t size1 = 0, size2 = 0;
    stream_handle_seek(h, static_cast<int64_t>(indices[i]), SEEK_SET);
    stream_handle_read(h, &size1, 4);
    stream_handle_read(h, &size2, 4);
    geno_block.resize(size1 - 4);
    stream_handle_read(h, geno_block.data(), size1 - 4);
    total += size1;
  }
  return total;
}

int main(int argc, char** argv) {
  if (argc < 2) {
    fprintf(stderr, "usage: bench_read <file.bgen> [reps]\n");
    return 2;
  }
  char const* path = argv[1];
  int reps = (argc > 2) ? atoi(argv[2]) : 5;

  // Collect genotype offsets the same way regenie does, via the metafile.
  bgen_file* bf = bgen_file_open(path);
  if (!bf) { fprintf(stderr, "cannot open %s\n", path); return 1; }

  std::string mpath = std::string(path) + ".metafile";
  bgen_metafile* mf = bgen_metafile_open(mpath.c_str());
  if (!mf) { fprintf(stderr, "cannot open %s\n", mpath.c_str()); return 1; }

  std::vector<uint64_t> indices;
  indices.reserve(bgen_metafile_nvariants(mf));
  uint32_t npart = bgen_metafile_npartitions(mf);
  for (uint32_t p = 0; p < npart; ++p) {
    bgen_partition const* part = bgen_metafile_read_partition(mf, p);
    uint32_t n = bgen_partition_nvariants(part);
    for (uint32_t i = 0; i < n; ++i)
      indices.push_back(bgen_partition_get_variant(part, i)->genotype_offset);
    bgen_partition_destroy(part);
  }
  bgen_metafile_close(mf);
  bgen_file_close(bf);

  printf("file    : %s\n", path);
  printf("variants: %zu\n", indices.size());
  printf("reps    : %d\n\n", reps);

  std::vector<double> t_if, t_hd;
  uint64_t bytes_if = 0, bytes_hd = 0;

  for (int r = 0; r < reps; ++r) {
    {
      std::ifstream f(path, std::ios::in | std::ios::binary);
      if (!f.good()) { fprintf(stderr, "ifstream open failed\n"); return 1; }
      double t0 = now_s();
      bytes_if = read_ifstream(&f, indices);
      t_if.push_back(now_s() - t0);
    }
    {
      stream_handle* h = stream_handle_open(path);
      if (!h) { fprintf(stderr, "stream_handle_open failed\n"); return 1; }
      double t0 = now_s();
      bytes_hd = read_handle(h, indices);
      t_hd.push_back(now_s() - t0);
      stream_handle_close(h);
    }
    printf("  rep %d: ifstream %7.3f s   stream_handle %7.3f s\n",
           r + 1, t_if.back(), t_hd.back());
  }

  if (bytes_if != bytes_hd) {
    printf("\nMISMATCH: byte totals differ (%llu vs %llu)\n",
           (unsigned long long)bytes_if, (unsigned long long)bytes_hd);
    return 1;
  }

  std::sort(t_if.begin(), t_if.end());
  std::sort(t_hd.begin(), t_hd.end());
  double best_if = t_if.front(), best_hd = t_hd.front();
  double med_if = t_if[t_if.size() / 2], med_hd = t_hd[t_hd.size() / 2];

  printf("\nbytes read : %llu (identical both paths)\n",
         (unsigned long long)bytes_if);
  printf("ifstream      best %7.3f s  median %7.3f s  %8.1f MB/s\n",
         best_if, med_if, bytes_if / best_if / 1e6);
  printf("stream_handle best %7.3f s  median %7.3f s  %8.1f MB/s\n",
         best_hd, med_hd, bytes_hd / best_hd / 1e6);
  printf("delta (median, handle vs ifstream): %+.2f%%\n",
         (med_hd - med_if) / med_if * 100.0);
  return 0;
}
