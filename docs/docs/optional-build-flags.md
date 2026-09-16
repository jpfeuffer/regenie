## Optional build flags

**regenie**'s CMake build exposes four optional flags that unlock
additional backends or features. None are required for a functional
build — the defaults give a portable binary that depends only on the
bundled/required libraries.

Pass any flag to CMake with `-D<FLAG>=ON`, e.g.:

```bash
cmake -DWITH_OPENBLAS=ON ..
```

---

### `WITH_MKL` — Intel Math Kernel Library

**Default:** OFF

Replaces Eigen's built-in linear algebra kernels with Intel's
[MKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html)
BLAS/LAPACK routines by setting `EIGEN_USE_BLAS` and
`EIGEN_USE_LAPACKE`.

Without this flag Eigen performs all dense matrix operations in pure
C++ template code. With MKL, those calls are dispatched to
hand-tuned, AVX-512–capable routines that can be **2–10× faster** on
Intel CPUs for the large matrix multiplications that dominate Step 1
ridge regression and Step 2 linear/logistic model fitting.

**Requires:** MKL installed and discoverable by CMake's `find_package(MKL)`.
Cannot be combined with `WITH_OPENBLAS`.

---

### `WITH_OPENBLAS` — OpenBLAS

**Default:** OFF

Identical effect to `WITH_MKL` (sets `EIGEN_USE_BLAS` and
`EIGEN_USE_LAPACKE`) but links against
[OpenBLAS](https://www.openblas.net/), the portable open-source BLAS
implementation. Suitable for AMD, Arm, or environments where MKL is
not available.

**Requires:** OpenBLAS installed and discoverable by
`find_package(OpenBLAS)`. Cannot be combined with `WITH_MKL`.

---

### `WITH_BOOST_IOSTREAM` — Boost IOstreams

**Default:** OFF

Links [Boost.IOstreams](https://www.boost.org/doc/libs/release/libs/iostreams/)
and enables the `HAS_BOOST_IOSTREAM` compile definition used in
`Files.cpp`.

Without this flag, compressed input/output (`.gz` files) goes through
the raw zlib API. With it, **regenie** can use Boost's higher-level
filtering stream layer, which provides:

- Transparent gzip/bzip2 decompression of phenotype and covariate
  files without manual chunk management.
- Compressed association result output.

**Requires:** Boost with the `iostreams` component installed.

---

### `WITH_HTSLIB` — HTSlib / REMETA

**Default:** OFF

The most impactful optional dependency. Enabling it:

1. Links [HTSlib](https://www.htslib.org/) — the C library underlying
   samtools and bcftools — so **regenie** can read `.vcf.gz` and
   `.bcf` genotype files directly in addition to BGEN/PLINK formats.
2. Builds and links the bundled
   [REMETA](https://rgcgithub.github.io/remeta/) library
   (`external_libs/remeta/`) for rare-variant meta-analysis summary
   statistics.

Without this flag the `--anno-file` / REMETA meta-analysis code paths
in `Masks.cpp` and `Joint_Tests.cpp` are compiled out and VCF/BCF
input is not supported.

**Requires:** HTSlib installed and detectable via `pkg-config`.
Optionally, set `HTSLIB_PATH` to the installation prefix if HTSlib is
not on the default `PKG_CONFIG_PATH`.

```bash
cmake -DWITH_HTSLIB=ON -DHTSLIB_PATH=/opt/htslib ..
```
