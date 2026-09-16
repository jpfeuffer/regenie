#!/usr/bin/env bash
# ─── Regenie Benchmark Matrix ────────────────────────────────────────────────
#
# Builds and benchmarks regenie with multiple BLAS backends and
# runs C++ vs Fortran micro-benchmarks for QUADPACK and MVTDST.
#
# Usage:
#   test/benchmark_matrix.sh [OPTIONS]
#
# Options:
#   --variants <list>   comma-separated build variants to run
#                       (default: eigen,openblas)
#                       available: eigen,eigen-fortran,eigen-threadsafe-qf,openblas,mkl
#                       Note: eigen uses C QUADPACK; eigen-fortran uses Fortran QUADPACK;
#                             eigen-threadsafe-qf uses C QUADPACK + thread-safe qf_mt()
#   --bsize <N>         block size for regenie (default: 1000)
#   --reps <N>          repetitions for benchmark BGEN generation (default: 20)
#   --skip-build        skip cmake build steps (use existing builds)
#   --skip-regenie      skip full regenie benchmark (only micro-benchmarks)
#   --skip-micro        skip C++ vs Fortran micro-benchmarks
#   --run-vc            also run SKATO VC benchmark (auto-enabled if eigen-fortran in variants)
#   --n-vc-sets <N>     synthetic gene sets for VC benchmark (default: 500)
#   --jobs <N>          parallel build jobs (default: nproc)
#   --outdir <dir>      output directory for results (default: benchmark_results)
#   --bgen-dir <path>   bgen install prefix (auto-detected from existing build)
#   --prefix-path <p>   extra CMAKE_PREFIX_PATH entries
#
# Environment:
#   OPENBLAS_DIR        path to OpenBLAS installation (default: /opt/homebrew/opt/openblas)
#   MKL_DIR             path to MKL installation

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# ── Defaults ──────────────────────────────────────────────────────────────────
VARIANTS="eigen,openblas"
BSIZE=1000
REPS=20
SKIP_BUILD=0
SKIP_REGENIE=0
SKIP_MICRO=0
RUN_VC=0
N_VC_SETS=10
N_VC_MASKS=50
THREADS=1
JOBS=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)
OUTDIR="${REPO_ROOT}/benchmark_results"
BGEN_DIR=""
BGEN_TOOLS_DIR=""
EXTRA_PREFIX_PATH=""
OPENBLAS_DIR="${OPENBLAS_DIR:-/opt/homebrew/opt/openblas}"
MKL_DIR="${MKL_DIR:-}"

# ── Parse arguments ──────────────────────────────────────────────────────────
while [[ "$#" -gt 0 ]]; do
  case $1 in
    --variants)     VARIANTS="$2";        shift ;;
    --bsize)        BSIZE="$2";           shift ;;
    --reps)         REPS="$2";            shift ;;
    --skip-build)   SKIP_BUILD=1 ;;
    --skip-regenie) SKIP_REGENIE=1 ;;
    --skip-micro)   SKIP_MICRO=1 ;;
    --run-vc)       RUN_VC=1 ;;
    --n-vc-sets)    N_VC_SETS="$2";  shift ;;
    --n-vc-masks)   N_VC_MASKS="$2"; shift ;;
    --threads)      THREADS="$2";    shift ;;
    --jobs)         JOBS="$2";       shift ;;
    --outdir)       OUTDIR="$2";          shift ;;
    --bgen-dir)     BGEN_DIR="$2";        shift ;;
    --bgen-tools-dir) BGEN_TOOLS_DIR="$2"; shift ;;
    --prefix-path)  EXTRA_PREFIX_PATH="$2"; shift ;;
    -h|--help)
      sed -n '/^# Usage/,/^[^#]/{ /^[^#]/!p; }' "$0"
      exit 0 ;;
    *) echo "ERROR: Unknown parameter: $1" >&2; exit 1 ;;
  esac
  shift
done

# ── Auto-detect bgen_DIR from existing build ──────────────────────────────────
if [[ -z "$BGEN_DIR" ]]; then
  if [[ -f "${REPO_ROOT}/build/CMakeCache.txt" ]]; then
    BGEN_DIR=$(grep '^bgen_DIR' "${REPO_ROOT}/build/CMakeCache.txt" \
      | cut -d= -f2- | sed 's|/lib/cmake/bgen$||')
  fi
fi
if [[ -z "$BGEN_DIR" || ! -d "$BGEN_DIR" ]]; then
  echo "ERROR: Cannot find bgen install directory." >&2
  echo "  Pass --bgen-dir <path> or ensure build/CMakeCache.txt exists." >&2
  exit 1
fi

BGEN_CMAKE_DIR="${BGEN_DIR}/lib/cmake/bgen"

# Locate bgen tools (cat-bgen, bgenix) — prefer explicit --bgen-tools-dir,
# then test the install dir, then scan for a working copy.
BGEN_TOOLS_DIR="${BGEN_TOOLS_DIR:-${BGEN_DIR}/bin}"
if [[ ! -x "${BGEN_TOOLS_DIR}/cat-bgen" ]] || \
   ! "${BGEN_TOOLS_DIR}/cat-bgen" -help &>/dev/null; then
  # The install/bin copy may have broken RPATH; find one that actually works.
  for _candidate in \
    "$(dirname "${BGEN_DIR}")/build/cat-bgen" \
    "${HOME}/git/bgen/build/cat-bgen"; do
    if [[ -x "$_candidate" ]] && "$_candidate" -help &>/dev/null; then
      BGEN_TOOLS_DIR="$(dirname "$_candidate")"
      break
    fi
  done
  # Last resort: find any working cat-bgen under ~/git/bgen
  if ! "${BGEN_TOOLS_DIR}/cat-bgen" -help &>/dev/null 2>&1; then
    PIXI_BGEN=$(find "${HOME}/git/bgen" -name cat-bgen -path "*build*" -type f 2>/dev/null | head -1)
    if [[ -n "$PIXI_BGEN" ]] && "$PIXI_BGEN" -help &>/dev/null; then
      BGEN_TOOLS_DIR="$(dirname "$PIXI_BGEN")"
    else
      echo "WARNING: no working cat-bgen found; benchmark data generation will fail" >&2
    fi
  fi
fi

# Auto-enable VC benchmark when a comparison variant is requested
if [[ ("$VARIANTS" == *"eigen-fortran"* || "$VARIANTS" == *"eigen-threadsafe-qf"*) ]] && \
   [[ "$SKIP_REGENIE" -eq 0 ]]; then
  RUN_VC=1
fi

echo "╔══════════════════════════════════════════════════════════════════╗"
echo "║           REGENIE BENCHMARK MATRIX                             ║"
echo "╚══════════════════════════════════════════════════════════════════╝"
echo ""
echo "  Variants    : ${VARIANTS}"
echo "  bgen dir    : ${BGEN_DIR}"
echo "  OpenBLAS    : ${OPENBLAS_DIR}"
echo "  MKL         : ${MKL_DIR:-not set}"
echo "  Jobs        : ${JOBS}"
echo "  bsize       : ${BSIZE}"
echo "  Output      : ${OUTDIR}"
if [[ "$RUN_VC" -eq 1 ]]; then
  echo "  VC sets     : ${N_VC_SETS} sets × ${N_VC_MASKS} masks (SKATO benchmark)"
fi
if [[ "$THREADS" -gt 1 ]]; then
  echo "  Threads     : ${THREADS}"
fi
echo ""

mkdir -p "${OUTDIR}"

# ── Common CMake args ─────────────────────────────────────────────────────────
BASE_CMAKE_ARGS=(
  -GNinja
  -DCMAKE_BUILD_TYPE=Release
  -Dbgen_DIR="${BGEN_CMAKE_DIR}"
  -DWITH_AWS_S3=OFF
  -DBUILD_STATIC=OFF
  -DCMAKE_POLICY_VERSION_MINIMUM=3.5
  -Wno-dev
)

# Help Apple Clang find libomp (Homebrew)
if [[ "$(uname)" == "Darwin" && -d /opt/homebrew/opt/libomp ]]; then
  BASE_CMAKE_ARGS+=(
    -DOpenMP_C_FLAGS="-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include"
    -DOpenMP_CXX_FLAGS="-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include"
    -DOpenMP_C_LIB_NAMES="omp"
    -DOpenMP_CXX_LIB_NAMES="omp"
    -DOpenMP_omp_LIBRARY="/opt/homebrew/opt/libomp/lib/libomp.dylib"
  )
fi
if [[ -n "$EXTRA_PREFIX_PATH" ]]; then
  BASE_CMAKE_ARGS+=(-DCMAKE_PREFIX_PATH="${EXTRA_PREFIX_PATH}")
fi

# ── Build function ────────────────────────────────────────────────────────────
build_variant() {
  local name="$1"; shift
  local build_dir="${REPO_ROOT}/build_bench_${name}"
  local extra_args=("$@")

  if [[ "$SKIP_BUILD" -eq 1 && -x "${build_dir}/bin/regenie" ]]; then
    echo "  [${name}] Using existing build at ${build_dir}"
    return 0
  fi

  echo "  [${name}] Configuring..."
  cmake -B "${build_dir}" "${REPO_ROOT}" \
    "${BASE_CMAKE_ARGS[@]}" "${extra_args[@]}" \
    > "${OUTDIR}/build_${name}_configure.log" 2>&1

  echo "  [${name}] Building (${JOBS} jobs)..."
  cmake --build "${build_dir}" -j "${JOBS}" \
    > "${OUTDIR}/build_${name}_build.log" 2>&1

  if [[ ! -x "${build_dir}/bin/regenie" ]]; then
    echo "  [${name}] ERROR: build failed — see ${OUTDIR}/build_${name}_build.log" >&2
    return 1
  fi
  echo "  [${name}] Build OK: ${build_dir}/bin/regenie"
}

# ── Generate benchmark data (once) ───────────────────────────────────────────
BENCH_DATADIR="${OUTDIR}/benchmark_data"
if [[ ! -f "${BENCH_DATADIR}/benchmark.bgen" ]]; then
  echo "==> Generating benchmark BGEN (${REPS} reps)..."
  bash "${REPO_ROOT}/scripts/create_benchmark_bgen.sh" \
    --bgen-tools-dir "${BGEN_TOOLS_DIR}" \
    --reps "${REPS}" \
    --outdir "${BENCH_DATADIR}"
else
  echo "==> Benchmark BGEN already exists at ${BENCH_DATADIR}"
fi
echo ""

# ── Build all variants ────────────────────────────────────────────────────────
echo "==> Building variants..."

IFS=',' read -ra VARIANT_LIST <<< "$VARIANTS"

for variant in "${VARIANT_LIST[@]}"; do
  case "$variant" in
    eigen)
      build_variant eigen \
        -DWITH_MKL=OFF \
        -DWITH_OPENBLAS=OFF \
        -DWITH_C_QUADPACK=ON \
        -DMVTDST_BUILD_BENCHMARKS=ON
      ;;
    eigen-fortran)
      # Same as eigen but uses Fortran QUADPACK instead of C reimplementation
      build_variant eigen-fortran \
        -DWITH_MKL=OFF \
        -DWITH_OPENBLAS=OFF \
        -DWITH_C_QUADPACK=OFF
      ;;
    eigen-threadsafe-qf)
      # Same as eigen but uses thread-safe qf_mt() (future: enables parallel mask loop)
      build_variant eigen-threadsafe-qf \
        -DWITH_MKL=OFF \
        -DWITH_OPENBLAS=OFF \
        -DWITH_C_QUADPACK=ON \
        -DWITH_THREADSAFE_QF=ON \
        -DMVTDST_BUILD_BENCHMARKS=ON
      ;;
    openblas)
      if [[ ! -d "$OPENBLAS_DIR" ]]; then
        echo "  [openblas] SKIPPED: ${OPENBLAS_DIR} not found" >&2
        continue
      fi
      # FindBLAS locates OpenBLAS via BLA_VENDOR=OpenBLAS (set in CMakeLists.txt
      # when WITH_OPENBLAS=ON). Point CMAKE_PREFIX_PATH at the OpenBLAS prefix so
      # FindBLAS can find the library.
      build_variant openblas \
        -DWITH_MKL=OFF \
        -DWITH_OPENBLAS=ON \
        -DWITH_C_QUADPACK=ON \
        -DCMAKE_PREFIX_PATH="${OPENBLAS_DIR}${EXTRA_PREFIX_PATH:+;${EXTRA_PREFIX_PATH}}"
      ;;
    mkl)
      if [[ -z "$MKL_DIR" || ! -d "$MKL_DIR" ]]; then
        echo "  [mkl] SKIPPED: MKL_DIR not set or not found" >&2
        continue
      fi
      # FindBLAS locates MKL via BLA_VENDOR=Intel10_64lp_dyn (set in CMakeLists.txt
      # when WITH_MKL=ON). Add the MKL root to CMAKE_PREFIX_PATH.
      build_variant mkl \
        -DWITH_MKL=ON \
        -DWITH_OPENBLAS=OFF \
        -DWITH_C_QUADPACK=ON \
        -DCMAKE_PREFIX_PATH="${MKL_DIR}${EXTRA_PREFIX_PATH:+;${EXTRA_PREFIX_PATH}}"
      ;;
    *)
      echo "  WARNING: Unknown variant '${variant}' — skipping" >&2
      ;;
  esac
done
echo ""

# ── Run regenie benchmarks ────────────────────────────────────────────────────
if [[ "$SKIP_REGENIE" -eq 0 ]]; then
  echo "==> Running regenie benchmarks..."
  echo ""

  SUMMARY_FILE="${OUTDIR}/benchmark_summary.txt"
  printf "%-15s  %10s  %10s  %10s  %12s  %10s  %12s\n" \
    "variant" "step1_s" "step2_s" "variants" "throughput" "vc_s" "step1+2_s" \
    > "$SUMMARY_FILE"
  printf "%s\n" "$(printf '─%.0s' {1..90})" >> "$SUMMARY_FILE"

  for variant in "${VARIANT_LIST[@]}"; do
    BUILD_DIR="${REPO_ROOT}/build_bench_${variant}"
    REGENIE_BIN="${BUILD_DIR}/bin/regenie"

    if [[ ! -x "$REGENIE_BIN" ]]; then
      echo "  [${variant}] SKIPPED: binary not found"
      continue
    fi

    echo "  [${variant}] Running step 1 + step 2..."
    VARIANT_DATADIR="${OUTDIR}/results_${variant}"
    mkdir -p "$VARIANT_DATADIR"

    # Copy benchmark data paths
    BENCH_EXTRA_ARGS=()
    [[ "$RUN_VC" -eq 1 ]] && BENCH_EXTRA_ARGS+=(--run-vc --n-vc-sets "$N_VC_SETS" --n-vc-masks "$N_VC_MASKS")
    [[ "$THREADS" -gt 1 ]] && BENCH_EXTRA_ARGS+=(--threads "$THREADS")

    bash "${REPO_ROOT}/test/benchmark.sh" \
      --regenie "$REGENIE_BIN" \
      --datadir "${BENCH_DATADIR}" \
      --bsize "$BSIZE" \
      "${BENCH_EXTRA_ARGS[@]}" \
      > "${VARIANT_DATADIR}/benchmark_output.log" 2>&1 || true

    # Copy timing results
    if [[ -f "${BENCH_DATADIR}/benchmark_timing.txt" ]]; then
      cp "${BENCH_DATADIR}/benchmark_timing.txt" "${VARIANT_DATADIR}/benchmark_timing.txt"

      # Parse results for summary
      S1_ELAPSED=$(grep '^step1_elapsed_s=' "${VARIANT_DATADIR}/benchmark_timing.txt" 2>/dev/null | cut -d= -f2 || echo "N/A")
      S2_ELAPSED=$(grep '^step2_elapsed_s=' "${VARIANT_DATADIR}/benchmark_timing.txt" 2>/dev/null | cut -d= -f2 || echo "N/A")
      VC_ELAPSED=$(grep '^vc_elapsed_s='    "${VARIANT_DATADIR}/benchmark_timing.txt" 2>/dev/null | cut -d= -f2 || echo "")
      NVARS=$(grep '^variants=' "${VARIANT_DATADIR}/benchmark_timing.txt" 2>/dev/null | cut -d= -f2 || echo "0")
      THROUGHPUT=$(grep '^throughput_vars_per_sec=' "${VARIANT_DATADIR}/benchmark_timing.txt" 2>/dev/null | cut -d= -f2 || echo "0")

      # Compute total
      TOTAL="N/A"
      if [[ "$S1_ELAPSED" != "N/A" && "$S2_ELAPSED" != "N/A" ]] && command -v bc &>/dev/null; then
        TOTAL=$(echo "$S1_ELAPSED + $S2_ELAPSED" | bc 2>/dev/null || echo "N/A")
      fi

      VC_DISP="${VC_ELAPSED:-N/A}"
      printf "%-15s  %10s  %10s  %10s  %12s  %10s  %12s\n" \
        "$variant" "$S1_ELAPSED" "$S2_ELAPSED" "$NVARS" "$THROUGHPUT" "$VC_DISP" "$TOTAL" \
        >> "$SUMMARY_FILE"

      echo "  [${variant}] step1=${S1_ELAPSED}s  step2=${S2_ELAPSED}s  throughput=${THROUGHPUT} v/s${VC_ELAPSED:+  vc(SKATO)=${VC_ELAPSED}s}"
    else
      echo "  [${variant}] WARNING: no timing results produced"
      printf "%-15s  %10s  %10s  %10s  %12s  %10s  %12s\n" \
        "$variant" "FAILED" "FAILED" "-" "-" "-" "-" \
        >> "$SUMMARY_FILE"
    fi
  done

  echo ""
  echo "==> Regenie Benchmark Summary"
  echo ""
  cat "$SUMMARY_FILE"
  echo ""

  # ── C (cquadpack) vs Fortran QUADPACK SKATO comparison ────────────────────
  if [[ "$RUN_VC" -eq 1 ]] && \
     [[ -f "${OUTDIR}/results_eigen/benchmark_timing.txt" ]] && \
     [[ -f "${OUTDIR}/results_eigen-fortran/benchmark_timing.txt" ]]; then
    VC_CPP=$(grep '^vc_elapsed_s=' "${OUTDIR}/results_eigen/benchmark_timing.txt" 2>/dev/null | cut -d= -f2 || echo "")
    VC_FORT=$(grep '^vc_elapsed_s=' "${OUTDIR}/results_eigen-fortran/benchmark_timing.txt" 2>/dev/null | cut -d= -f2 || echo "")
    if [[ -n "$VC_CPP" && -n "$VC_FORT" ]] && command -v bc &>/dev/null; then
      SPEEDUP=$(echo "scale=3; $VC_FORT / $VC_CPP" | bc 2>/dev/null || echo "N/A")
      echo "==> C reimplementation vs Fortran — SKATO (QUADPACK) end-to-end"
      echo ""
      printf "  %-25s : %s s\n" "eigen (C cquadpack)"    "$VC_CPP"
      printf "  %-25s : %s s\n" "eigen-fortran (Fortran)" "$VC_FORT"
      printf "  %-25s : %sx\n"  "Speedup (Fortran/C)"     "$SPEEDUP"
      echo "  Sets: ${N_VC_SETS} sets × ${N_VC_MASKS} masks  ×  8 rho  ×  2 phenotypes = $(( N_VC_SETS * N_VC_MASKS * 8 * 2 )) QUADPACK calls"
      echo ""
    fi
  fi

  # ── thread-safe qf_mt() vs original qf() comparison ─────────────────────
  if [[ "$RUN_VC" -eq 1 ]] && \
     [[ -f "${OUTDIR}/results_eigen/benchmark_timing.txt" ]] && \
     [[ -f "${OUTDIR}/results_eigen-threadsafe-qf/benchmark_timing.txt" ]]; then
    VC_QF=$(grep '^vc_elapsed_s=' "${OUTDIR}/results_eigen/benchmark_timing.txt" 2>/dev/null | cut -d= -f2 || echo "")
    VC_QFMT=$(grep '^vc_elapsed_s=' "${OUTDIR}/results_eigen-threadsafe-qf/benchmark_timing.txt" 2>/dev/null | cut -d= -f2 || echo "")
    if [[ -n "$VC_QF" && -n "$VC_QFMT" ]] && command -v bc &>/dev/null; then
      SPEEDUP_MT=$(echo "scale=3; $VC_QF / $VC_QFMT" | bc 2>/dev/null || echo "N/A")
      echo "==> Original qf() vs thread-safe qf_mt() — SKATO end-to-end"
      echo ""
      printf "  %-30s : %s s\n" "eigen (global-state qf, 1T)"   "$VC_QF"
      printf "  %-30s : %s s\n" "eigen-threadsafe-qf (qf_mt)"   "$VC_QFMT"
      printf "  %-30s : %sx\n"  "Speedup"                        "$SPEEDUP_MT"
      if [[ "$THREADS" -gt 1 ]]; then
        echo "  (qf_mt running with ${THREADS} OpenMP threads)"
      else
        echo "  (ratio ~1 expected at 1 thread — struct overhead only)"
      fi
      echo ""
    fi
  fi
fi

# ── C++ vs Fortran micro-benchmarks ──────────────────────────────────────────
if [[ "$SKIP_MICRO" -eq 0 ]]; then
  echo "==> C++ vs Fortran micro-benchmarks"
  echo ""

  # QUADPACK benchmark (from the 'eigen' build which has WITH_C_QUADPACK=ON)
  QP_BUILD="${REPO_ROOT}/build_bench_eigen"
  QP_BIN="${QP_BUILD}/bin/bench_quadpack"
  if [[ -x "$QP_BIN" ]]; then
    echo "  ── QUADPACK: C (cquadpack) vs Fortran (libquad.a) ──"
    "$QP_BIN" | tee "${OUTDIR}/bench_quadpack_output.txt"
    echo ""
  else
    echo "  [quadpack] Building bench_quadpack target..."
    if cmake --build "$QP_BUILD" --target bench_quadpack -j "$JOBS" \
         > "${OUTDIR}/build_bench_quadpack.log" 2>&1; then
      echo "  ── QUADPACK: C (cquadpack) vs Fortran (libquad.a) ──"
      "${QP_BUILD}/bin/bench_quadpack" | tee "${OUTDIR}/bench_quadpack_output.txt"
      echo ""
    else
      echo "  [quadpack] Build failed — see ${OUTDIR}/build_bench_quadpack.log" >&2
    fi
  fi

  # MVTDST benchmark (from the 'eigen' build which has MVTDST_BUILD_BENCHMARKS=ON)
  MVTDST_BIN="${QP_BUILD}/bin/benchmark_mvtdst"
  if [[ -x "$MVTDST_BIN" ]]; then
    echo "  ── MVTDST: C++ port vs Fortran original ──"
    "$MVTDST_BIN" | tee "${OUTDIR}/bench_mvtdst_output.txt"
    echo ""
  else
    echo "  [mvtdst] Building benchmark_mvtdst target..."
    if cmake --build "$QP_BUILD" --target benchmark_mvtdst -j "$JOBS" \
         > "${OUTDIR}/build_bench_mvtdst.log" 2>&1; then
      echo "  ── MVTDST: C++ port vs Fortran original ──"
      "${QP_BUILD}/extern/mvtdst/benchmark_mvtdst" | tee "${OUTDIR}/bench_mvtdst_output.txt"
      echo ""
    else
      echo "  [mvtdst] Build failed — see ${OUTDIR}/build_bench_mvtdst.log" >&2
    fi
  fi
fi

# ── Final summary ─────────────────────────────────────────────────────────────
echo "╔══════════════════════════════════════════════════════════════════╗"
echo "║           BENCHMARK COMPLETE                                   ║"
echo "╚══════════════════════════════════════════════════════════════════╝"
echo ""
echo "  Results in: ${OUTDIR}/"
echo ""
ls -la "${OUTDIR}"/*.txt 2>/dev/null | awk '{print "    " $NF}'
echo ""
