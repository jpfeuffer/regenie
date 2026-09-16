#!/usr/bin/env bash
# Runs both steps of regenie on a pre-generated benchmark BGEN file and
# records per-step timing, throughput, and peak memory usage.
#
#   Step 1: ridge-regression training on the benchmark BGEN → LOCO predictions
#   Step 2: GWAS scan using those LOCO predictions
#
# Prerequisite: generate the benchmark data first with
#   scripts/create_benchmark_bgen.sh --bgen-tools-dir <dir>
#
# Usage:
#   test/benchmark.sh [OPTIONS]
#
# Options:
#   --regenie    <path>  path to regenie binary (auto-detected if omitted)
#   --datadir    <dir>   directory containing benchmark.bgen/.sample
#                          (default: <repo-root>/benchmark_data)
#   --path       <dir>   repo root (default: parent of this script's directory)
#   --bsize      <N>     BGEN read-block size (default: 1000)
#   --skip-step1         skip step 1; run step 2 with --ignore-pred instead
#   --run-vc             also run SKATO variance-component benchmark (exercises QUADPACK)
#   --n-vc-sets  <N>     number of synthetic gene sets for VC benchmark (default: 10)
#   --n-vc-masks <N>     number of mask definitions per gene set (default: 50)
#                          total QUADPACK calls = n-vc-sets × n-vc-masks × 8 rho × 2 pheno
#   --threads    <N>     number of threads passed to regenie --threads (default: 1)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

REGENIE_BIN=""
DATADIR="${REPO_ROOT}/benchmark_data"
BSIZE=1000
SKIP_STEP1=0
RUN_VC=0
N_VC_SETS=10
N_VC_MASKS=50
THREADS=1

# ── Parse arguments ──────────────────────────────────────────────────────────
while [[ "$#" -gt 0 ]]; do
  case $1 in
    --regenie)    REGENIE_BIN="$2"; shift ;;
    --datadir)    DATADIR="$2";     shift ;;
    --path)       REPO_ROOT="$2";   shift ;;
    --bsize)      BSIZE="$2";       shift ;;
    --skip-step1) SKIP_STEP1=1 ;;
    --run-vc)      RUN_VC=1 ;;
    --n-vc-sets)   N_VC_SETS="$2";  shift ;;
    --n-vc-masks)  N_VC_MASKS="$2"; shift ;;
    --threads)     THREADS="$2";    shift ;;
    -h|--help)
      sed -n '/^# Usage/,/^$/p' "$0"
      exit 0 ;;
    *) echo "ERROR: Unknown parameter: $1" >&2; exit 1 ;;
  esac
  shift
done

# ── Locate the regenie binary ────────────────────────────────────────────────
if [[ -z "$REGENIE_BIN" ]]; then
  for candidate in \
    "${REPO_ROOT}/build/bin/regenie" \
    "${REPO_ROOT}/bin/regenie" \
    "${REPO_ROOT}/regenie"; do
    if [[ -x "$candidate" ]]; then
      REGENIE_BIN="$candidate"
      break
    fi
  done
fi

if [[ -z "$REGENIE_BIN" || ! -x "$REGENIE_BIN" ]]; then
  echo "ERROR: regenie binary not found." >&2
  echo "  Pass --regenie <path> or build the project first." >&2
  exit 1
fi

# ── Input paths ───────────────────────────────────────────────────────────────
BENCH_BGEN="${DATADIR}/benchmark.bgen"
BENCH_SAMPLE="${DATADIR}/benchmark.sample"
COVAR_FILE="${REPO_ROOT}/example/covariates.txt"
PHENO_FILE="${REPO_ROOT}/example/phenotype.txt"
RESULTS_FILE="${DATADIR}/benchmark_timing.txt"

for f in "$BENCH_BGEN" "$BENCH_SAMPLE" "$COVAR_FILE" "$PHENO_FILE"; do
  if [[ ! -f "$f" ]]; then
    echo "ERROR: Required file not found: $f" >&2
    exit 1
  fi
done

# ── Time command ──────────────────────────────────────────────────────────────
# Linux: /usr/bin/time -v gives "Maximum resident set size"
# macOS with GNU coreutils: gtime -v
TIME_CMD=""
if /usr/bin/time -v true 2>&1 | grep -q "Maximum resident"; then
  TIME_CMD="/usr/bin/time -v"
elif command -v gtime &>/dev/null && gtime -v true 2>&1 | grep -q "Maximum resident"; then
  TIME_CMD="gtime -v"
fi

TMPDIR_BENCH=$(mktemp -d)
trap 'rm -rf "${TMPDIR_BENCH}"' EXIT

REGENIE_VERSION=$("$REGENIE_BIN" --version 2>&1 | head -1 || true)

echo "==> regenie full-pipeline BGEN benchmark"
echo "    Binary  : ${REGENIE_BIN}"
echo "    Version : ${REGENIE_VERSION}"
echo "    Input   : ${BENCH_BGEN}"
echo "    bsize   : ${BSIZE}"
[[ "$SKIP_STEP1" -eq 1 ]] && echo "    Mode    : step 2 only (--ignore-pred)"
[[ "$RUN_VC"    -eq 1 ]] && echo "    VC sets : ${N_VC_SETS} synthetic gene sets (SKATO)"
[[ "$THREADS" -gt 1 ]] && echo "    Threads : ${THREADS}"
echo ""

# ── Helper: run one regenie step, emit key=value pairs ────────────────────────
# Prints: <prefix>__wall=N  <prefix>__elapsed_s=X  <prefix>__peak_rss_kb=N
run_timed() {
  local prefix="$1"; shift   # e.g. "step1"
  local log_file="$1"; shift # regenie's own .log path
  local time_log="$1"; shift # /usr/bin/time verbose output path
  # remaining args: the regenie command

  local start end wall elapsed_s peak_rss_kb
  start=$(date +%s)

  if [[ -n "$TIME_CMD" ]]; then
    $TIME_CMD "$@" 2>&1 | tee "${time_log}"
  else
    "$@"
  fi

  end=$(date +%s)
  wall=$(( end - start ))

  elapsed_s=$(grep '^Elapsed time' "${log_file}" 2>/dev/null \
    | grep -oE '[0-9]+[.][0-9]+' | head -1 || echo "")

  peak_rss_kb=""
  if [[ -f "$time_log" ]]; then
    peak_rss_kb=$(grep -oP 'Maximum resident set size \(kbytes\): \K[0-9]+' \
      "$time_log" 2>/dev/null || true)
  fi

  echo "${prefix}__wall=${wall}"
  echo "${prefix}__elapsed_s=${elapsed_s}"
  echo "${prefix}__peak_rss_kb=${peak_rss_kb}"
}

# ── Step 1: ridge-regression training ─────────────────────────────────────────
STEP1_WALL=0; STEP1_ELAPSED_S=""; STEP1_PEAK_RSS_KB=""

if [[ "$SKIP_STEP1" -eq 0 ]]; then
  echo "==> Step 1: ridge-regression training"
  while IFS='=' read -r k v; do
    case $k in
      step1__wall)        STEP1_WALL=$v ;;
      step1__elapsed_s)   STEP1_ELAPSED_S=$v ;;
      step1__peak_rss_kb) STEP1_PEAK_RSS_KB=$v ;;
    esac
  done < <(run_timed step1 \
    "${TMPDIR_BENCH}/bench_step1_out.log" \
    "${DATADIR}/benchmark_time_step1_verbose.txt" \
    "$REGENIE_BIN" \
      --step 1 \
      --bgen    "${BENCH_BGEN}" \
      --sample  "${BENCH_SAMPLE}" \
      --covarFile "${COVAR_FILE}" \
      --phenoFile "${PHENO_FILE}" \
      --bsize   "${BSIZE}" \
      --threads "${THREADS}" \
      --lowmem \
      --lowmem-prefix "${TMPDIR_BENCH}/bench_s1_lowmem" \
      --out     "${TMPDIR_BENCH}/bench_step1_out")

  S1_DISP="${STEP1_WALL}s"; [[ -n "$STEP1_ELAPSED_S" ]] && S1_DISP="${STEP1_ELAPSED_S}s"
  echo ""
  echo "    Step 1 elapsed : ${S1_DISP}"
  [[ -n "$STEP1_PEAK_RSS_KB" ]] && echo "    Step 1 peak RSS: $(( STEP1_PEAK_RSS_KB / 1024 )) MB"
  echo ""
fi

# ── Step 2: GWAS scan ─────────────────────────────────────────────────────────
echo "==> Step 2: GWAS scan"

STEP2_EXTRA_ARGS=()
if [[ "$SKIP_STEP1" -eq 0 && -f "${TMPDIR_BENCH}/bench_step1_out_pred.list" ]]; then
  STEP2_EXTRA_ARGS+=(--pred "${TMPDIR_BENCH}/bench_step1_out_pred.list")
else
  STEP2_EXTRA_ARGS+=(--ignore-pred)
fi

STEP2_WALL=0; STEP2_ELAPSED_S=""; STEP2_PEAK_RSS_KB=""

while IFS='=' read -r k v; do
  case $k in
    step2__wall)        STEP2_WALL=$v ;;
    step2__elapsed_s)   STEP2_ELAPSED_S=$v ;;
    step2__peak_rss_kb) STEP2_PEAK_RSS_KB=$v ;;
  esac
done < <(run_timed step2 \
  "${TMPDIR_BENCH}/bench_step2_out.log" \
  "${DATADIR}/benchmark_time_step2_verbose.txt" \
  "$REGENIE_BIN" \
    --step 2 \
    --bgen    "${BENCH_BGEN}" \
    --sample  "${BENCH_SAMPLE}" \
    --covarFile "${COVAR_FILE}" \
    --phenoFile "${PHENO_FILE}" \
    --bsize   "${BSIZE}" \
    --threads "${THREADS}" \
    --force-qt \
    "${STEP2_EXTRA_ARGS[@]}" \
    --out     "${TMPDIR_BENCH}/bench_step2_out")

# Count variants tested (header excluded)
NVARS=0
if [[ -f "${TMPDIR_BENCH}/bench_step2_out_Y1.regenie" ]]; then
  NVARS=$(( $(wc -l < "${TMPDIR_BENCH}/bench_step2_out_Y1.regenie") - 1 ))
fi

S2_DISP="${STEP2_WALL}s"; [[ -n "$STEP2_ELAPSED_S" ]] && S2_DISP="${STEP2_ELAPSED_S}s"

THROUGHPUT=0
if [[ -n "$STEP2_ELAPSED_S" ]] && command -v bc &>/dev/null && (( NVARS > 0 )); then
  THROUGHPUT=$(echo "scale=0; $NVARS / $STEP2_ELAPSED_S" | bc 2>/dev/null || echo 0)
elif (( STEP2_WALL > 0 && NVARS > 0 )); then
  THROUGHPUT=$(( NVARS / STEP2_WALL ))
fi

# ── Step 3: SKATO variance-component benchmark (exercises QUADPACK) ───────────
VC_WALL=0; VC_ELAPSED_S=""; VC_PEAK_RSS_KB=""

if [[ "$RUN_VC" -eq 1 ]]; then
  EXAMPLE_BGEN="${REPO_ROOT}/example/example_3chr.bgen"
  EXAMPLE_SAMPLE="${REPO_ROOT}/example/example_3chr.sample"

  if [[ ! -f "$EXAMPLE_BGEN" || ! -f "$EXAMPLE_SAMPLE" ]]; then
    echo "WARNING: ${EXAMPLE_BGEN} not found — skipping VC benchmark" >&2
    RUN_VC=0
  else
    echo "==> Step 3: SKATO VC benchmark (${N_VC_SETS} gene sets × ${N_VC_MASKS} masks)"
    echo "    Total QUADPACK calls: $(( N_VC_SETS * N_VC_MASKS * 8 * 2 )) (sets × masks × 8 rho × 2 pheno)"
    echo ""

    # ── Generate synthetic gene-set files ────────────────────────────────────
    VC_SETLIST="${TMPDIR_BENCH}/bench_vc.setlist"
    VC_ANNO="${TMPDIR_BENCH}/bench_vc.annotations"
    VC_MASKS_FILE="${TMPDIR_BENCH}/bench_vc.masks"

    # setlist: N_VC_SETS sets on chr1, each with mog_0..mog_39
    awk -v n="${N_VC_SETS}" 'BEGIN {
      vlist = ""; for (j=0; j<40; j++) vlist = vlist (j?",":"") "mog_" j
      for (i=1; i<=n; i++) printf "SYNTHSET%d 1 %d %s\n", i, i*1000, vlist
    }' > "$VC_SETLIST"

    # annotations: one category "func" per (variant, set) pair — regenie
    # allows only one category per variant per gene set; to get N_VC_MASKS
    # parallel masks we define N_VC_MASKS mask lines all referencing "func"
    awk -v n="${N_VC_SETS}" 'BEGIN {
      for (i=1; i<=n; i++)
        for (j=0; j<40; j++) printf "mog_%d SYNTHSET%d func\n", j, i
    }' > "$VC_ANNO"

    # N_VC_MASKS mask definitions, all using "func" — identical masks that
    # each trigger a full SKATO computation (saturates the OMP parallel for)
    awk -v m="${N_VC_MASKS}" 'BEGIN {
      for (k=1; k<=m; k++) printf "M%d func\n", k
    }' > "$VC_MASKS_FILE"

    # ── Run SKATO step ────────────────────────────────────────────────────────
    while IFS='=' read -r k v; do
      case $k in
        vc__wall)        VC_WALL=$v ;;
        vc__elapsed_s)   VC_ELAPSED_S=$v ;;
        vc__peak_rss_kb) VC_PEAK_RSS_KB=$v ;;
      esac
    done < <(run_timed vc \
      "${TMPDIR_BENCH}/bench_vc_out.log" \
      "${TMPDIR_BENCH}/benchmark_time_vc_verbose.txt" \
      "$REGENIE_BIN" \
        --step 2 \
        --bgen      "${EXAMPLE_BGEN}" \
        --sample    "${EXAMPLE_SAMPLE}" \
        --covarFile "${COVAR_FILE}" \
        --phenoFile "${PHENO_FILE}" \
        --bsize     200 \
        --ignore-pred \
        --force-qt \
        --set-list  "${VC_SETLIST}" \
        --anno-file "${VC_ANNO}" \
        --mask-def  "${VC_MASKS_FILE}" \
        --vc-tests  skato \
        --aaf-bins  0.2 \
        --threads   "${THREADS}" \
        --out       "${TMPDIR_BENCH}/bench_vc_out")

    VC_DISP="${VC_WALL}s"; [[ -n "$VC_ELAPSED_S" ]] && VC_DISP="${VC_ELAPSED_S}s"
    echo ""
    echo "    VC elapsed      : ${VC_DISP}"
    [[ -n "$VC_PEAK_RSS_KB" ]] && echo "    VC peak RSS     : $(( VC_PEAK_RSS_KB / 1024 )) MB"
    echo ""
  fi
fi

# ── Report ────────────────────────────────────────────────────────────────────
echo ""
echo "==> Benchmark results"
if [[ "$SKIP_STEP1" -eq 0 ]]; then
  S1_DISP="${STEP1_WALL}s"; [[ -n "$STEP1_ELAPSED_S" ]] && S1_DISP="${STEP1_ELAPSED_S}s"
  echo "    Step 1 elapsed  : ${S1_DISP}"
  [[ -n "$STEP1_PEAK_RSS_KB" ]] && echo "    Step 1 peak RSS : $(( STEP1_PEAK_RSS_KB / 1024 )) MB"
fi
echo "    Step 2 elapsed  : ${S2_DISP}"
echo "    Variants tested : ${NVARS}"
echo "    Throughput      : ~${THROUGHPUT} variants/sec"
[[ -n "$STEP2_PEAK_RSS_KB" ]] && echo "    Step 2 peak RSS : $(( STEP2_PEAK_RSS_KB / 1024 )) MB"
if [[ "$RUN_VC" -eq 1 ]]; then
  VC_DISP="${VC_WALL}s"; [[ -n "$VC_ELAPSED_S" ]] && VC_DISP="${VC_ELAPSED_S}s"
  echo "    VC (SKATO) elapsed : ${VC_DISP}  [${N_VC_SETS} sets × 8 rho × 2 pheno = $(( N_VC_SETS * 8 * 2 )) QUADPACK calls]"
  [[ -n "$VC_PEAK_RSS_KB" ]] && echo "    VC peak RSS     : $(( VC_PEAK_RSS_KB / 1024 )) MB"
fi

# ── Write machine-readable results ───────────────────────────────────────────
mkdir -p "$(dirname "$RESULTS_FILE")"
{
  echo "timestamp=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "git_sha=$(git -C "${REPO_ROOT}" rev-parse --short HEAD 2>/dev/null || echo unknown)"
  echo "regenie_version=${REGENIE_VERSION}"
  if [[ "$SKIP_STEP1" -eq 0 ]]; then
    echo "step1_wall_s=${STEP1_WALL}"
    [[ -n "$STEP1_ELAPSED_S" ]]   && echo "step1_elapsed_s=${STEP1_ELAPSED_S}"
    [[ -n "$STEP1_PEAK_RSS_KB" ]] && echo "step1_peak_rss_kb=${STEP1_PEAK_RSS_KB}"
  fi
  echo "step2_wall_s=${STEP2_WALL}"
  [[ -n "$STEP2_ELAPSED_S" ]]   && echo "step2_elapsed_s=${STEP2_ELAPSED_S}"
  [[ -n "$STEP2_PEAK_RSS_KB" ]] && echo "step2_peak_rss_kb=${STEP2_PEAK_RSS_KB}"
  echo "variants=${NVARS}"
  echo "throughput_vars_per_sec=${THROUGHPUT}"
  echo "bsize=${BSIZE}"
  echo "bgen_file=${BENCH_BGEN}"
  if [[ "$RUN_VC" -eq 1 ]]; then
    echo "vc_n_sets=${N_VC_SETS}"
    echo "vc_wall_s=${VC_WALL}"
    [[ -n "$VC_ELAPSED_S" ]]   && echo "vc_elapsed_s=${VC_ELAPSED_S}"
    [[ -n "$VC_PEAK_RSS_KB" ]] && echo "vc_peak_rss_kb=${VC_PEAK_RSS_KB}"
  fi
} > "$RESULTS_FILE"

echo ""
echo "    Results saved to: ${RESULTS_FILE}"
