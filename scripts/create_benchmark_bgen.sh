#!/usr/bin/env bash
# Creates a medium-sized BGEN benchmark dataset by concatenating
# example/example_3chr.bgen N times using the bgen library's cat-bgen tool.
#
# The result is a ~20 MB BGEN file with 100 000 variants across 500 samples,
# derived entirely from committed example data – no external downloads needed.
#
# Usage:
#   scripts/create_benchmark_bgen.sh \
#       --bgen-tools-dir <path-to-bgen-apps-dir> \
#       [--reps <N>]               (default: 200  → 100 000 variants)
#       [--outdir <dir>]           (default: benchmark_data/)
#
# Requirements: cat-bgen (and optionally bgenix) from the bgen library.
# Both are produced by `python3 waf` inside the bgen source tree at
#   <bgen-version>/build/apps/

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

BGEN_TOOLS_DIR=""
REPS=200
OUTDIR="${REPO_ROOT}/benchmark_data"

# ── Parse arguments ──────────────────────────────────────────────────────────
while [[ "$#" -gt 0 ]]; do
  case $1 in
    --bgen-tools-dir) BGEN_TOOLS_DIR="$2"; shift ;;
    --reps)           REPS="$2";           shift ;;
    --outdir)         OUTDIR="$2";         shift ;;
    -h|--help)
      sed -n '/^# Usage/,/^$/p' "$0"
      exit 0 ;;
    *) echo "ERROR: Unknown parameter: $1" >&2; exit 1 ;;
  esac
  shift
done

CAT_BGEN="${BGEN_TOOLS_DIR}/cat-bgen"
BGENIX="${BGEN_TOOLS_DIR}/bgenix"

if [[ ! -x "$CAT_BGEN" ]]; then
  echo "ERROR: cat-bgen not found or not executable at: $CAT_BGEN" >&2
  echo "  Pass --bgen-tools-dir pointing to the bgen library apps directory." >&2
  echo "  Build the bgen library first:  cd <bgen-version> && python3 waf configure && python3 waf" >&2
  exit 1
fi

SOURCE_BGEN="${REPO_ROOT}/example/example_3chr.bgen"
SOURCE_SAMPLE="${REPO_ROOT}/example/example_3chr.sample"

if [[ ! -f "$SOURCE_BGEN" ]]; then
  echo "ERROR: Source BGEN not found: $SOURCE_BGEN" >&2
  exit 1
fi

# ── Generate concatenated BGEN ────────────────────────────────────────────────
mkdir -p "$OUTDIR"
OUT_BGEN="${OUTDIR}/benchmark.bgen"
OUT_SAMPLE="${OUTDIR}/benchmark.sample"

SOURCE_NVARS=$(awk 'END{print NR}' "${REPO_ROOT}/example/example_3chr.bim")
TOTAL_VARS=$(( REPS * SOURCE_NVARS ))
echo "==> Generating benchmark BGEN"
echo "    Source    : ${SOURCE_BGEN}"
echo "    Source vars: ${SOURCE_NVARS} variants x ${REPS} reps = ${TOTAL_VARS} total variants"

# Build the -g argument list (repeat the source file REPS times)
INPUT_ARGS=()
for (( i = 0; i < REPS; i++ )); do
  INPUT_ARGS+=("${SOURCE_BGEN}")
done

"$CAT_BGEN" -g "${INPUT_ARGS[@]}" -og "${OUT_BGEN}"

FILE_SIZE=$(du -sh "${OUT_BGEN}" | cut -f1)
echo "    Created   : ${OUT_BGEN} (${FILE_SIZE})"

# ── Index for efficient access ────────────────────────────────────────────────
if [[ -x "$BGENIX" ]]; then
  echo "==> Indexing benchmark BGEN with bgenix..."
  "$BGENIX" -g "${OUT_BGEN}" -index
  echo "    Created   : ${OUT_BGEN}.bgi"
else
  echo "    (bgenix not found – skipping index creation)"
fi

# ── Copy sample file ──────────────────────────────────────────────────────────
cp "${SOURCE_SAMPLE}" "${OUT_SAMPLE}"

NSAMP=$(( $(awk 'END{print NR}' "${SOURCE_SAMPLE}") - 2 ))   # subtract 2 header lines

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "==> Benchmark dataset ready"
echo "    Variants : ${TOTAL_VARS}"
echo "    Samples  : ${NSAMP}"
echo "    Directory: ${OUTDIR}/"
echo ""
echo "    BGEN  : ${OUT_BGEN}"
echo "    Sample: ${OUT_SAMPLE}"
