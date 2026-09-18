#!/usr/bin/env bash
# Times regenie step 2 over the large local BGEN, repeated, reporting per-run
# and best/median throughput.
#
# Exists to guard the local read path against regressions: it is the hot loop
# for every run, so a change there needs a before/after number, not a guess.
#
# Usage: scripts/bench_step2.sh <regenie-binary> <label> [reps] [threads]

set -euo pipefail

BIN="${1:?usage: bench_step2.sh <regenie-binary> <label> [reps] [threads]}"
LABEL="${2:?missing label}"
REPS="${3:-5}"
THREADS="${4:-4}"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BGEN="${REPO_ROOT}/benchmark_data_large/benchmark.bgen"
SAMPLE="${REPO_ROOT}/benchmark_data_large/benchmark.sample"
PHENO="${REPO_ROOT}/example/phenotype_bin.txt"
COVAR="${REPO_ROOT}/example/covariates.txt"

[[ -f "$BGEN" ]] || { echo "missing $BGEN" >&2; exit 1; }

NVAR=$(python3 - "$BGEN" <<'EOF'
import struct, sys
with open(sys.argv[1], 'rb') as f:
    f.read(4)
    lh, m = struct.unpack('<II', f.read(8))
print(m)
EOF
)

echo "=== $LABEL ==="
echo "binary   : $BIN"
echo "bgen     : $BGEN ($NVAR variants)"
echo "threads  : $THREADS   reps: $REPS"

times=()
for i in $(seq 1 "$REPS"); do
  out=$(mktemp -d)
  start=$(python3 -c 'import time; print(time.time())')
  "$BIN" --step 2 --ignore-pred --bt \
    --bgen "$BGEN" --sample "$SAMPLE" \
    --phenoFile "$PHENO" --covarFile "$COVAR" \
    --bsize 1000 --threads "$THREADS" \
    --out "$out/bench" > "$out/log" 2>&1
  end=$(python3 -c 'import time; print(time.time())')
  el=$(python3 -c "print(f'{$end - $start:.3f}')")
  tp=$(python3 -c "print(f'{$NVAR / ($end - $start):.0f}')")
  times+=("$el")
  printf '  run %d: %8s s  %10s vars/s\n' "$i" "$el" "$tp"
  rm -rf "$out"
done

python3 - "$NVAR" "${times[@]}" <<'EOF'
import statistics, sys
nvar = int(sys.argv[1])
ts = sorted(float(x) for x in sys.argv[2:])
best, med = ts[0], statistics.median(ts)
print(f"  best   : {best:8.3f} s  {nvar/best:10.0f} vars/s")
print(f"  median : {med:8.3f} s  {nvar/med:10.0f} vars/s")
print(f"  spread : {(ts[-1]-ts[0])/ts[0]*100:.2f}%")
EOF
