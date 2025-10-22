#!/usr/bin/env bash
set -euo pipefail

# Arguments parsing
usage() {
  echo "Usage: $0 --d1 DSET1 --d2 DSET2 --nperm N --n1sub N --seed N --delta X.YY|--deltas X,Y,Z [--kernel gaussian|lognormal|step|lorentzian] [--base NAME]"
  exit 2
}
DSET1=""; DSET2=""; N_PERM=""; N1_SUB=""; SEED=""; DELTA=""; DELTAS=""; BASE="energy_scan"; KERNEL="gaussian"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --d1)    DSET1="$2"; shift 2;;
    --d2)    DSET2="$2"; shift 2;;
    --nperm) N_PERM="$2"; shift 2;;
    --n1sub) N1_SUB="$2"; shift 2;;
    --seed)  SEED="$2"; shift 2;;
    --delta) DELTA="$2"; shift 2;;
    --deltas) DELTAS="$2"; shift 2;;
    --kernel) KERNEL="$2"; shift 2;;
    --base)  BASE="$2"; shift 2;;
    *) echo "Unknown arg: $1"; usage;;
  esac
done
# Require either DELTA or DELTAS, and not kernel
[[ -n "$DSET1" && -n "$DSET2" && -n "$N_PERM" && -n "$N1_SUB" && -n "$SEED" && ( -n "$DELTA" || -n "$DELTAS" ) ]] || usage
# ----------------------------------------------------------------------

# Diagnostics
[[ -x ./energy_test ]] || { echo "energy_test not found or not executable"; exit 3; }
[[ -f "$DSET1" && -f "$DSET2" ]] || { echo "Input files missing: $DSET1 / $DSET2"; exit 4; }

process_delta() {
  local DEL=$1
  ./energy_test \
    --n-permutations "$N_PERM" \
    --max-permutation-events-1 "$N1_SUB" \
    --delta-value "$DEL" \
    --kernel "$KERNEL" \
    --seed "$SEED" \
    "$DSET1" "$DSET2"

  # Extract "<delta> <pvalue>" from pvalues.txt
  read -r PV_DELTA PVAL < <(tail -n 1 pvalues.txt)

  OUTFILE="${BASE}.csv"
  printf "%s,%s,%.12g,%.12g,%d,%d,%d,%s\n" \
    "$(basename "$DSET1")" \
    "$(basename "$DSET2")" \
    "$PV_DELTA" "$PVAL" "$N_PERM" "$N1_SUB" "$SEED" "$KERNEL" >> "$OUTFILE"
}

if [[ -n "$DELTAS" ]]; then
  IFS=',' read -r -a arr <<< "$DELTAS"
  for d in "${arr[@]}"; do
    process_delta "$d"
  done
else
  process_delta "$DELTA"
fi