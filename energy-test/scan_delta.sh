# ---- inputs ---------------------------------------------------------------
DSET1="${1:-Sample3_20000.txt}"
DSET2="${2:-Sample2_20000.txt}"

# core controls
N_PERM="${N_PERM:-100000}"         # --n-permutations
N1_SUB="${N1_SUB:-150}"           # --max-permutation-events-1
SEED="${SEED:-12345}"             # --seed
# Generate DELTAS to run as a sequence using seq (start, step, end)
DELTAS=()
for ((i=15000; i<=15100; i+=10)); do
  DELTAS+=( "$(printf "%.2f" "$(bc -l <<< "$i/100")")" )
done
# --------------------------------------------------------------------------

# prompt for a filename base (for CSV and plots)
DEFAULT_BASE="$(basename "${DSET1%.txt}")_vs_$(basename "${DSET2%.txt}")_$(date +%Y%m%d_%H%M%S)"
read -r -p "[default: ${DEFAULT_BASE}] Enter a filename base (for output files): " BASE
BASE="${BASE:-$DEFAULT_BASE}"
# make it filesystem-friendly (spaces -> underscores)
BASE="${BASE// /_}"

# makes output locations if they don't exist
RESULTS_DIR="${RESULTS_DIR:-results}"
PLOTS_DIR="${PLOTS_DIR:-plots}"
mkdir -p "$RESULTS_DIR" "$PLOTS_DIR"

OUTCSV="${OUTCSV:-${RESULTS_DIR}/${BASE}_pvalues_scan.csv}"
OUTPNG="${OUTPNG:-${PLOTS_DIR}/${BASE}_pvalues_vs_delta.png}"
OUTPDF="${OUTPDF:-${PLOTS_DIR}/${BASE}_pvalues_vs_delta.pdf}"

# sanity checks
[[ -x ./energy_test ]] || { echo "energy_test not found or not executable."; exit 2; }
[[ -f "$DSET1" && -f "$DSET2" ]] || { echo "Input TXT files not found: $DSET1 / $DSET2"; exit 2; }

# CSV header
echo "dataset1,dataset2,delta,p_value,n_perm,n1_sub,seed" > "$OUTCSV"

# Archive old pvalues.txt if it exists, then clear it for new results
if [[ -f pvalues.txt ]]; then
  # Convert old pvalues.txt to CSV and append to archive
  awk '{printf "%s,%s\n", $1, $2}' pvalues.txt >> pvalues_archive.txt
  > pvalues.txt
fi

for DELTA in "${DELTAS[@]}"; do
  echo "[*] Running delta=${DELTA}  ($DSET1 vs $DSET2)"
  ./energy_test \
    --n-permutations "$N_PERM" \
    --max-permutation-events-1 "$N1_SUB" \
    --delta-value "$DELTA" \
    --seed "$SEED" \
    "$DSET1" "$DSET2"

  # last line in pvalues.txt: "<delta> <pvalue>"
  read -r PV_DELTA PVAL < <(tail -n 1 pvalues.txt)
  printf "%s,%s,%.6g,%.6g,%d,%d,%d\n" \
    "$(basename "$DSET1")" "$(basename "$DSET2")" \
    "$PV_DELTA" "$PVAL" "$N_PERM" "$N1_SUB" "$SEED" >> "$OUTCSV"
done

python3 plot_pvalues.py "$OUTCSV" "$OUTPNG" "$OUTPDF"

echo "[*] Done."
echo "    Results:  $OUTCSV"
echo "    Figures:  $OUTPNG  and  $OUTPDF"