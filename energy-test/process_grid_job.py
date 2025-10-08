#!/usr/bin/env python3
from pathlib import Path
import argparse
import sys

HEADER = "dataset1,dataset2,delta,p_value,n_perm,n1_sub,seed"
EXPECTED_COLS = 7

def is_header_line(line: str) -> bool:
    s = line.strip()
    if not s:
        return True
    if s.lower().startswith("dataset1,"):
        return True

ap = argparse.ArgumentParser()
ap.add_argument("job_number", help="Job number under LocalXML (e.g. 9)")
ap.add_argument("--base", default="/afs/cern.ch/user/f/fmilhenc/gangadir/workspace/fmilhenc/LocalXML")
ap.add_argument("--outname", default=None)
args = ap.parse_args()

job_dir = Path(args.base) / str(args.job_number)
if not job_dir.is_dir():
    print(f"Error: job directory not found: {job_dir}", file=sys.stderr)
    sys.exit(1)

csvs = sorted(job_dir.glob("*/output/*.csv"), key=lambda p: p.name)
if not csvs:
    print(f"No CSV files under {job_dir}/*/output/", file=sys.stderr)
    sys.exit(2)

# Results dir relative to this script's location
script_dir = Path(__file__).resolve().parent
results_dir = script_dir / "results"
results_dir.mkdir(parents=True, exist_ok=True)

out_file = args.outname or f"job{args.job_number}_processed.csv"
out_path = results_dir / out_file

rows = 0
with out_path.open("w", newline="") as fout:
    fout.write(HEADER + "\n")
    for f in csvs:
        with f.open("r") as fin:
            for line in fin:
                line = line.strip()
                if not line or is_header_line(line):
                    continue
                if line.count(",") != (EXPECTED_COLS - 1):
                    continue
                fout.write(line + "\n")
                rows += 1

print(f"Wrote {rows} rows to {out_path}")