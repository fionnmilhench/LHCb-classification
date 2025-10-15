import argparse
import csv, os, sys
import numpy as np
from scipy.stats import chi2 as chi2_dist

def compute_histNd(X, nbins, ranges):
    # uniform binning in all dimensions
    H, _ = np.histogramdd(X, bins=[nbins] * X.shape[1], range=ranges)
    return H.astype(np.int64)

p = argparse.ArgumentParser()
p.add_argument("dataset1", type=str)
p.add_argument("dataset2", type=str)
p.add_argument("--nbins", type=int, default=40)
p.add_argument("--out-csv", type=str, default="miranda_results.csv")
args = p.parse_args()

# Load data
d1 = np.loadtxt(args.dataset1, ndmin=2)
d2 = np.loadtxt(args.dataset2, ndmin=2)

if d1.size == 0 or d2.size == 0 or d1.shape[0] == 0 or d2.shape[0] == 0:
    print("Error: one of the datasets is empty after parsing.", file=sys.stderr)
    sys.exit(2)
if d1.shape[1] != d2.shape[1]:
    print(f"Error: datasets have different dimensionality ({d1.shape[1]} vs {d2.shape[1]}).", file=sys.stderr)
    sys.exit(2)

D = d1.shape[1] # dimensionality
# We need min/max per dimension from all data for histogram binning
combined = np.vstack([d1, d2])
ranges = []
for j in range(D):
    vmin = combined[:, j].min()
    vmax = combined[:, j].max()
    if not (np.isfinite(vmin) and np.isfinite(vmax)):
        print("Error: non-finite binning ranges.", file=sys.stderr)
        sys.exit(2)
    if vmax <= vmin:
        print(f"Error: invalid range in dim {j} (max <= min).", file=sys.stderr)
        sys.exit(2)
    ranges.append((float(vmin), float(vmax)))

# N-D histograms, then flatten
H1 = compute_histNd(d1, args.nbins, ranges)
H2 = compute_histNd(d2, args.nbins, ranges)

n1 = H1.ravel()
n2 = H2.ravel()
total_per_bin = n1 + n2

# Keep bins with at least one count
keep = total_per_bin > 0

bins = int(keep.sum())
low_count_bins = (total_per_bin > 0) & (total_per_bin < 20)
low_count_bins = int(low_count_bins.sum())
if bins == 0:
    print("Error: all bins are empty.", file=sys.stderr)
    sys.exit(2)

# "Miranda" procedure test statistic implementation (details: https://doi.org/10.1103/PhysRevD.80.096006)
# test statistic S_i = (n1_i - n2_i) / sqrt(n1_i + n2_i) where n1_i, n2_i are counts in bin i
# chi2 = sum_i S_i^2

num = n1[keep] - n2[keep]
den = np.sqrt(n1[keep] + n2[keep])

with np.errstate(divide="ignore", invalid="ignore"):
    S2 = np.zeros_like(den, dtype=float)
    mask = den > 0
    S2[mask] = (num[mask] / den[mask]) ** 2

chi2 = float(S2.sum())

# No fitted parameters, so dof = bins
dof = bins
if dof <= 0:
    print("Error: dof <= 0, cannot compute p-value.", file=sys.stderr)
    sys.exit(2)

p_value = float(chi2_dist.sf(chi2, dof))


# Write/append CSV, report status
header = ["dataset1","dataset2","nbins","chi2","p_value","bins","low_count_bins"]
row = [os.path.basename(args.dataset1), os.path.basename(args.dataset2), 
    args.nbins, chi2, p_value, bins, low_count_bins,
]

file_exists = os.path.isfile(args.out_csv)
with open(args.out_csv, "a", newline="") as f:
    w = csv.writer(f)
    if not file_exists:
        w.writerow(header)
    w.writerow(row)

print(f"Miranda test completed: chi2 = {chi2:.2f}, p-value = {p_value:.3e}, bins = {bins}, low count bins = {low_count_bins}")