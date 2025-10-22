import argparse
import csv, os, sys
import numpy as np
from scipy.stats import chi2 as chi2_dist
import time as t
t1 = t.time()

def compute_histNd(X, nbins, ranges):
    # uniform binning in all dimensions
    H, _ = np.histogramdd(X, bins=[nbins] * X.shape[1], range=ranges)
    return H.astype(np.int64)

p = argparse.ArgumentParser()
p.add_argument("dataset1", type=str)
p.add_argument("dataset2", type=str)
p.add_argument("--nbins", type=int, default=40)
p.add_argument("--out", type=str, default="miranda_results.csv")
p.add_argument("--nperm", type=int, default=0, help="(0: skip emp. p-value + use SciPy p-value only)")
p.add_argument("--seed", type=int, default=2025)
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

p_value_scipy = float(chi2_dist.sf(chi2, dof))

# Pooled-null permutation "empirical" p-value as discussed with Will
# Combine counts from both datasets, then randomly assign counts to two new datasets
# of same sizes as original datasets, compute Miranda chi2 for each permutation
# Count how many times permuted chi2 >= observed chi2, 
# divide by total permutations to get empirical p-value
if args.nperm and args.nperm > 0:
    pooled_counts = n1 + n2
    nonempty_mask = pooled_counts > 0
    pooled_nonempty = pooled_counts[nonempty_mask].astype(np.int64)

    # Observed Miranda statistic calculated on non-empty bins
    denom = np.sqrt(pooled_nonempty) # shape (nbins_nonempty,)
    T_obs = chi2 # as calculated previously

    n1_total = int(n1.sum())
    n_total = int(pooled_nonempty.sum())
    if not (0 <= n1_total <= n_total):
        print("Error: invalid totals for permutation.", file=sys.stderr)
        sys.exit(2)

    rng = np.random.default_rng(args.seed)
    n_perms = int(args.nperm)

    # draw: shape (n_perms, K)
    labels1_draws = rng.multivariate_hypergeometric(pooled_nonempty, n1_total, size=n_perms)

    # Compute Miranda T for every permutation
    diff_draws = labels1_draws - (pooled_nonempty - labels1_draws)   # diff between labels1, and remainder
    T_draws = ((diff_draws / denom) ** 2).sum(axis=1)

    # Empirical tail with add-one smoothing
    ge = np.count_nonzero(T_draws >= T_obs) 
    p_value_empirical = float((ge + 1.0) / (n_perms + 1.0)) # add-one because the observed is a permutation too + prevents p=0

t2 = t.time()
time_taken = float(t2 - t1)

# Write/append CSV, report status
header = ["dataset1","dataset2","nbins","chi2","p_value_scipy","p_value_empirical","bins","low_count_bins"]
row = [os.path.basename(args.dataset1), os.path.basename(args.dataset2), 
    args.nbins, chi2, p_value_scipy, p_value_empirical, bins, low_count_bins,
]

file_exists = os.path.isfile(args.out)
with open(args.out, "a", newline="") as f:
    w = csv.writer(f)
    if not file_exists:
        w.writerow(header)
    w.writerow(row)

if args.nperm and args.nperm > 0:
    print(f"Miranda test completed in {time_taken:.2f} seconds: chi2 = {chi2:.2f}, p-value (SciPy) = {p_value_scipy:.3e}, p-value (empirical) = {p_value_empirical:.3e}, bins = {bins}, low count bins = {low_count_bins}")
else:
    print(f"Miranda test completed in {time_taken:.2f} seconds: chi2 = {chi2:.2f}, p-value (SciPy) = {p_value_scipy:.3e}, bins = {bins}, low count bins = {low_count_bins}")