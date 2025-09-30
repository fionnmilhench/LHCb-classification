# This file is called with python3 plot_pvalues.py "$OUTCSV" "$OUTPNG" "$OUTPDF"
import sys, csv
import matplotlib.pyplot as plt

# csv header: dataset1,dataset2,delta,p_value,n_perm,n1_sub,seed
csv_path, out_png, out_pdf = sys.argv[1], sys.argv[2], sys.argv[3]
xs, ys, fieldnames = [], [], []
with open(csv_path, newline="") as f:
    r = csv.DictReader(f)
    for row in r:
        if fieldnames == []:
            fieldnames = [row["dataset1"], row["dataset2"], row["n_perm"], row["n1_sub"], row["seed"]]
            # Fieldnames formatted as "Sample1_20000.txt" --> we want only "Sample1"
            fieldnames = [fn.split("_")[0] for fn in fieldnames]
        xs.append(float(row["delta"]))
        ys.append(float(row["p_value"]))

# sort by delta for a clean line
pairs = sorted(zip(xs, ys), key=lambda t: t[0])
xs = [p[0] for p in pairs]
ys = [p[1] for p in pairs]

plt.loglog(xs, ys, marker="o", label="{} vs {}".format(r.fieldnames and "", ""))  # simple label
plt.axhline(0.003, linestyle="--", linewidth=1, color="red", label="SL (0.3%)")
plt.xlabel(r"$\delta$")
plt.ylabel("p-value")
title = r"Energy test: $p$-value vs $\delta$ for {} vs {}".format(fieldnames[0], fieldnames[1])
plt.title(title)
legend_title = r"$N_{{\text{{perm}}}}:$ {:,}  $N_{{1,\text{{sub}}}}:$ {:,}".format(int(fieldnames[2]), int(fieldnames[3])) 
plt.legend(title=legend_title, loc ="lower right")
plt.tight_layout()
plt.savefig(out_png, dpi=300)
plt.savefig(out_pdf)