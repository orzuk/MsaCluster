#!/usr/bin/env python3
"""Figure and summary for the concordance false-positive calibration.

Left: the per-family concordance of the real seven predictors against the
fold-free pseudo-methods, which carry sequence composition but no structural
information. Right: how many families each flags at BH q < 0.05.

Usage: plot_fp_calibration.py FPCAL_PER_FAMILY_CSV FPCAL_REPS_CSV OUT_STEM
"""
import os
import sys

os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib.pyplot as plt       # noqa: E402
import numpy as np                    # noqa: E402
import pandas as pd                   # noqa: E402

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REAL = os.path.join(REPO, "docs",
                    "fold_diversity_concordance_corrected_perm_7method_bh.csv")

fam_csv, reps_csv, stem = sys.argv[1], sys.argv[2], sys.argv[3]
real = pd.read_csv(REAL)
real = real[real.n_clusters >= 3]
fam = pd.read_csv(fam_csv)
reps = pd.read_csv(reps_csv)

fig, ax = plt.subplots(1, 2, figsize=(9.0, 3.4))

bins = np.linspace(0.40, 0.72, 33)
ax[0].hist(fam.observed, bins=bins, density=True, color="0.70",
           label=f"fold-free pseudo-methods ({len(fam)} draws)")
ax[0].hist(real.observed_mean_concordance_corrected, bins=bins, density=True,
           histtype="step", color="C3", lw=1.8,
           label=f"seven predictors ({len(real)} families)")
ax[0].axvline(0.5, color="k", ls=":", lw=1.0)
ax[0].set_xlabel("mean pairwise sign concordance (corrected)")
ax[0].set_ylabel("density")
ax[0].legend(fontsize=7, frameon=False)

ax[1].hist(reps.hits, bins=np.arange(-0.5, max(reps.hits.max() + 2, 14), 1),
           color="0.70", label="fold-free replicates")
ax[1].axvline(33, color="C3", lw=1.8, label="observed (33)")
ax[1].set_xlabel("families at BH $q<0.05$ (of 87)")
ax[1].set_ylabel("replicates")
ax[1].legend(fontsize=7, frameon=False)

for a in ax:
    a.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(f"{stem}.{ext}", dpi=200)

ro = real.observed_mean_concordance_corrected
print(f"real     : n={len(real)} median={ro.median():.4f} "
      f"mean={ro.mean():.4f} frac>0.5={(ro > 0.5).mean():.3f}")
print(f"fold-free: n={len(fam)} median={fam.observed.median():.4f} "
      f"mean={fam.observed.mean():.4f} frac>0.5={(fam.observed > 0.5).mean():.3f}")
print(f"hits     : observed 33   fold-free mean {reps.hits.mean():.1f} "
      f"median {reps.hits.median():.0f} range {reps.hits.min()}-{reps.hits.max()}")
print(f"wrote {stem}.png/.pdf")
