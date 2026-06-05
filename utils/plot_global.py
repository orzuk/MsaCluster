"""Cross-pair global statistics plots (per-fold support panel + Moran's plots)."""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from config import TABLES_RES
from utils.method_labels import METHOD_ORDER, disp, order_methods, metric_label

SURVEY_CSV = os.path.join(TABLES_RES, "fold_diversity_survey.csv")

# Git-TRACKED output dir for the global figures (TABLES_RES is docs/). Writing
# here -- not the gitignored Pipeline/Results/Figures -- means a plain
# `git add docs/ && commit && push` updates the GitHub-Pages global page.
# Paths in GLOBAL_FIGURE_SPEC are relative to docs/ and point at figures/global/.
GLOBAL_FIG_DIR = os.path.join(TABLES_RES, "figures", "global")

# Experimentally validated fold-switchers, for highlighting on global plots.
KNOWN_POSITIVES = {
    "2n54B_2hdmA": "XCL1",
    "1zk9A_3jv6A": "RelB",
    "2namA_1uxmK": "SOD1",
    "5jytA_2qkeE": "KaiB",
}

# ---------------------------------------------------------------------
# Ordered + captioned spec for the global analysis HTML page. ONE place
# defines which global figures appear, in what order, and the explanatory
# text shown before each (Rmarkdown-style). Paths are relative to docs/.
# ---------------------------------------------------------------------
GLOBAL_FIGURE_SPEC = [
    {
        "file": "figures/global/global_per_fold_support_panel.png",
        "title": "Per-pair cluster distribution in the fold-1/fold-2 plane",
        "caption":
            "Each panel is one method; each ELLIPSE is one fold-switch pair — a "
            "1σ 2D-Gaussian fit to that pair's sequence clusters, plotting support "
            "for fold 1 (x) vs fold 2 (y) in the method's native metric (TM-score "
            "for AF2/AF3/ESMFold2/Boltz-2; −ΣΔΔG kcal/mol for DDG; SS H/E-composition "
            "similarity for S4PRED; contact recall for MSA-Transformer/CCMpred). "
            "The ellipse's ORIENTATION is the key: elongated PERPENDICULAR to the "
            "y=x diagonal means the pair's clusters split between the two folds — "
            "the fold-switching signal — whereas elongation ALONG the diagonal is "
            "shared variation with no fold discrimination, and a small round "
            "ellipse means low cluster diversity. Ellipses are colored by their "
            "perpendicular spread (brighter = stronger fold-switch signal), and the "
            "few highest-signal pairs per method are labelled. Pairs with fewer "
            "than 5 clusters are omitted (no stable covariance). Note: TM axes are "
            "comparable across pairs (0–1); for DDG/recall the absolute position is "
            "pair-specific, so read orientation, not location.",
    },
    {
        "file": "figures/global/global_concordance_volcano.png",
        "title": "Concordance axis — cross-method agreement per pair",
        "caption":
            "The second analysis axis (the one that yields hits): for each pair, "
            "how strongly do the orthogonal methods AGREE on a per-cluster fold "
            "preference. x = mean cross-method concordance (sequence-divergence "
            "corrected); y = −log10 permutation p. Pairs passing full-pool BH are "
            "ringed; the four experimentally validated fold-switchers are starred. "
            "KaiB is strongly concordant; XCL1 and RelB are the within-class "
            "(stratified-BH) recoveries — see the trigger-class figure below.",
    },
    {
        "file": "figures/global/global_raw_vs_corrected.png",
        "title": "Three roles of the sequence-divergence correction",
        "caption":
            "Per pair, significance before (x) vs after (y) correcting for "
            "cluster-to-query sequence divergence. Relative to the p=0.05 lines: "
            "top-right = significant both ways (PRESERVE — a true signal survives, "
            "e.g. XCL1), bottom-right = significant only before (FILTER — a "
            "divergence-driven false positive removed, e.g. SOD1), top-left = "
            "significant only after (REVEAL — a real signal unmasked by the "
            "correction, e.g. RelB). These three canonical cases anchor the "
            "correction's behaviour.",
    },
    {
        "file": "figures/global/global_trigger_stratification.png",
        "title": "Concordance signal by fold-switch trigger class",
        "caption":
            "Corrected concordance distribution within each biological trigger "
            "class (ligand / oligomerization / protein-binding / equilibrium / "
            "mutation), the basis for the stratified (within-class) BH test. "
            "Validated positives are starred in their classes — XCL1 in mutation, "
            "RelB in protein-binding, SOD1 in oligomerization — showing how "
            "stratifying by trigger recovers hits that the full-pool test misses.",
    },
    {
        "file": "figures/global/global_morans_method_x_pair.png",
        "title": "Phylogenetic clade signal — method × pair",
        "caption":
            "Sequence-divergence-corrected clade signal (Moran's I_norm) for every "
            "method (columns) and fold-switch pair (rows). Red = clusters that "
            "prefer the same fold are phylogenetically grouped (clade-structured); "
            "blue = anti-structured; white ≈ none. Rows are sorted by the median "
            "signal across methods (robust to the sparse, noisy contact-based "
            "methods), so the most clade-structured pairs sit at the top. Signal "
            "is weak overall — few cells exceed I_norm ≈ 0.2.",
    },
    {
        "file": "figures/global/global_method_correlation.png",
        "title": "Method × method agreement on clade signal",
        "caption":
            "Spearman correlation between methods of their per-pair clade signal. "
            "High values mean two methods light up on the same pairs. Only the "
            "structure predictors (AF2/AF3/Boltz-2/ESMFold2) correlate "
            "appreciably; ΔΔG, MSA-Transformer, CCMpred and S4PRED are largely "
            "orthogonal, i.e. they contribute independent evidence.",
    },
    {
        "file": "figures/global/global_morans_effectsize_by_method.png",
        "title": "Clade-signal effect size, per method",
        "caption":
            "Distribution across pairs of the clade-signal effect size (Moran's "
            "I_norm) for each method. Boxes near zero indicate little systematic "
            "clade structure; ΔΔG shows the widest positive tail. Because n ≈ 100 "
            "clusters makes the permutation test overpowered, we rank methods by "
            "effect size here rather than by significance counts.",
    },
    {
        "file": "figs/cross_method_correlation.png",
        "title": "Cross-method correlation of raw fold preference",
        "caption":
            "Spearman correlation of the raw per-cluster fold preference "
            "(centered TM difference), pooled across all clusters and pairs — a "
            "signal-level companion to the clade-signal correlation above. AF2 and "
            "AF3 share training bias (highest correlation); every other method pair "
            "is near zero, confirming the methods supply largely independent "
            "evidence (most orthogonal: S4PRED).",
    },
]


# =====================================================================
# Per-fold support panel — one scatter per method (fold1 vs fold2)
# =====================================================================

def plot_per_fold_support_panel(output_dir: str, agg: str = "max",
                                n_sigma: float = 1.0, min_clusters: int = 5,
                                label_top: int = 3) -> None:
    """2x4 panel: per method, ONE 1-sigma ellipse per pair in the (fold 1, fold 2)
    support plane.

    Each pair's clusters are summarized by a 2D Gaussian (mean + covariance); the
    ellipse shows where the pair sits, how diverse its clusters are, and -- via its
    orientation -- WHICH WAY. Elongation PERPENDICULAR to y=x means the pair's
    clusters split between the two folds (fold-switching signal); elongation ALONG
    y=x means shared variation with no fold discrimination; a small round ellipse
    means low diversity. Ellipses are colored by their perpendicular spread (the
    fold-switch signal), and the `label_top` pairs with the largest perpendicular
    spread are annotated per method. Support is each method's native per-fold
    metric (TM-score / -sum ddG / SS-similarity / contact recall).
    """
    from matplotlib.patches import Ellipse
    from matplotlib import cm, colors as mcolors

    if not (os.path.isfile(SURVEY_CSV) and os.path.getsize(SURVEY_CSV) > 0):
        print(f"[global] no {SURVEY_CSV}; skipping per-fold support panel.")
        return
    df = pd.read_csv(SURVEY_CSV)
    x_col, y_col = f"TM1_{agg}", f"TM2_{agg}"
    if x_col not in df.columns or y_col not in df.columns:
        print(f"[global] survey lacks {x_col}/{y_col}; skipping per-fold panel.")
        return

    methods = order_methods(df["method"].unique())
    ncol = 4
    nrow = int(np.ceil(len(methods) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 4.2 * nrow),
                             constrained_layout=True)
    axes = np.atleast_1d(axes).ravel()
    cmap = cm.get_cmap("plasma")

    # unit vector perpendicular to y=x; projection length = fold-1 - fold-2 signal
    perp = np.array([1.0, -1.0]) / np.sqrt(2.0)

    for ax, m in zip(axes, methods):
        s = df[df["method"] == m]
        # per pair: mean, covariance, perpendicular spread
        ells = []          # (mx, my, cov, perp_std, pair_id)
        allx, ally = [], []
        for pid, g in s.groupby("pair_id"):
            x = pd.to_numeric(g[x_col], errors="coerce")
            y = pd.to_numeric(g[y_col], errors="coerce")
            ok = x.notna() & y.notna()
            x, y = x[ok].values, y[ok].values
            if len(x) < min_clusters or (x.std() == 0 and y.std() == 0):
                continue
            P = np.column_stack([x, y])
            cov = np.cov(P, rowvar=False)
            perp_std = float(np.sqrt(max((perp @ cov @ perp), 0.0)))
            ells.append((x.mean(), y.mean(), cov, perp_std, pid))
            allx += list(x); ally += list(y)

        if not ells:
            ax.set_title(f"{disp(m)}  (no pairs ≥{min_clusters} clusters)", fontsize=9)
            ax.axis("off"); continue

        pstd = np.array([e[3] for e in ells])
        norm = mcolors.Normalize(vmin=float(pstd.min()), vmax=float(pstd.max()) or 1.0)

        for mx, my, cov, ps, pid in ells:
            vals, vecs = np.linalg.eigh(cov)        # ascending
            vals = np.clip(vals, 0, None)
            order = vals.argsort()[::-1]
            vals, vecs = vals[order], vecs[:, order]
            w, h = 2.0 * n_sigma * np.sqrt(vals)     # full axis lengths
            ang = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
            col = cmap(norm(ps))
            ax.add_patch(Ellipse((mx, my), max(w, 1e-9), max(h, 1e-9), angle=ang,
                                  fill=False, edgecolor=col, lw=1.1, alpha=0.8))

        # robust limits from the clusters that contributed
        both = np.array(allx + ally, dtype=float)
        lo, hi = float(np.nanpercentile(both, 1)), float(np.nanpercentile(both, 99))
        pad = 0.06 * (hi - lo or 1.0)
        ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad],
                color="gray", lw=0.8, ls="--", alpha=0.7)
        ax.set_xlim(lo - pad, hi + pad); ax.set_ylim(lo - pad, hi + pad)

        # annotate the few highest-perpendicular-spread (fold-switch-like) pairs
        for mx, my, cov, ps, pid in sorted(ells, key=lambda e: e[3], reverse=True)[:label_top]:
            ax.annotate(pid, (mx, my), fontsize=5.5, ha="center", va="center",
                        color="black")

        mm = metric_label(m)
        ax.set_title(f"{disp(m)}  ({len(ells)} pairs)", fontsize=10)
        ax.set_xlabel(f"{mm} (fold 1)"); ax.set_ylabel(f"{mm} (fold 2)")
        ax.grid(True, alpha=0.25)

    for ax in axes[len(methods):]:
        ax.axis("off")

    sm = cm.ScalarMappable(cmap=cmap, norm=mcolors.Normalize(0, 1))
    cbar = fig.colorbar(sm, ax=axes.tolist(), shrink=0.5, pad=0.01)
    cbar.set_label("perpendicular spread → fold-switch signal\n(relative, per method)",
                   fontsize=9)
    fig.suptitle(f"Per-pair cluster distribution in the fold-1/fold-2 plane "
                 f"({int(n_sigma)}σ ellipse per pair, {agg} over models)", fontsize=13)
    out = os.path.join(output_dir, "global_per_fold_support_panel.png")
    fig.savefig(out, dpi=170); plt.close(fig)
    print(f"[plot] wrote {out}")


# =====================================================================
# Phylogenetic-signal (Moran's) global plots — read docs/phylo_placement*.csv
# =====================================================================

def _load_phylo(labels: str = "corrected"):
    """Load the Moran's phylo-placement table; return (df, metric) or (None, None).

    metric is the best available effect/signal column: I_norm > z > I.
    """
    fn = ("phylo_placement_corrected.csv" if labels == "corrected"
          else "phylo_placement.csv")
    path = os.path.join(TABLES_RES, fn)
    if not (os.path.isfile(path) and os.path.getsize(path) > 0):
        print(f"[global] no {fn}; skipping Moran's plots.")
        return None, None
    d = pd.read_csv(path)
    # I_norm is a fraction of achievable clade signal and is bounded to [-1, 1].
    # Older CSVs (pre ceiling-clamp fix in phylo_placement.py) can carry a few
    # |I_norm|>1 values from the heuristic max_perm underestimating the true
    # ceiling on tiny/sparse pairs (MSAT/CCMpred). Clip defensively so figures
    # never show impossible values; harmless once the source CSV is regenerated.
    if "morans_I_norm" in d.columns:
        d["morans_I_norm"] = pd.to_numeric(d["morans_I_norm"], errors="coerce").clip(-1.0, 1.0)
    for m in ("morans_I_norm", "morans_z", "morans_I"):
        if m in d.columns and d[m].notna().any():
            return d, m
    print(f"[global] {fn} has no morans_* columns (old thresholded file?); "
          f"re-run phylo_placement. Skipping Moran's plots.")
    return None, None


def plot_method_clade_heatmap(output_dir: str, labels: str = "corrected",
                              row_sort: str = "median") -> None:
    """Methods x pairs heatmap of clade-signal strength (effect size / z).

    Rows (pairs) are sorted by a method-agnostic summary of the signal across
    methods (descending), so the most clade-structured pairs rise to the top.
    `row_sort` = "median" (default, robust to a couple of noisy methods like
    MSA-Transformer/CCMpred), "mean", or "max".
    """
    d, metric = _load_phylo(labels)
    if d is None:
        return
    wide = d.pivot_table(index="pair_id", columns="method", values=metric, aggfunc="first")
    cols = order_methods(list(wide.columns))
    wide = wide[cols]
    # Per-pair winner tally (printed for the "is DDG strongest?" question)
    win = wide.idxmax(axis=1).value_counts()
    print(f"[global] strongest method per pair ({metric}): "
          + ", ".join(f"{m}:{int(win.get(m,0))}" for m in cols))
    # Method-agnostic ordering: median (robust) / mean / max across methods.
    agg_fn = {"median": wide.median, "mean": wide.mean, "max": wide.max}.get(
        row_sort, wide.median)
    w = wide.assign(_key=agg_fn(axis=1, skipna=True)).sort_values(
        "_key", ascending=False).drop(columns="_key")
    M = w.to_numpy(dtype=float)
    vlim = float(np.nanpercentile(np.abs(M), 98)) or 1.0
    fig, ax = plt.subplots(figsize=(6, max(4, 0.16 * len(w))))
    im = ax.imshow(M, cmap="RdBu_r", vmin=-vlim, vmax=vlim, aspect="auto")
    ax.set_xticks(range(len(cols))); ax.set_xticklabels([disp(c) for c in cols],
                                                        rotation=45, ha="right")
    ax.set_yticks(range(len(w))); ax.set_yticklabels(w.index, fontsize=4)
    ax.set_title(f"Clade signal per method x pair ({metric}, {labels})\n"
                 f"rows = pairs, sorted by {row_sort} signal across methods",
                 fontsize=10)
    fig.colorbar(im, ax=ax, shrink=0.6, label=metric)
    fig.tight_layout()
    out = os.path.join(output_dir, "global_morans_method_x_pair.png")
    fig.savefig(out, dpi=180); plt.close(fig)
    print(f"[plot] wrote {out}")


def plot_method_correlation_heatmap(output_dir: str, labels: str = "corrected") -> None:
    """Method x method Spearman correlation of the clade signal across pairs."""
    d, metric = _load_phylo(labels)
    if d is None:
        return
    wide = d.pivot_table(index="pair_id", columns="method", values=metric, aggfunc="first")
    cols = order_methods(list(wide.columns))
    wide = wide[cols]
    labs = [disp(c) for c in cols]
    C = wide.corr(method="spearman")
    fig, ax = plt.subplots(figsize=(6.2, 5.4))
    im = ax.imshow(C.to_numpy(), cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(len(cols))); ax.set_xticklabels(labs, rotation=45, ha="right")
    ax.set_yticks(range(len(cols))); ax.set_yticklabels(labs)
    for i in range(len(cols)):
        for j in range(len(cols)):
            v = C.iloc[i, j]
            if pd.notna(v):
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=6,
                        color="black" if abs(v) < 0.6 else "white")
    ax.set_title(f"Method x method correlation of clade signal\n(Spearman of {metric}, {labels})",
                 fontsize=10)
    fig.colorbar(im, ax=ax, shrink=0.7, label="Spearman ρ")
    fig.tight_layout()
    out = os.path.join(output_dir, "global_method_correlation.png")
    fig.savefig(out, dpi=180); plt.close(fig)
    print(f"[plot] wrote {out}")


def plot_morans_effectsize_by_method(output_dir: str, labels: str = "corrected") -> None:
    """Distribution of per-pair clade effect size per method (boxplot)."""
    d, metric = _load_phylo(labels)
    if d is None:
        return
    cols = order_methods(list(d["method"].unique()))
    data = [pd.to_numeric(d.loc[d["method"] == m, metric], errors="coerce").dropna().values
            for m in cols]
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.boxplot(data, labels=[disp(c) for c in cols], showfliers=True)
    ax.axhline(0, color="gray", lw=0.8, ls="--")
    ax.set_ylabel(metric)
    ax.set_title(f"Clade signal effect size per method across pairs ({metric}, {labels})",
                 fontsize=10)
    ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()
    out = os.path.join(output_dir, "global_morans_effectsize_by_method.png")
    fig.savefig(out, dpi=180); plt.close(fig)
    print(f"[plot] wrote {out}")


def _neglog10p(p, n_perm=1000):
    """-log10(p) with a permutation-resolution floor (p can be 0 at 1/(n_perm+1))."""
    floor = 1.0 / (n_perm + 1)
    return -np.log10(np.clip(pd.to_numeric(p, errors="coerce"), floor, 1.0))


def plot_concordance_volcano(output_dir: str) -> None:
    """Volcano of the cross-method concordance axis (the test that yields hits).

    x = per-pair mean cross-method concordance (sequence-divergence corrected);
    y = -log10 permutation p. Full-pool BH-significant pairs are ringed; the four
    experimentally validated fold-switchers are labelled.
    """
    f = os.path.join(TABLES_RES, "fold_diversity_concordance_corrected_perm_bh.csv")
    if not os.path.isfile(f):
        print(f"[global] no {f}; skipping concordance volcano."); return
    d = pd.read_csv(f)
    x = pd.to_numeric(d["observed_mean_concordance_corrected"], errors="coerce")
    y = _neglog10p(d["perm_mean_concordance_p"])
    bh = pd.to_numeric(d.get("p_bh", pd.Series([np.nan] * len(d))), errors="coerce")
    sig = bh < 0.05

    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    ax.scatter(x[~sig], y[~sig], s=18, c="#7a7a7a", alpha=0.6, edgecolors="none",
               label="n.s. (full-pool BH)")
    ax.scatter(x[sig], y[sig], s=46, facecolors="none", edgecolors="#d62728",
               linewidths=1.6, label="BH-significant (full pool)")
    ax.axhline(_neglog10p(pd.Series([0.05]))[0], color="gray", ls="--", lw=0.8)
    ax.text(ax.get_xlim()[0], _neglog10p(pd.Series([0.05]))[0], " p=0.05",
            fontsize=7, va="bottom", color="gray")
    for pid, name in KNOWN_POSITIVES.items():
        r = d[d["pair_id"] == pid]
        if not len(r):
            continue
        xi = float(pd.to_numeric(r["observed_mean_concordance_corrected"].iloc[0]))
        yi = float(_neglog10p(r["perm_mean_concordance_p"]).iloc[0])
        ax.scatter([xi], [yi], s=70, marker="*", color="#1f77b4", zorder=5)
        ax.annotate(name, (xi, yi), fontsize=8, fontweight="bold",
                    xytext=(4, 4), textcoords="offset points")
    ax.set_xlabel("mean cross-method concordance (corrected)")
    ax.set_ylabel("-log10  permutation p")
    ax.set_title("Concordance axis: cross-method agreement per pair", fontsize=11)
    ax.legend(fontsize=8, loc="lower left"); ax.grid(True, alpha=0.25)
    fig.tight_layout()
    out = os.path.join(output_dir, "global_concordance_volcano.png")
    fig.savefig(out, dpi=180); plt.close(fig)
    print(f"[plot] wrote {out}")


def plot_raw_vs_corrected(output_dir: str) -> None:
    """The three roles of the sequence-divergence correction.

    Per pair: -log10 p before (x) vs after (y) the correction. Quadrants vs the
    p=0.05 lines: top-right = significant both (PRESERVE), bottom-right =
    significant only raw (FILTER, false positive removed), top-left = significant
    only corrected (REVEAL, false negative unmasked). XCL1/SOD1/RelB are the
    canonical preserve/filter/reveal cases.
    """
    fr = os.path.join(TABLES_RES, "fold_diversity_concordance_perm_bh.csv")
    fc = os.path.join(TABLES_RES, "fold_diversity_concordance_corrected_perm_bh.csv")
    if not (os.path.isfile(fr) and os.path.isfile(fc)):
        print("[global] missing concordance perm_bh CSV(s); skipping raw-vs-corrected.")
        return
    raw = pd.read_csv(fr)[["pair_id", "perm_mean_concordance_p"]].rename(
        columns={"perm_mean_concordance_p": "p_raw"})
    cor = pd.read_csv(fc)[["pair_id", "perm_mean_concordance_p"]].rename(
        columns={"perm_mean_concordance_p": "p_cor"})
    m = raw.merge(cor, on="pair_id")
    x = _neglog10p(m["p_raw"]); y = _neglog10p(m["p_cor"])
    thr = _neglog10p(pd.Series([0.05]))[0]

    fig, ax = plt.subplots(figsize=(6.8, 6.4))
    ax.scatter(x, y, s=18, c="#7a7a7a", alpha=0.55, edgecolors="none")
    lim = max(float(np.nanmax(x)), float(np.nanmax(y))) * 1.08
    ax.plot([0, lim], [0, lim], color="gray", lw=0.8, ls=":")
    ax.axvline(thr, color="gray", lw=0.8, ls="--")
    ax.axhline(thr, color="gray", lw=0.8, ls="--")
    ax.text(thr, lim, " p=0.05 raw", rotation=90, va="top", fontsize=7, color="gray")
    ax.text(lim, thr, "p=0.05 corrected ", ha="right", va="bottom", fontsize=7, color="gray")
    roles = {"2n54B_2hdmA": "XCL1 (preserve)", "2namA_1uxmK": "SOD1 (filter)",
             "1zk9A_3jv6A": "RelB (reveal)"}
    for pid, lab in roles.items():
        r = m[m["pair_id"] == pid]
        if not len(r):
            continue
        xi = float(_neglog10p(r["p_raw"]).iloc[0]); yi = float(_neglog10p(r["p_cor"]).iloc[0])
        ax.scatter([xi], [yi], s=70, marker="*", color="#d62728", zorder=5)
        ax.annotate(lab, (xi, yi), fontsize=8, fontweight="bold",
                    xytext=(4, 4), textcoords="offset points")
    ax.set_xlim(0, lim); ax.set_ylim(0, lim)
    ax.set_xlabel("-log10 p  (raw)")
    ax.set_ylabel("-log10 p  (sequence-divergence corrected)")
    ax.set_title("Three roles of the correction: preserve / filter / reveal", fontsize=11)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    out = os.path.join(output_dir, "global_raw_vs_corrected.png")
    fig.savefig(out, dpi=180); plt.close(fig)
    print(f"[plot] wrote {out}")


def plot_trigger_stratification(output_dir: str) -> None:
    """Per-pair concordance signal grouped by fold-switch trigger class.

    Shows the corrected concordance distribution within each trigger class
    (the stratified-BH axis), with the validated positives labelled in their
    classes (XCL1 → mutation, RelB → protein_binding, SOD1 → oligomerization).
    """
    fc = os.path.join(TABLES_RES, "fold_diversity_concordance_corrected_perm_bh.csv")
    ft = os.path.join(TABLES_RES, "triggers_from_pdb.csv")
    if not (os.path.isfile(fc) and os.path.isfile(ft)):
        print("[global] missing concordance or triggers CSV; skipping trigger plot.")
        return
    d = pd.read_csv(fc)[["pair_id", "observed_mean_concordance_corrected"]]
    t = pd.read_csv(ft)[["pair_id", "trigger_class"]]
    m = d.merge(t, on="pair_id").dropna(subset=["trigger_class"])
    m["val"] = pd.to_numeric(m["observed_mean_concordance_corrected"], errors="coerce")
    m = m.dropna(subset=["val"])
    order = (m.groupby("trigger_class")["val"].median().sort_values(ascending=False).index.tolist())

    fig, ax = plt.subplots(figsize=(8, 5))
    data = [m.loc[m["trigger_class"] == c, "val"].values for c in order]
    ax.boxplot(data, labels=[f"{c}\n(n={len(v)})" for c, v in zip(order, data)],
               showfliers=False)
    rng = np.random.RandomState(0)
    for i, c in enumerate(order, start=1):
        vals = m.loc[m["trigger_class"] == c, "val"].values
        ax.scatter(np.full(len(vals), i) + rng.uniform(-0.12, 0.12, len(vals)),
                   vals, s=12, c="#7a7a7a", alpha=0.5, edgecolors="none", zorder=2)
    for pid, name in KNOWN_POSITIVES.items():
        r = m[m["pair_id"] == pid]
        if not len(r):
            continue
        c = r["trigger_class"].iloc[0]
        if c not in order:
            continue
        xi = order.index(c) + 1
        yi = float(r["val"].iloc[0])
        ax.scatter([xi], [yi], s=70, marker="*", color="#d62728", zorder=5)
        ax.annotate(name, (xi, yi), fontsize=8, fontweight="bold",
                    xytext=(5, 0), textcoords="offset points")
    ax.set_ylabel("mean cross-method concordance (corrected)")
    ax.set_title("Concordance signal by fold-switch trigger class", fontsize=11)
    ax.tick_params(axis="x", labelsize=8)
    ax.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    out = os.path.join(output_dir, "global_trigger_stratification.png")
    fig.savefig(out, dpi=180); plt.close(fig)
    print(f"[plot] wrote {out}")


def write_global_html(output_html: str = os.path.join("docs", "pairs_global_analysis.html"),
                      title: str = "Global Analysis — Fold-Switch Evolution") -> str:
    """Render the global analysis page from GLOBAL_FIGURE_SPEC (ordered + captioned).

    Image paths in the spec are relative to the page's own directory; spec entries
    whose file is missing are skipped. Pure stdlib + the spec, so it imports
    without the heavy pipeline dependencies.
    """
    import html as _html

    out_dir = os.path.abspath(os.path.dirname(output_html))
    os.makedirs(out_dir, exist_ok=True)

    sections = []
    for spec in GLOBAL_FIGURE_SPEC:
        rel = spec["file"]
        if not os.path.isfile(os.path.join(out_dir, rel)):
            print(f"[html] skip (missing): {rel}")
            continue
        sections.append(
            '<section class="fig">\n'
            f'  <h2>{_html.escape(spec["title"])}</h2>\n'
            f'  <p class="caption">{_html.escape(spec["caption"])}</p>\n'
            f'  <img loading="lazy" src="{_html.escape(rel)}" '
            f'alt="{_html.escape(spec["title"])}"/>\n'
            '</section>')

    body = ("\n".join(sections) if sections
            else "<p>No figures available — run the global plots first.</p>")

    html_doc = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<title>{_html.escape(title)}</title>
<style>
  :root {{ color-scheme: dark; }}
  body {{
    background:#121212; color:#E0E0E0;
    font-family:'Segoe UI',Tahoma,Verdana,sans-serif;
    margin:0; padding:28px 20px; line-height:1.5;
  }}
  h1 {{ margin:0 0 6px; font-size:1.6rem; }}
  .intro {{ max-width:60rem; opacity:.8; margin:0 0 28px; }}
  section.fig {{ max-width:60rem; margin:0 0 40px; }}
  section.fig h2 {{ margin:0 0 6px; font-size:1.2rem; color:#9ecbff; }}
  .caption {{ margin:0 0 12px; opacity:.85; }}
  section.fig img {{ max-width:100%; border:1px solid #333; border-radius:4px; background:#fff; }}
</style>
</head>
<body>
<h1>{_html.escape(title)}</h1>
<p class="intro">Cross-pair summary of the 8-method fold-switch analysis
(AF2, AF3, ESMFold2, Boltz-2, ΔΔG, MSA-Transformer, CCMpred, S4PRED) over the
fold-switching pair set. Each figure is described above it.</p>
{body}
</body>
</html>"""
    with open(output_html, "w", encoding="utf-8") as f:
        f.write(html_doc)
    print(f"[html] wrote: {output_html}")
    return output_html


def make_all_global_plots(output_dir: str | None = None) -> None:
    """Single entry point: generate EVERY global plot into output_dir.

    Per-fold support panel (one scatter per method) + the Moran's phylogenetic-
    signal plots (methods x pairs heatmap, method x method correlation, effect-size
    distribution). The global HTML page (gen_html_for_global_plots) renders these
    in a fixed, captioned order.
    """
    if output_dir is None:
        output_dir = GLOBAL_FIG_DIR
    os.makedirs(output_dir, exist_ok=True)
    for fn in (plot_per_fold_support_panel,
               plot_concordance_volcano,
               plot_raw_vs_corrected,
               plot_trigger_stratification,
               plot_method_clade_heatmap,
               plot_method_correlation_heatmap,
               plot_morans_effectsize_by_method):
        try:
            fn(output_dir)
        except Exception as e:
            print(f"[global] WARN {fn.__name__}: {e}")
