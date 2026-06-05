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
               plot_method_clade_heatmap,
               plot_method_correlation_heatmap,
               plot_morans_effectsize_by_method):
        try:
            fn(output_dir)
        except Exception as e:
            print(f"[global] WARN {fn.__name__}: {e}")
