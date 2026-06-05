"""Cross-pair global statistics plots (per-fold support panel + Moran's plots)."""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from config import FIGURE_RES_DIR, TABLES_RES
from utils.method_labels import METHOD_ORDER, disp, order_methods

SURVEY_CSV = os.path.join(TABLES_RES, "fold_diversity_survey.csv")

# ---------------------------------------------------------------------
# Ordered + captioned spec for the global analysis HTML page. ONE place
# defines which global figures appear, in what order, and the explanatory
# text shown before each (Rmarkdown-style). Paths are relative to docs/.
# ---------------------------------------------------------------------
GLOBAL_FIGURE_SPEC = [
    {
        "file": "figures/global/global_per_fold_support_panel.png",
        "title": "Per-fold support, by method",
        "caption":
            "Each panel is one method. For every sequence cluster we plot its "
            "predicted support for fold 1 (x) against fold 2 (y). For AF2, AF3, "
            "ESMFold2 and Boltz-2 the axes are TM-scores to each reference "
            "structure; for ΔΔG, MSA-Transformer, CCMpred and S4PRED they are the "
            "unified per-fold proxy (stability / contact recall / secondary-"
            "structure similarity) mapped onto a common scale. A cluster near the "
            "diagonal reaches both folds (fold-switcher); off-diagonal clusters "
            "prefer one fold. Spread along both axes means the method is resolving "
            "fold-specific signal.",
    },
    {
        "file": "figures/global/global_morans_method_x_pair.png",
        "title": "Phylogenetic clade signal — method × pair",
        "caption":
            "Sequence-divergence-corrected clade signal (Moran's I_norm) for every "
            "method (columns) and fold-switch pair (rows). Red = clusters that "
            "prefer the same fold are phylogenetically grouped (clade-structured); "
            "blue = anti-structured; white ≈ none. Rows are sorted by the mean "
            "signal across all methods, so the most clade-structured pairs sit at "
            "the top. Signal is weak overall — few cells exceed I_norm ≈ 0.2.",
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

def plot_per_fold_support_panel(output_dir: str, agg: str = "max") -> None:
    """2x4 panel: per method, support for fold 1 (x) vs fold 2 (y), per cluster.

    Reads the unified survey (docs/fold_diversity_survey.csv). Every method has a
    per-fold support pair in the TM1/TM2 slots: for AF2/AF3/ESMFold2/Boltz-2 these
    are literal TM-scores to each reference fold; for ΔΔG/CCMpred/MSAT/S4PRED they
    are the unified per-fold proxy (ΔΔG / recall / SS-similarity) mapped onto the
    same axes. A point near (high, low) prefers fold 1, (high, high) reaches both
    (fold-switcher), (low, low) neither.
    """
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
    fig, axes = plt.subplots(nrow, ncol, figsize=(4 * ncol, 3.7 * nrow))
    axes = np.atleast_1d(axes).ravel()

    for ax, m in zip(axes, methods):
        s = df[df["method"] == m]
        x = pd.to_numeric(s[x_col], errors="coerce")
        y = pd.to_numeric(s[y_col], errors="coerce")
        good = x.notna() & y.notna()
        x, y = x[good], y[good]
        ax.scatter(x, y, s=8, alpha=0.35, edgecolors="none")
        if len(x):
            lo = float(min(x.min(), y.min()))
            hi = float(max(x.max(), y.max()))
            pad = 0.04 * (hi - lo or 1.0)
            ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad],
                    color="gray", lw=0.8, ls="--", alpha=0.7)
            ax.set_xlim(lo - pad, hi + pad)
            ax.set_ylim(lo - pad, hi + pad)
        ax.set_title(f"{disp(m)}  (n={int(good.sum())})", fontsize=10)
        ax.set_xlabel("support fold 1")
        ax.set_ylabel("support fold 2")
        ax.grid(True, alpha=0.25)

    for ax in axes[len(methods):]:
        ax.axis("off")

    fig.suptitle(f"Per-cluster support for each fold, by method ({agg} over models)",
                 fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
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
    for m in ("morans_I_norm", "morans_z", "morans_I"):
        if m in d.columns and d[m].notna().any():
            return d, m
    print(f"[global] {fn} has no morans_* columns (old thresholded file?); "
          f"re-run phylo_placement. Skipping Moran's plots.")
    return None, None


def plot_method_clade_heatmap(output_dir: str, labels: str = "corrected") -> None:
    """Methods x pairs heatmap of clade-signal strength (effect size / z).

    Rows (pairs) are sorted by the MEAN signal across methods (descending), so the
    most clade-structured pairs rise to the top independent of any single method.
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
    # Method-agnostic ordering: mean signal across methods, strongest at top.
    w = wide.assign(_mean=wide.mean(axis=1, skipna=True)).sort_values(
        "_mean", ascending=False).drop(columns="_mean")
    M = w.to_numpy(dtype=float)
    vlim = float(np.nanpercentile(np.abs(M), 98)) or 1.0
    fig, ax = plt.subplots(figsize=(6, max(4, 0.16 * len(w))))
    im = ax.imshow(M, cmap="RdBu_r", vmin=-vlim, vmax=vlim, aspect="auto")
    ax.set_xticks(range(len(cols))); ax.set_xticklabels([disp(c) for c in cols],
                                                        rotation=45, ha="right")
    ax.set_yticks(range(len(w))); ax.set_yticklabels(w.index, fontsize=4)
    ax.set_title(f"Clade signal per method x pair ({metric}, {labels})\n"
                 f"rows = pairs, sorted by mean signal across methods",
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
        output_dir = FIGURE_RES_DIR
    os.makedirs(output_dir, exist_ok=True)
    for fn in (plot_per_fold_support_panel,
               plot_method_clade_heatmap,
               plot_method_correlation_heatmap,
               plot_morans_effectsize_by_method):
        try:
            fn(output_dir)
        except Exception as e:
            print(f"[global] WARN {fn.__name__}: {e}")
