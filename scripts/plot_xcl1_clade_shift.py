#!/usr/bin/env python3
"""Publication figure for the clade-shift result of one fold-switching pair
(paper Fig. 5: XCL1, pair 2n54B_2hdmA; F1 = 2n54B alternative fold [blue],
F2 = 2hdmA chemokine fold [red]).

Panels
------
  (a) COARSE cluster-level tree (one leaf per ShallowMsa cluster; the same
      tree scripts/phylo_placement.py tests on, built with
      utils.ancestral_utils.build_coarse_tree) drawn with
      utils.phytree_utils.draw_tree_aligned. The highlighted clade (default:
      the AF2 clade-shift best clade) is coloured by its shift direction
      (blue = toward F1, red = toward F2) and shaded across the panels.
  (b) per-cluster heatmap of the centered fold preference (TMdiff_centered
      from docs/fold_diversity_survey.csv), one column per method. Units
      differ (TM-score / kcal/mol / SSS), so each column is divided by its own
      robust range r (98th percentile of |value|), printed under the column;
      blue = +, toward F1; red = -, toward F2.
  (c) annotation columns: clade membership and mean % identity to the
      query (F1 chain) per cluster, over query match columns
      (scripts/seq_divergence_correction.gapaware_pid), computed from the
      cluster a3m files, falling back to docs/per_cluster_pid_symmetric.csv.
  (d) per method: values of the highlighted clusters vs all other clusters,
      with the clade-shift scan p (clade_p of that method's row in the
      phylo_placement clade CSV) and whether that method's own best clade
      is the highlighted one.

Usage (cluster, from the repo root):
    python3 scripts/plot_xcl1_clade_shift.py \
        --clade-csv /sci/labs/orzuk/orzuk/github/clade_shift_results/final/phylo_clade_raw.csv \
        --out docs/figures/fig5_xcl1_clade_shift
    # any pair's best clade (SI):
    python3 scripts/plot_xcl1_clade_shift.py --pair 2namA_1uxmK \
        --clade-csv <csv> --out docs/figures/si_clade_2namA_1uxmK
    # explicit clade:
    python3 scripts/plot_xcl1_clade_shift.py --clade-csv <csv> \
        --highlight ShallowMsa_005,ShallowMsa_037,ShallowMsa_038,ShallowMsa_075

Writes <out>.png (300 dpi), <out>.pdf and <out>_data.csv (plotted values).
"""

from __future__ import annotations

import os
# Headless rendering env vars must be set before importing matplotlib / ete3.
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("MPLBACKEND", "Agg")

import argparse
import re
import sys
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from utils.method_labels import disp  # noqa: E402
from scripts.seq_divergence_correction import (  # noqa: E402
    normalize_cluster, load_a3m_seqs, gapaware_pid)

# Heatmap / strip-panel methods, left to right (matches the paper's order).
FIG_METHODS = ["AF2", "AF3", "ESM", "Boltz2", "BioEmu", "DDG", "S4PRED"]
# Native unit of TMdiff_centered per method (for the strip-panel titles).
_UNIT = {"DDG": "kcal/mol", "S4PRED": "ΔSSS"}

F1_COLOR = "#1f77b4"   # project colours: blue = fold 1, red = fold 2
F2_COLOR = "#d62728"
REST_COLOR = "#8c8c8c"
FIG_WIDTH_IN = 17.0 / 2.54
FONT_PT = 7


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------

def short_label(tag: str) -> str:
    """ShallowMsa_019 -> S19, DeepMsa -> D (as utils.ancestral_utils)."""
    if str(tag).lower().startswith("deep"):
        return "D"
    m = re.search(r"(\d+)$", str(tag))
    return f"S{int(m.group(1)):02d}" if m else str(tag)


def parse_members(s) -> List[str]:
    """'ShallowMsa_005;ShallowMsa_037' (or comma list / S05) -> tags."""
    if s is None or (isinstance(s, float) and np.isnan(s)):
        return []
    out = []
    for tok in re.split(r"[;,\s]+", str(s).strip()):
        if not tok:
            continue
        m = re.fullmatch(r"S(\d+)", tok)
        out.append(f"ShallowMsa_{int(m.group(1)):03d}" if m
                   else normalize_cluster(tok))
    return out


def load_values(survey_csv: str, pair: str, methods: Sequence[str],
                value_col: str = "TMdiff_centered") -> pd.DataFrame:
    """Per-cluster value matrix (index = cluster tag, columns = methods)."""
    df = pd.read_csv(survey_csv)
    sub = df[df["pair_id"] == pair].copy()
    if sub.empty:
        raise SystemExit(f"[error] pair {pair} not in {survey_csv}")
    sub["_tag"] = sub["cluster"].map(normalize_cluster)
    sub[value_col] = pd.to_numeric(sub[value_col], errors="coerce")
    mat = sub.pivot_table(index="_tag", columns="method", values=value_col,
                          aggfunc="first")
    return mat.reindex(columns=[m for m in methods])


def load_clade_rows(clade_csv: Optional[str], pair: str) -> Dict[str, dict]:
    """{method: row-dict} of the clade-shift CSV rows for ``pair``."""
    if not clade_csv:
        return {}
    if not os.path.isfile(clade_csv):
        print(f"[warn] clade CSV not found: {clade_csv}")
        return {}
    df = pd.read_csv(clade_csv)
    sub = df[df["pair_id"] == pair]
    return {str(r["method"]): r.to_dict() for _, r in sub.iterrows()}


def resolve_highlight(arg: Optional[str], clade_rows: Dict[str, dict],
                      ref_method: str = "AF2") -> List[str]:
    if arg:
        return parse_members(arg)
    row = clade_rows.get(ref_method)
    if row is None:
        raise SystemExit(f"[error] no --highlight given and no {ref_method} row "
                         f"in the clade CSV")
    tags = parse_members(row.get("clade_best_members"))
    if not tags:
        raise SystemExit(f"[error] {ref_method} row has empty clade_best_members "
                         f"(note={row.get('clade_note')!r})")
    return tags


def cluster_pid_from_a3m(pair_dir: str) -> Dict[str, float]:
    """Mean % identity of each cluster's sequences to the query (first
    sequence of output_get_msa/DeepMsa.a3m) over query match columns; same
    computation as scripts/seq_divergence_correction.compute_pair_pids."""
    deep = os.path.join(pair_dir, "output_get_msa", "DeepMsa.a3m")
    msa_dir = os.path.join(pair_dir, "output_msa_cluster")
    if not (os.path.isfile(deep) and os.path.isdir(msa_dir)):
        return {}
    deep_seqs = load_a3m_seqs(deep)
    if not deep_seqs:
        return {}
    query = deep_seqs[0]
    out = {}
    for fn in sorted(os.listdir(msa_dir)):
        m = re.fullmatch(r"ShallowMsa_(\d+)\.a3m", fn)
        if not m:
            continue
        seqs = load_a3m_seqs(os.path.join(msa_dir, fn))
        pids = [gapaware_pid(query, s) for s in seqs[1:]]
        pids = [p for p in pids if not np.isnan(p)]
        if pids:
            out[f"ShallowMsa_{int(m.group(1)):03d}"] = 100.0 * float(np.mean(pids))
    return out


def cluster_pid_from_cache(cache_csv: str, pair: str) -> Dict[str, float]:
    if not (cache_csv and os.path.isfile(cache_csv)):
        return {}
    df = pd.read_csv(cache_csv)
    df = df[df["pair_id"] == pair]
    col = "mean_pid_to_F1" if "mean_pid_to_F1" in df.columns else "mean_pid_to_query"
    return {normalize_cluster(c): 100.0 * float(v)
            for c, v in zip(df["cluster"], df[col]) if pd.notna(v)}


def build_tree(pair: str, tree_nwk: Optional[str] = None):
    """Return (tree, tag_to_label). Default: the coarse tree exactly as
    phylo_placement.py builds it; --tree-nwk loads a saved coarse Newick
    whose leaves are 'S05'-style labels or 'ShallowMsa_005' tags."""
    if tree_nwk:
        from ete3 import Tree
        tree = Tree(tree_nwk, format=1)
        tag_to_label = {}
        for lf in tree.iter_leaves():
            tag = parse_members(lf.name)[0] if lf.name != "D" else "DeepMsa"
            tag_to_label[tag] = short_label(tag)
            lf.name = tag_to_label[tag]
        return tree, tag_to_label
    from utils.ancestral_utils import build_coarse_tree
    tree, tag_to_label, _ = build_coarse_tree(pair)
    return tree, tag_to_label


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def robust_range(v: np.ndarray, q: float = 98.0) -> float:
    a = np.abs(np.asarray(v, dtype=float))
    a = a[np.isfinite(a)]
    if a.size == 0:
        return 1.0
    r = float(np.percentile(a, q))
    if r <= 0:
        r = float(a.max())
    return r if r > 0 else 1.0


def _fmt_p(p) -> str:
    try:
        p = float(p)
    except (TypeError, ValueError):
        return "p = n/a"
    if not np.isfinite(p):
        return "p = n/a"
    return f"p = {p:.1e}" if p < 1e-3 else f"p = {p:.3f}"


def _fmt_r(r: float) -> str:
    return f"±{r:.2g}" if r < 10 else f"±{r:.0f}"


def _runs(mask: Sequence[bool]) -> List[tuple]:
    """Contiguous runs of True as (start, end) inclusive row indices."""
    runs, start = [], None
    for i, v in enumerate(list(mask) + [False]):
        if v and start is None:
            start = i
        elif not v and start is not None:
            runs.append((start, i - 1))
            start = None
    return runs


def _ensure_internal_names(tree) -> None:
    k = 0
    for nd in tree.traverse("preorder"):
        if not nd.is_leaf() and not nd.name:
            nd.name = f"_int{k}"
            k += 1


def _restyle_tree_axis(ax) -> None:
    """Adapt draw_tree_aligned's web styling to the print figure: drop its
    legend, darken neutral branches, shrink tip markers to the row pitch, keep
    only the scale-bar distance and enforce the minimum font size."""
    import matplotlib.colors as mcolors
    leg = ax.get_legend()
    if leg is not None:
        leg.remove()
    neutral = mcolors.to_rgba("#cccccc")
    for ln in ax.lines:
        if ln.get_marker() in ("o", "X"):
            ln.set_markersize(3.2)
            ln.set_markeredgewidth(0.4)
        elif mcolors.to_rgba(ln.get_color()) == neutral:
            ln.set_color("#555555")
        ln.set_linewidth(min(ln.get_linewidth(), 0.8))
    # scale bar: drawn above the tree (axes-fraction y ~1.01) where it would
    # collide with the heatmap's column labels; move it below the tree.
    bar_x = None
    for ln in ax.lines:
        if ln.get_transform() is not ax.transData and len(ln.get_ydata()) == 2:
            ln.set_ydata([-0.012, -0.012])
            ln.set_linewidth(1.0)
            bar_x = float(np.max(ln.get_xdata()))
    for t in ax.texts:
        t.set_text(t.get_text().split("·")[0].strip())
        t.set_fontsize(FONT_PT)
        if t.get_transform() is not ax.transData:
            t.set_y(-0.02)
            t.set_va("top")
            if bar_x is not None:     # right-align: keep clear of panel b
                t.set_x(bar_x)
                t.set_ha("right")


def make_figure(tree, tag_to_label: Dict[str, str], values: pd.DataFrame,
                highlight: Sequence[str], clade_rows: Dict[str, dict],
                pid: Dict[str, float], pair: str, out_base: str,
                mark_label: str = "clade",
                fold_names: Sequence[str] = ("F1", "F2"),
                ref_method: str = "AF2", methods: Sequence[str] = FIG_METHODS,
                height_in: float = 7.2, strips: bool = True) -> pd.DataFrame:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import gridspec
    import matplotlib.patches  # noqa: F401
    from matplotlib.colors import ListedColormap
    from utils.phytree_utils import draw_tree_aligned

    plt.rcParams.update({"font.size": FONT_PT, "axes.labelsize": FONT_PT,
                         "xtick.labelsize": FONT_PT, "ytick.labelsize": FONT_PT,
                         "axes.linewidth": 0.6, "pdf.fonttype": 42,
                         "ps.fonttype": 42})

    _ensure_internal_names(tree)
    label_to_tag = {lab: tag for tag, lab in tag_to_label.items()}
    leaf_order = [lf.name for lf in tree.iter_leaves()]
    tags = [label_to_tag.get(l, l) for l in leaf_order]
    n = len(leaf_order)
    hl = set(highlight)
    missing = hl - set(tags)
    if missing:
        print(f"[warn] highlighted clusters not in tree: {sorted(missing)}")
    in_hl = np.array([t in hl for t in tags])
    methods = [m for m in methods if m in values.columns
               and values[m].notna().any()]
    mat = values.reindex(index=tags, columns=methods).astype(float)

    # Shift direction: the reference method's scan direction, else the sign
    # of the clade-minus-rest difference of the reference (first) method.
    direction = np.nan
    row = clade_rows.get(ref_method, {})
    if set(parse_members(row.get("clade_best_members"))) == hl:
        direction = float(row.get("clade_best_direction", np.nan))
    if not np.isfinite(direction):
        col = mat[ref_method] if ref_method in mat else mat.iloc[:, 0]
        direction = np.sign(np.nanmean(col[in_hl]) - np.nanmean(col[~in_hl]))
    hl_color = F1_COLOR if direction >= 0 else F2_COLOR

    # strips=False drops panel (d), leaving tree, heatmap and annotation only
    if strips:
        fig = plt.figure(figsize=(FIG_WIDTH_IN, height_in))
        outer = gridspec.GridSpec(1, 2, width_ratios=[2.05, 1.0], wspace=0.30,
                                  left=0.01, right=0.985, top=0.90, bottom=0.075)
    else:
        fig = plt.figure(figsize=(FIG_WIDTH_IN * 0.66, height_in))
        outer = gridspec.GridSpec(1, 1, left=0.015, right=0.985, top=0.90,
                                  bottom=0.075)
    left = gridspec.GridSpecFromSubplotSpec(
        2, 5, subplot_spec=outer[0], height_ratios=[1.0, 0.03], hspace=0.17,
        wspace=0.04, width_ratios=[2.3, 0.21 * len(methods), 0.21, 0.21, 0.55])
    ax_tree = fig.add_subplot(left[0, 0])
    ax_heat = fig.add_subplot(left[0, 1], sharey=ax_tree)
    ax_mark = fig.add_subplot(left[0, 2], sharey=ax_tree)
    ax_pid = fig.add_subplot(left[0, 3], sharey=ax_tree)
    ax_lab = fig.add_subplot(left[0, 4], sharey=ax_tree)

    # (a) tree -- shared drawing code from utils.phytree_utils
    clade_nodes = {}
    for nd in tree.traverse("postorder"):
        if nd.is_leaf():
            continue
        lv = [lf.name for lf in nd.iter_leaves()]
        if lv and all(label_to_tag.get(l, l) in hl for l in lv):
            clade_nodes[nd.name] = "F1" if direction >= 0 else "F2"
    draw_tree_aligned(ax_tree, tree, leaf_order,
                      leaf_colors={l: hl_color for l, h in zip(leaf_order, in_hl) if h},
                      internal_node_states=clade_nodes, fold_labels=None)
    _restyle_tree_axis(ax_tree)
    for ln in ax_tree.lines:   # thicker clade branches
        if ln.get_marker() in (None, "None", "") and \
                matplotlib.colors.to_hex(ln.get_color()) == hl_color:
            ln.set_linewidth(1.6)
    ax_tree.set_ylim(n - 0.5, -0.5)

    runs = _runs(in_hl)
    for a in (ax_tree, ax_mark, ax_pid):
        for s, e in runs:
            a.axhspan(s - 0.5, e + 0.5, color=hl_color, alpha=0.13, lw=0,
                      zorder=0)

    # (b) heatmap, each column scaled by its own robust range
    ranges = {m: robust_range(mat[m].values) for m in methods}
    scaled = np.column_stack([np.clip(mat[m].values / ranges[m], -1, 1)
                              for m in methods])
    cmap = plt.get_cmap("RdBu").copy()      # low = red (F2), high = blue (F1)
    cmap.set_bad("#e6e6e6")
    im = ax_heat.imshow(np.ma.masked_invalid(scaled), cmap=cmap, vmin=-1,
                        vmax=1, aspect="auto", interpolation="nearest")
    ax_heat.set_xticks(range(len(methods)))
    ax_heat.set_xticklabels([disp(m) for m in methods], rotation=90)
    ax_heat.xaxis.tick_top()
    ax_heat.tick_params(axis="x", length=0, pad=2)
    ax_heat.tick_params(axis="y", left=False, labelleft=False)
    for j, m in enumerate(methods):
        ax_heat.text(j, n - 0.5 + 0.6, _fmt_r(ranges[m]), rotation=90,
                     ha="center", va="top", fontsize=FONT_PT)
    for s, e in runs:
        ax_heat.add_patch(matplotlib.patches.Rectangle(
            (-0.5, s - 0.5), len(methods), e - s + 1, fill=False,
            edgecolor="black", lw=0.9, zorder=5, clip_on=False))

    # (c) annotation: clade membership + % identity to the query
    ax_mark.imshow(in_hl.astype(float)[:, None], aspect="auto",
                   cmap=ListedColormap(["white", "#222222"]), vmin=0, vmax=1,
                   interpolation="nearest", zorder=1, alpha=None)
    pid_v = np.array([pid.get(t, np.nan) for t in tags], dtype=float)
    pid_im = None
    if np.isfinite(pid_v).any():
        pcmap = plt.get_cmap("Greens").copy()
        pcmap.set_bad("#e6e6e6")
        lo = np.floor(np.nanmin(pid_v) / 10) * 10
        hi = np.ceil(np.nanmax(pid_v) / 10) * 10
        pid_im = ax_pid.imshow(np.ma.masked_invalid(pid_v)[:, None],
                               aspect="auto", cmap=pcmap, vmin=lo,
                               vmax=max(hi, lo + 10), interpolation="nearest")
    for a, lab in ((ax_mark, mark_label), (ax_pid, "% id to query")):
        a.set_xticks([0])
        a.set_xticklabels([lab], rotation=90)
        a.xaxis.tick_top()
        a.tick_params(axis="x", length=0, pad=2)
        a.tick_params(axis="y", left=False, labelleft=False)
        for sp in a.spines.values():
            sp.set_linewidth(0.4)
    for sp in ax_heat.spines.values():
        sp.set_linewidth(0.4)

    # clade labels: one bracket per contiguous run
    ax_lab.axis("off")
    ax_lab.set_xlim(0, 1)
    for s, e in runs:
        ax_lab.plot([0.08, 0.18, 0.18, 0.08], [s - 0.4, s - 0.4, e + 0.4, e + 0.4],
                    color="black", lw=0.7, clip_on=False)
        labs = [short_label(t) for t in tags[s:e + 1]]
        txt = f"{mark_label}\n" + "\n".join(
            " ".join(labs[i:i + 2]) for i in range(0, len(labs), 2))
        # anchor toward the middle of the panel so the text never runs off
        if (s + e) / 2.0 > 0.8 * n:
            yy, va = e + 0.4, "bottom"
        elif (s + e) / 2.0 < 0.2 * n:
            yy, va = s - 0.4, "top"
        else:
            yy, va = (s + e) / 2.0, "center"
        ax_lab.text(0.26, yy, txt, ha="left", va=va, fontsize=FONT_PT,
                    linespacing=1.05)
    if "D" in leaf_order:
        i = leaf_order.index("D")
        ax_lab.text(0.08, i, "Deep MSA", ha="left", va="center",
                    fontsize=FONT_PT)

    # colour bars
    cax = fig.add_subplot(left[1, 1])
    cb = fig.colorbar(im, cax=cax, orientation="horizontal", ticks=[-1, 0, 1])
    cb.ax.set_xticklabels(["−r", "0", "+r"])
    cb.outline.set_linewidth(0.4)
    cb.ax.tick_params(length=2, pad=1)
    f1n, f2n = fold_names
    cb.set_label(f"centered pref. / column scale r\n"
                 f"← {f2n}        {f1n} →", labelpad=1)
    if pid_im is not None:
        pcax = fig.add_subplot(left[1, 3:5])
        pcb = fig.colorbar(pid_im, cax=pcax, orientation="horizontal",
                           ticks=[pid_im.norm.vmin, pid_im.norm.vmax])
        pcb.ax.set_xticklabels([f"{pid_im.norm.vmin:.0f}",
                                f"{pid_im.norm.vmax:.0f}"])
        pcb.outline.set_linewidth(0.4)
        pcb.ax.tick_params(length=2, pad=1)
        pcb.set_label("% id", labelpad=1)

    # (d) per-method clade vs rest
    rs = np.random.RandomState(0)
    n_in, n_out = int(in_hl.sum()), int((~in_hl).sum())
    first_strip = None
    any_other = False
    if strips:
        right = gridspec.GridSpecFromSubplotSpec(len(methods), 1,
                                                 subplot_spec=outer[1], hspace=1.05)
    for k, m in enumerate(methods if strips else []):
        ax = fig.add_subplot(right[k])
        if first_strip is None:
            first_strip = ax
        v = mat[m].values
        vin, vout = v[in_hl & np.isfinite(v)], v[~in_hl & np.isfinite(v)]
        if vout.size:
            ax.boxplot([vout], positions=[0], vert=False, widths=0.55,
                       showfliers=False, patch_artist=True,
                       boxprops=dict(facecolor="#eeeeee", lw=0.5),
                       medianprops=dict(color="black", lw=0.8),
                       whiskerprops=dict(lw=0.5), capprops=dict(lw=0.5))
            ax.scatter(vout, rs.uniform(-0.22, 0.22, vout.size), s=3,
                       color=REST_COLOR, lw=0, zorder=3)
        if vin.size:
            ax.scatter(vin, 1 + rs.uniform(-0.12, 0.12, vin.size), s=11,
                       color=hl_color, edgecolor="black", lw=0.35, zorder=4)
        ax.axvline(0, color="#999999", lw=0.5, ls=":", zorder=1)
        ax.set_ylim(-0.6, 1.6)
        ax.set_yticks([0, 1])
        ax.set_yticklabels([f"other ({n_out})", f"{mark_label} ({n_in})"])
        ax.tick_params(length=2, pad=1, width=0.5)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        crow = clade_rows.get(m)
        if crow is not None:
            same = set(parse_members(crow.get("clade_best_members"))) == hl
            ptxt = _fmt_p(crow.get("clade_p")) + ("" if same else "†")
            any_other = any_other or not same
        else:
            ptxt = "no scan"
        ax.set_title(f"{disp(m)} ({_UNIT.get(m, 'ΔTM')})   {ptxt}",
                     loc="left", fontsize=FONT_PT, pad=2)

    # panel letters
    for a, letter in ((ax_tree, "a"), (ax_heat, "b"), (ax_mark, "c"),
                      (first_strip, "d")):
        if a is None:
            continue
        bb = a.get_position()
        fig.text(bb.x0 - (0.0 if letter == "a" else 0.008), 0.975, letter,
                 fontsize=9, fontweight="bold", ha="left", va="top")
    if any_other:
        fig.text(0.985, 0.012, f"† clade-shift p of that method's own best "
                 f"clade (not the {mark_label} clade)", ha="right",
                 va="bottom", fontsize=FONT_PT)
    comb = clade_rows.get("COMBINED")
    if comb is not None:
        fig.text(0.985, 0.975, f"combined scan: {_fmt_p(comb.get('clade_p'))}",
                 ha="right", va="top", fontsize=FONT_PT)

    os.makedirs(os.path.dirname(os.path.abspath(out_base)), exist_ok=True)
    fig.savefig(out_base + ".png", dpi=300)
    fig.savefig(out_base + ".pdf")
    plt.close(fig)

    data = mat.copy()
    data.insert(0, "leaf", leaf_order)
    data["in_clade"] = in_hl
    data["pid_to_query"] = pid_v
    data.index.name = "cluster"
    data.to_csv(out_base + "_data.csv")
    print(f"[ok] wrote {out_base}.png / .pdf / _data.csv "
          f"({n} leaves, clade = {sorted(hl)})")
    return data


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None):
    from config import DATA_DIR, TABLES_RES
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--pair", default="2n54B_2hdmA")
    p.add_argument("--survey-csv",
                   default=os.path.join(TABLES_RES, "fold_diversity_survey.csv"))
    p.add_argument("--value-col", default="TMdiff_centered")
    p.add_argument("--clade-csv", default=None,
                   help="phylo_placement.py output with clade_* columns")
    p.add_argument("--highlight", default=None,
                   help="comma list of clusters (ShallowMsa_005 or S05); "
                        "default: clade_best_members of the --ref-method row")
    p.add_argument("--ref-method", default="AF2")
    p.add_argument("--methods", default=",".join(FIG_METHODS))
    p.add_argument("--mark-label", default=None,
                   help="name of the highlighted group (default XCL for "
                        "2n54B_2hdmA, else 'clade')")
    p.add_argument("--tree-nwk", default=None,
                   help="optional saved coarse-tree Newick (default: rebuild "
                        "with utils.ancestral_utils.build_coarse_tree)")
    p.add_argument("--pid-cache",
                   default=os.path.join(TABLES_RES, "per_cluster_pid_symmetric.csv"),
                   help="fallback per-cluster PID table if a3m files are absent")
    p.add_argument("--height-cm", type=float, default=18.3)
    p.add_argument("--no-strips", action="store_true",
                   help="omit panel (d), the per-method clade-vs-rest strips")
    p.add_argument("--out", default=None,
                   help="output base path (no extension); default "
                        "docs/figures/clade_shift_<pair>")
    a = p.parse_args(argv)

    out = a.out or os.path.join(TABLES_RES, "figures", f"clade_shift_{a.pair}")
    out = re.sub(r"\.(png|pdf)$", "", out)
    methods = [m.strip() for m in a.methods.split(",") if m.strip()]
    clade_rows = load_clade_rows(a.clade_csv, a.pair)
    highlight = resolve_highlight(a.highlight, clade_rows, a.ref_method)
    values = load_values(a.survey_csv, a.pair, methods, a.value_col)
    tree, tag_to_label = build_tree(a.pair, a.tree_nwk)
    pid = cluster_pid_from_a3m(os.path.join(DATA_DIR, a.pair))
    if not pid:
        print("[info] cluster a3m files not found; using PID cache")
        pid = cluster_pid_from_cache(a.pid_cache, a.pair)
    f1, f2 = (a.pair.split("_") + ["F2"])[:2]
    mark = a.mark_label or ("XCL" if a.pair == "2n54B_2hdmA" else "clade")
    make_figure(tree, tag_to_label, values, highlight, clade_rows, pid, a.pair,
                out, mark_label=mark, fold_names=(f"F1 {f1}", f"F2 {f2}"),
                ref_method=a.ref_method, methods=methods,
                height_in=a.height_cm / 2.54, strips=not a.no_strips)


if __name__ == "__main__":
    main()
