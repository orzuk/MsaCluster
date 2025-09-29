#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gen_pair_html.py — config-aware per‑pair HTML renderer

Fixes requested:
1) Smaller phylogenetic tree figure (max-height)
2) Deep vs Best CMAPs shown side‑by‑side with matched visual dimensions
3) "All clusters" now shows both the mosaic (if present) and a grid of per‑cluster CMAPs (if available)
4) 3D structure panels: attempt to show all requested overlays (AF2/AF3/ESM2/ESM3, per fold). If a figure is missing, emit a bold callout
5) Narrower table layout (max width, smaller font, wrapping cells)
6) Order: 3D structures first, then CMAPs, then table, finally tree

Modes:
  --mode inline  : embed images as base64 (self‑contained HTML; best for GitHub Pages)
  --mode copy    : copy figures to publish dir (config.TABLES_RES/config.PER_PAIR_PUBLISH_SUBDIR) and link them
  --mode link    : just link to already‑published figures

Usage:
  python TableResults/gen_pair_html.py --pairs 2qqjA_4qdsA --mode inline
  python TableResults/gen_pair_html.py --pairs ALL --mode copy
"""
import os, sys, html, base64, argparse, shutil, re
from pathlib import Path

# --- Import config from repo root ---
REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))
try:
    import config as CFG
except Exception as e:
    print("[gen_pair_html] ERROR: cannot import config.py from repo root:", e, file=sys.stderr)
    sys.exit(2)

# --- Paths from config.py ---
DATA_DIR = Path(getattr(CFG, "DATA_DIR"))
OUTPUT_PATH_NOTEBOOKS = Path(getattr(CFG, "OUTPUT_PATH_NOTEBOOKS"))
TABLES_RES = Path(getattr(CFG, "TABLES_RES"))
PAIR_DIR_RE = getattr(CFG, "PAIR_DIR_RE", re.compile(r"^.+_.+$"))
if not hasattr(PAIR_DIR_RE, "match"):
    PAIR_DIR_RE = re.compile(str(PAIR_DIR_RE))

PER_PAIR_PUBLISH_SUBDIR = getattr(CFG, "PER_PAIR_PUBLISH_SUBDIR", os.path.join("HTML","figs"))
PUBLISH_FIGURES = bool(getattr(CFG, "PUBLISH_FIGURES", True))
PUBLISH_BASE = TABLES_RES / Path(PER_PAIR_PUBLISH_SUBDIR)  # e.g., docs/HTML/figs

# Optional knobs
MAX_TREE_HEIGHT = 520    # px
MAX_TWOCOL_HEIGHT = 520  # px, deep vs best height cap

# Optional cluster CMAP filename patterns (loose, to pick any naming)
CLUSTER_TILE_PATTERNS = [
    "*cluster*_*cmap*.png",
    "*cluster*_*CMAP*.png",
    "*cluster*_*contact*.png",
]

# 3D structure figure patterns by label
STRUCT_PATTERNS = {
    "two_structures": [
        "*two*structures*.png", "*two_structures*.png", "*2struct*.png",
        "*true*vs*true*.png", "*native*vs*native*.png"
    ],
    "two_structures_aligned": [
        "*aligned_two_structures*.png", "*two*structures*aligned*.png",
        "*overlay*true*true*.png", "*overlay*native*native*.png", "*3d_aligned*.png"
    ],
    # Fold1/Fold2 vs best model variants
    "fold1_vs_best_AF2": ["*fold1*best*AF2*.png", "*AF2*best*fold1*.png", "*fold1*AF2*best*.png"],
    "fold2_vs_best_AF2": ["*fold2*best*AF2*.png", "*AF2*best*fold2*.png", "*fold2*AF2*best*.png"],
    "fold1_vs_best_AF3": ["*fold1*best*AF3*.png", "*AF3*best*fold1*.png", "*fold1*AF3*best*.png"],
    "fold2_vs_best_AF3": ["*fold2*best*AF3*.png", "*AF3*best*fold2*.png", "*fold2*AF3*best*.png"],
    "fold1_vs_best_ESM2": ["*fold1*best*ESM2*.png", "*ESM2*best*fold1*.png", "*fold1*ESM2*best*.png"],
    "fold2_vs_best_ESM2": ["*fold2*best*ESM2*.png", "*ESM2*best*fold2*.png", "*fold2*ESM2*best*.png"],
    "fold1_vs_best_ESM3": ["*fold1*best*ESM3*.png", "*ESM3*best*fold1*.png", "*fold1*ESM3*best*.png"],
    "fold2_vs_best_ESM3": ["*fold2*best*ESM3*.png", "*ESM3*best*fold2*.png", "*fold2*ESM3*best*.png"],
}

# Interactive 3D viewer patterns (.html saved by the plotting code)
HTML_STRUCT_PATTERNS = {
    "two_structures_aligned_html": ["*two*structures*aligned*.html", "*aligned*true*true*.html", "*aligned*.html"],
    "fold1_vs_best_AF2_html": ["*fold1*best*AF2*.html"],
    "fold2_vs_best_AF2_html": ["*fold2*best*AF2*.html"],
    "fold1_vs_best_AF3_html": ["*fold1*best*AF3*.html"],
    "fold2_vs_best_AF3_html": ["*fold2*best*AF3*.html"],
    "fold1_vs_best_ESM2_html": ["*fold1*best*ESM2*.html"],
    "fold2_vs_best_ESM2_html": ["*fold2*best*ESM2*.html"],
    "fold1_vs_best_ESM3_html": ["*fold1*best*ESM3*.html"],
    "fold2_vs_best_ESM3_html": ["*fold2*best*ESM3*.html"],
}



# --- Optional per-pair table helper ---
build_pair_cluster_table = None
try:
    from Analysis.postprocess_unified import build_pair_cluster_table  # type: ignore
except Exception:
    try:
        from postprocess_unified import build_pair_cluster_table  # type: ignore
    except Exception:
        build_pair_cluster_table = None



def _pair_tokens(pair_id: str) -> tuple[str, str]:
    if "_" in pair_id:
        a, b = pair_id.split("_", 1)
    else:
        a, b = pair_id, ""
    return a, b


def _fig_dir_for_pair(pair_id: str) -> Path:
    # DATA_DIR/<pair>/output_figs
    return Path(DATA_DIR) / pair_id / "output_figs"


def _find_first(fig_dir: Path, patterns: list[str]) -> Path | None:
    import glob
    for pat in patterns:
        cands = sorted(glob.glob(str(fig_dir / pat)))
        if cands:
            return Path(cands[0])
    return None


def _find_all(fig_dir: Path, patterns: list[str]) -> list[Path]:
    import glob
    seen = []
    for pat in patterns:
        for f in sorted(glob.glob(str(fig_dir / pat))):
            p = Path(f)
            if p.is_file() and p not in seen:
                seen.append(p)
    return seen


def _data_uri_for_file(p: Path) -> str | None:
    try:
        data = p.read_bytes()
        b64 = base64.b64encode(data).decode("ascii")
        if p.suffix.lower() in (".jpg", ".jpeg"):
            mime = "image/jpeg"
        elif p.suffix.lower() == ".svg":
            mime = "image/svg+xml"
        else:
            mime = "image/png"
        return f"data:{mime};base64,{b64}"
    except Exception:
        return None


def _data_uri_for_html_file(p: Path) -> str | None:
    try:
        b64 = base64.b64encode(p.read_bytes()).decode("ascii")
        return f"data:text/html;base64,{b64}"
    except Exception:
        return None



def _to_html_src(html_path: Path, mode: str, publish_dir_for_pair: Path, html_out_dir: Path) -> str | None:
    # Force file-based src for interactivity (data: URIs can disable JS)
    publish_dir_for_pair.mkdir(parents=True, exist_ok=True)
    dest_path = publish_dir_for_pair / html_path.name

    if mode in ("inline", "copy"):
        try:
            shutil.copy2(html_path, dest_path)
        except Exception:
            return None
    # link mode just points to existing published file
    rel_src = os.path.relpath(dest_path, start=html_out_dir).replace("\\", "/")
    return rel_src


class PageBuilder:
    def __init__(self, title: str):
        self.title = title
        self.fig_idx = 0
        self.parts: list[str] = []

    def css(self) -> str:
        return f"""
<style>
  :root {{ color-scheme: dark; }}
  body {{ background:#121212; color:#E0E0E0; font-family:'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; margin:20px; }}
  h1 {{ margin: 0 0 4px; }}
  h2 {{ margin: 18px 0 10px; }}
  .sub {{ color:#B0BEC5; margin-bottom:10px; }}

  /* Two-up panels: enforce a shared max-height so Deep/Best look the same size */
  .two-up {{ display:grid; grid-template-columns: 1fr 1fr; gap:16px; align-items:start; }}
  .two-up figure, .two-up .imgbox {{ height:{MAX_TWOCOL_HEIGHT}px; display:flex; align-items:center; justify-content:center; overflow:hidden; background:transparent; margin:0; }}
  .two-up img {{ max-height:{MAX_TWOCOL_HEIGHT}px; width:auto; height:auto; object-fit:contain; display:block; border:1px solid #333; }}

  /* Tree: cap its height and center it, keep width auto */
  .tree figure, .tree .imgbox {{ max-height:{MAX_TREE_HEIGHT}px; display:flex; align-items:center; justify-content:center; overflow:hidden; margin: 0 auto 16px auto; }}
  .tree img {{ max-height:{MAX_TREE_HEIGHT}px; width:auto; height:auto; border:1px solid #333; display:block; }}

  figure {{ margin: 0 0 12px 0; }}
  figcaption {{ text-align:center; font-style:italic; color:#9aa5ad; margin-top:6px; }}
  img {{ max-width:100%; height:auto; display:block; border:1px solid #333; }}

  .warn {{ color:#ff6f6f; }}
  .callout-missing {{ color:#ff8080; font-weight:700; font-size:1.05em; margin:10px 0; }}

  /* Table container: limit width so it doesn't span full page */
  .table-wrap {{ max-width: 960px; margin: 10px auto 18px auto; }}
  table {{ width: 100%; border-collapse: collapse; font-size: 0.95em; }}
  th, td {{ border: 1px solid #333; padding: 6px 8px; text-align: left; white-space: normal; }}
  th {{ background:#b71c1c; color:#fff; position: sticky; top: 0; }}
  tr:nth-child(even) {{ background: #2b2b2b; }}
  tr:hover {{ background: #3a3a3a; }}

 .iframebox {{ width: 100%; height: 820px; border:1px solid #333; background:#fff; }}
  iframe.viewer {{ width: 100%; height: 100%; border: 0; }}

  /* Grid of per-cluster thumbnails */
  .thumb-grid {{ display:grid; grid-template-columns: repeat(auto-fill, minmax(220px, 1fr)); gap:12px; }}
  .thumb-grid figure {{ margin:0; }}
  .thumb-grid img {{ width:100%; height:auto; border:1px solid #333; }}
</style>
"""

    def header(self, pair_id: str) -> str:
        a, b = _pair_tokens(pair_id)
        return f"""
<h1>Analysis of pair {html.escape(pair_id)}</h1>
<div class="sub"><b>{html.escape(a)}</b> vs <b>{html.escape(b)}</b></div>
"""

    def fig_html(self, img_src: str, caption: str, extra_class: str = "") -> str:
        self.fig_idx += 1
        cap = f"Figure {self.fig_idx}. {html.escape(caption)}"
        cls = f' class="{extra_class}"' if extra_class else ""
        return f"""<div{cls}><figure class="imgbox"><img src="{img_src}"/></figure><figcaption>{cap}</figcaption></div>"""

    def missing_html(self, label: str) -> str:
        return f"""<div class="warn">[missing] {html.escape(label)}</div>"""

    def callout_missing(self, label: str) -> str:
        return f"""<div class="callout-missing">{html.escape(label)}</div>"""

    def two_up_html(self, left_html: str, right_html: str) -> str:
        return f"""<div class="two-up">{left_html}{right_html}</div>"""

    def push(self, html_chunk: str): self.parts.append(html_chunk)

    def build(self) -> str:
        return f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<title>{html.escape(self.title)}</title>{self.css()}</head><body>
{''.join(self.parts)}
</body></html>"""

    def fig_iframe(self, iframe_src: str, caption: str) -> str:
        self.fig_idx += 1
        cap = f"Figure {self.fig_idx}. {html.escape(caption)}"
        return f"""<div><div class="iframebox">
<iframe class="viewer" src="{iframe_src}"></iframe>
</div><figcaption>{cap}</figcaption></div>"""


def _iframe_from_patterns(pb: PageBuilder, fig_dir: Path, patterns: list[str], caption: str,
                          mode: str, publish_dir_for_pair: Path, html_out_dir: Path):
    p = _find_first(fig_dir, patterns)
    if not p or not p.is_file():
        pb.push(pb.missing_html(caption)); return
    src = _to_html_src(p, mode, publish_dir_for_pair, html_out_dir)
    if src is None:
        pb.push(pb.missing_html(f"{caption} (unreadable: {p.name})")); return
    pb.push(pb.fig_iframe(src, caption))



def _to_img_src(png_path: Path, mode: str, publish_dir_for_pair: Path, html_out_dir: Path) -> str | None:
    if mode == "inline":
        return _data_uri_for_file(png_path)
    publish_dir_for_pair.mkdir(parents=True, exist_ok=True)
    dest_path = publish_dir_for_pair / png_path.name
    if mode == "copy":
        try:
            shutil.copy2(png_path, dest_path)
        except Exception:
            return None
    rel_src = os.path.relpath(dest_path, start=html_out_dir).replace("\\", "/")
    return rel_src


def _figure_from_patterns(pb: PageBuilder, fig_dir: Path, patterns: list[str], caption: str,
                          mode: str, publish_dir_for_pair: Path, html_out_dir: Path,
                          extra_class: str = ""):
    import glob
    png_path = None
    for pat in patterns:
        cands = sorted(glob.glob(str(fig_dir / pat)))
        if cands:
            png_path = Path(cands[0])
            break
    if not png_path or not png_path.is_file():
        pb.push(pb.missing_html(caption)); return
    src = _to_img_src(png_path, mode, publish_dir_for_pair, html_out_dir)
    if src is None:
        pb.push(pb.missing_html(f"{caption} (unreadable: {png_path.name})")); return
    pb.push(pb.fig_html(src, caption, extra_class=extra_class))


def _two_figures(pb: PageBuilder, fig_dir: Path, left_patterns: list[str], left_caption: str,
                 right_patterns: list[str], right_caption: str, fallback_right: list[str] | None,
                 fallback_left: list[str] | None, mode: str, publish_dir_for_pair: Path, html_out_dir: Path):
    import glob
    def find_one(pats):
        for pat in (pats or []):
            cands = sorted(glob.glob(str(fig_dir / pat)))
            if cands: return Path(cands[0])
        return None

    lp = find_one(left_patterns) or find_one(fallback_left)
    rp = find_one(right_patterns) or find_one(fallback_right)

    if lp is None and rp is None:
        pb.push(pb.missing_html(left_caption + " / " + right_caption)); return

    left_html = pb.missing_html(left_caption)
    right_html = pb.missing_html(right_caption)
    if lp is not None:
        lsrc = _to_img_src(lp, mode, publish_dir_for_pair, html_out_dir)
        if lsrc: left_html = pb.fig_html(lsrc, left_caption)
    if rp is not None:
        rsrc = _to_img_src(rp, mode, publish_dir_for_pair, html_out_dir)
        if rsrc: right_html = pb.fig_html(rsrc, right_caption)

    pb.push(pb.two_up_html(left_html, right_html))


def _thumb_grid(pb: PageBuilder, fig_dir: Path, patterns: list[str], mode: str,
                publish_dir_for_pair: Path, html_out_dir: Path, title: str):
    thumbs = _find_all(fig_dir, patterns)
    if not thumbs:
        return False
    pb.push(f"<h2>{html.escape(title)}</h2>")
    items = []
    for p in thumbs:
        name = p.name.lower()
        # Skip panels that are already shown elsewhere
        if any(k in name for k in ("best", "deep", "all")):
            continue
        src = _to_img_src(p, mode, publish_dir_for_pair, html_out_dir)
        if not src:
            continue
        items.append(f'<figure><img src="{src}" alt="{html.escape(p.name)}"></figure>')
    if not items:
        return False
    pb.push(f'<div class="thumb-grid">{"".join(items)}</div>')
    return True


def render_pair_html(pair_id: str, output_dir: Path, mode: str = "inline") -> Path:
    assert mode in ("inline", "copy", "link")
    fig_dir = _fig_dir_for_pair(pair_id)
    output_dir.mkdir(parents=True, exist_ok=True)
    out_html = output_dir / f"{pair_id}.html"

    publish_dir_for_pair = PUBLISH_BASE / pair_id  # e.g., docs/HTML/figs/<pair>
    pb = PageBuilder(title=f"Analysis of {pair_id}")
    pb.push(pb.header(pair_id))

    # ===== ORDER 1: 3D STRUCTURES FIRST =====
    # Two-structures (un-aligned) and aligned
    _figure_from_patterns(pb, fig_dir, STRUCT_PATTERNS["two_structures"], "Two structures (Fold1 & Fold2)",
                          mode, publish_dir_for_pair, output_dir)
    _figure_from_patterns(pb, fig_dir, STRUCT_PATTERNS["two_structures_aligned"], "Two structures aligned",
                          mode, publish_dir_for_pair, output_dir)

    # Fold1/Fold2 vs best models (AF2/AF3/ESM2/ESM3)
    for tag, nice in [
        ("fold1_vs_best_AF2", "Fold1 vs best AF2"),
        ("fold2_vs_best_AF2", "Fold2 vs best AF2"),
        ("fold1_vs_best_AF3", "Fold1 vs best AF3"),
        ("fold2_vs_best_AF3", "Fold2 vs best AF3"),
        ("fold1_vs_best_ESM2", "Fold1 vs best ESM2"),
        ("fold2_vs_best_ESM2", "Fold2 vs best ESM2"),
        ("fold1_vs_best_ESM3", "Fold1 vs best ESM3"),
        ("fold2_vs_best_ESM3", "Fold2 vs best ESM3"),
    ]:
        p = _find_first(fig_dir, STRUCT_PATTERNS[tag])
        # Skip static if interactive version exists
        html_tag = f"{tag}_html"
        if html_tag in HTML_STRUCT_PATTERNS:
            interp = _find_first(fig_dir, HTML_STRUCT_PATTERNS[html_tag])
            if interp and interp.is_file():
                continue
        if p and p.is_file():
            src = _to_img_src(p, mode, publish_dir_for_pair, output_dir)
            if src: pb.push(pb.fig_html(src, nice))
        else:
            # Bold callout for missing items (as requested)
            pb.push(pb.callout_missing(f"Missing: {nice}"))

    # 3D INTERACTIVE VIEWERS (if present)
    _iframe_from_patterns(pb, fig_dir, HTML_STRUCT_PATTERNS["two_structures_aligned_html"],
                          "Two structures aligned — interactive", mode, publish_dir_for_pair, output_dir)

    for tag, nice in [
        ("fold1_vs_best_AF2_html", "Fold1 vs best AF2 — interactive"),
        ("fold2_vs_best_AF2_html", "Fold2 vs best AF2 — interactive"),
        ("fold1_vs_best_AF3_html", "Fold1 vs best AF3 — interactive"),
        ("fold2_vs_best_AF3_html", "Fold2 vs best AF3 — interactive"),
        ("fold1_vs_best_ESM2_html", "Fold1 vs best ESM2 — interactive"),
        ("fold2_vs_best_ESM2_html", "Fold2 vs best ESM2 — interactive"),
        ("fold1_vs_best_ESM3_html", "Fold1 vs best ESM3 — interactive"),
        ("fold2_vs_best_ESM3_html", "Fold2 vs best ESM3 — interactive"),
    ]:
        _iframe_from_patterns(pb, fig_dir, HTML_STRUCT_PATTERNS[tag], nice, mode, publish_dir_for_pair, output_dir)

    # ===== ORDER 2: CMAPS (Deep vs Best side-by-side, then ALL) =====
    _two_figures(pb, fig_dir,
                 left_patterns=[f"{pair_id}_deep_cmap.png", f"{pair_id}_deep_clusters_cmap.png",
                                "*deep*_cmap*.png", "*deep*clusters*cmap*.png", "*Deep*cmaps*.png"],
                 left_caption="Deep-MSA contact map panel",
                 right_patterns=[f"{pair_id}_best_clusters_cmap.png", "*best*clusters*cmap*.png", "*best*cmaps*.png"],
                 right_caption="Best clusters contact map panel",
                 fallback_left=[f"{pair_id}_all_clusters_cmap.png", "*all*clusters*cmap*.png"],
                 fallback_right=[f"{pair_id}_all_clusters_cmap.png", "*all*clusters*cmap*.png"],
                 mode=mode, publish_dir_for_pair=publish_dir_for_pair, html_out_dir=output_dir)

    # All clusters mosaic if present
    _figure_from_patterns(pb, fig_dir,
                          [f"{pair_id}_all_clusters_cmap.png", "*all*clusters*cmap*.png"],
                          "All clusters contact-map mosaic (small multiples).",
                          mode, publish_dir_for_pair, output_dir)

    # Grid of per-cluster tiles (if saved as many small images)
    _thumb_grid(pb, fig_dir, CLUSTER_TILE_PATTERNS, mode, publish_dir_for_pair, output_dir,
                title="Per-cluster contact maps")

    # ===== ORDER 3: TABLE (narrow container) =====
    if build_pair_cluster_table:
        try:
            df = build_pair_cluster_table(pair_id)
        except Exception as e:
            df = None
            pb.push(f"""<div class="warn">[warn] Could not build per-pair cluster table: {html.escape(str(e))}</div>""")
        if df is not None and not df.empty:
            for c in df.columns:
                try: df[c] = df[c].astype(float).round(2)
                except Exception: pass
            thead = "<tr>" + "".join(f"<th>{html.escape(str(c))}</th>" for c in df.columns) + "</tr>"
            rows = []
            for _, r in df.iterrows():
                tds = []
                for c in df.columns:
                    v = r[c]
                    s = "-" if v is None or (isinstance(v, float) and (v != v)) else str(v)
                    tds.append(f"<td>{html.escape(s)}</td>")
                rows.append("<tr>" + "".join(tds) + "</tr>")
            pb.push(f'<div class="table-wrap"><h2>Cluster metrics table</h2><table>{thead}{"".join(rows)}</table></div>')
        else:
            pb.push('<div class="warn">No per-pair cluster table available.</div>')
    else:
        pb.push('<div class="warn">build_pair_cluster_table not importable; table omitted.</div>')

    # ===== ORDER 4: TREE LAST (smaller) =====
    _figure_from_patterns(pb, fig_dir,
                          [f"{pair_id}_phytree_cluster.png", "*phytree*cluster*.png", "*tree*heatmap*.png"],
                          "Phylogenetic tree with per-cluster heatmap.",
                          mode, publish_dir_for_pair, output_dir,
                          extra_class="tree")

    out_html.write_text(pb.build(), encoding="utf-8")
    return out_html


def _discover_pairs() -> list[str]:
    pairs = []
    pipedir = Path(DATA_DIR)
    if not pipedir.is_dir():
        return pairs
    for d in pipedir.iterdir():
        if not d.is_dir(): continue
        if not (d / "output_figs").is_dir(): continue
        if PAIR_DIR_RE.match(d.name):
            pairs.append(d.name)
    return sorted(pairs)


def main():
    ap = argparse.ArgumentParser(description="Generate per-pair HTML pages using config paths.")
    ap.add_argument("--pairs", type=str, required=True,
                    help="Comma-separated pair IDs (e.g., 2qqjA_4qdsA,1abcX_2defY) or 'ALL'.")
    ap.add_argument("--mode", choices=["inline","copy","link"], default="inline",
                    help="How to include images in HTML. 'inline' embeds base64; 'copy' copies to docs/HTML/figs/<pair>; 'link' assumes already copied.")
    ap.add_argument("--output_dir", type=str, default=str(OUTPUT_PATH_NOTEBOOKS),
                    help="Directory to write HTML pages (default: config.OUTPUT_PATH_NOTEBOOKS).")
    args = ap.parse_args()

    output_dir = Path(args.output_dir).resolve()

    if args.pairs.strip().upper() == "ALL":
        pairs = _discover_pairs()
    else:
        pairs = [p.strip() for p in args.pairs.split(",") if p.strip()]

    if not pairs:
        print("[gen_pair_html] No pairs to render.", file=sys.stderr); sys.exit(1)

    if args.mode in ("copy","link") and not PUBLISH_FIGURES:
        print("[gen_pair_html] WARNING: config.PUBLISH_FIGURES is False; proceeding anyway.", file=sys.stderr)

    output_dir.mkdir(parents=True, exist_ok=True)

    ok, bad = [], []
    for pair_id in pairs:
        try:
            out = render_pair_html(pair_id, output_dir, mode=args.mode)
            print(f"[gen_pair_html] wrote: {out}")
            ok.append(pair_id)
        except Exception as e:
            print(f"[gen_pair_html] ERROR [{pair_id}]: {e}", file=sys.stderr)
            bad.append(pair_id)

    if bad:
        print(f"[gen_pair_html] Completed with {len(bad)} errors: {', '.join(bad)}", file=sys.stderr)
        sys.exit(2)

if __name__ == "__main__":
    main()
