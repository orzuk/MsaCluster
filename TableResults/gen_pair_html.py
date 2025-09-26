#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gen_pair_html.py
Build per-pair HTML pages using paths from config.py.

Images can be embedded inline (base64) or linked/copied into a publish directory
defined by config.TABLES_RES/config.PER_PAIR_PUBLISH_SUBDIR.

Usage examples:
  # Inline-embed (self-contained HTML):
  python TableResults/gen_pair_html.py --pairs 2qqjA_4qdsA --mode inline

  # Copy images to publish dir (docs/HTML/figs/<pair>) and link them:
  python TableResults/gen_pair_html.py --pairs 2qqjA_4qdsA --mode copy

  # Link only (assumes images already copied):
  python TableResults/gen_pair_html.py --pairs ALL --mode link
"""
import os, sys, html, base64, argparse, shutil, re
from pathlib import Path

# --- Locate repo root (script expected under TableResults/) and import config.py ---
REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))
try:
    import config as CFG
except Exception as e:
    print("[gen_pair_html] ERROR: cannot import config.py from repo root:", e, file=sys.stderr)
    sys.exit(2)

# --- Pull paths from config.py ---
# Required:
DATA_DIR = Path(getattr(CFG, "DATA_DIR"))
OUTPUT_PATH_NOTEBOOKS = Path(getattr(CFG, "OUTPUT_PATH_NOTEBOOKS"))
TABLES_RES = Path(getattr(CFG, "TABLES_RES"))

# Optional (with sensible defaults):
_PER_PAIR_PUBLISH_SUBDIR = getattr(CFG, "PER_PAIR_PUBLISH_SUBDIR", os.path.join("HTML","figs"))
PUBLISH_FIGURES_DEFAULT = bool(getattr(CFG, "PUBLISH_FIGURES", False))

# PAIR_DIR_RE can be a compiled regex or a string pattern; default matches "a_b"
_pair_re = getattr(CFG, "PAIR_DIR_RE", r"^.+_.+$")
if hasattr(_pair_re, "match"):
    PAIR_DIR_RE = _pair_re
else:
    PAIR_DIR_RE = re.compile(str(_pair_re))

PUBLISH_BASE = TABLES_RES / Path(_PER_PAIR_PUBLISH_SUBDIR)  # e.g., docs/HTML/figs

# --- Optional per-pair table helper (import from either path) ---
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


class PageBuilder:
    def __init__(self, title: str):
        self.title = title
        self.fig_idx = 0
        self.parts: list[str] = []

    def css(self) -> str:
        return """
<style>
  :root { color-scheme: dark; }
  body { background:#121212; color:#E0E0E0; font-family:'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; margin:20px; }
  h1 { margin: 0 0 4px; }
  h2 { margin: 18px 0 10px; }
  .sub { color:#B0BEC5; margin-bottom:10px; }
  .two-up { display:grid; grid-template-columns: 1fr 1fr; gap:16px; align-items:start; }
  figure { margin: 0 0 16px 0; }
  figcaption { text-align:center; font-style:italic; color:#9aa5ad; margin-top:6px; }
  img { max-width:100%; height:auto; display:block; border:1px solid #333; }
  .warn { color:#ff6f6f; }
  table { width: 100%; border-collapse: collapse; margin-top: 8px; }
  th, td { border: 1px solid #333; padding: 8px 10px; text-align: left; white-space: nowrap; }
  th { background:#b71c1c; color:#fff; position: sticky; top: 0; }
  tr:nth-child(even) { background: #2b2b2b; }
  tr:hover { background: #3a3a3a; }
  a { color:#64B5F6; text-decoration:none; } a:hover { text-decoration:underline; }
</style>
"""

    def header(self, pair_id: str) -> str:
        a, b = _pair_tokens(pair_id)
        return f"""
<h1>Analysis of pair {html.escape(pair_id)}</h1>
<div class="sub"><b>{html.escape(a)}</b> vs <b>{html.escape(b)}</b></div>
"""

    def fig_html(self, img_src: str, caption: str) -> str:
        """img_src is either a data: URI or a relative path (for published mode)."""
        self.fig_idx += 1
        cap = f"Figure {self.fig_idx}. {html.escape(caption)}"
        return f"""<figure><img src="{img_src}"/><figcaption>{cap}</figcaption></figure>"""

    def missing_html(self, label: str) -> str:
        return f"""<div class="warn">[missing] {html.escape(label)}</div>"""

    def two_up_html(self, left_html: str, right_html: str) -> str:
        return f"""<div class="two-up">{left_html}{right_html}</div>"""

    def push(self, html_chunk: str): self.parts.append(html_chunk)

    def build(self) -> str:
        return f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<title>{html.escape(self.title)}</title>{self.css()}</head><body>
{''.join(self.parts)}
</body></html>"""


def _to_img_src(png_path: Path, mode: str, publish_dir_for_pair: Path, html_out_dir: Path) -> str | None:
    """
    Return img src for the chosen mode:
      - inline: data: URI (base64)
      - copy:   copy to publish_dir_for_pair and return relative path from html_out_dir
      - link:   assume already present in publish_dir_for_pair; just return relative path
    """
    if mode == "inline":
        return _data_uri_for_file(png_path)
    publish_dir_for_pair.mkdir(parents=True, exist_ok=True)
    dest_path = publish_dir_for_pair / png_path.name
    if mode == "copy":
        try:
            shutil.copy2(png_path, dest_path)
        except Exception:
            return None
    # relative URL from HTML output dir to dest_path
    rel_src = os.path.relpath(dest_path, start=html_out_dir).replace("\\", "/")
    return rel_src


def _figure_from_patterns(pb: PageBuilder, fig_dir: Path, patterns: list[str], caption: str,
                          mode: str, publish_dir_for_pair: Path, html_out_dir: Path):
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
    pb.push(pb.fig_html(src, caption))


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
        if lsrc: left_html = pb.fig_html(lsrc, left_caption + ("" if lp.name.find("all_clusters") == -1 else " (fallback)"))
    if rp is not None:
        rsrc = _to_img_src(rp, mode, publish_dir_for_pair, html_out_dir)
        if rsrc: right_html = pb.fig_html(rcsrc := rsrc, right_caption + ("" if rp.name.find("all_clusters") == -1 else " (fallback)"))

    pb.push(pb.two_up_html(left_html, right_html))


def render_pair_html(pair_id: str, output_dir: Path, mode: str = "inline") -> Path:
    """
    mode: 'inline' (default), 'copy', or 'link'
    """
    assert mode in ("inline", "copy", "link")
    fig_dir = _fig_dir_for_pair(pair_id)
    output_dir.mkdir(parents=True, exist_ok=True)
    out_html = output_dir / f"{pair_id}.html"

    publish_dir_for_pair = PUBLISH_BASE / pair_id  # e.g., docs/HTML/figs/<pair>
    pb = PageBuilder(title=f"Analysis of {pair_id}")
    pb.push(pb.header(pair_id))

    # 1) Tree + heatmap (exactly once)
    _figure_from_patterns(pb, fig_dir,
                          [f"{pair_id}_phytree_cluster.png", "*phytree*cluster*.png", "*tree*heatmap*.png"],
                          "Phylogenetic tree with per-cluster heatmap.",
                          mode, publish_dir_for_pair, output_dir)

    # 2) Side-by-side: Deep-only cmaps vs Best clusters (fallbacks to 'all')
    _two_figures(pb, fig_dir,
                 left_patterns=[f"{pair_id}_deep_clusters_cmap.png", "*deep*clusters*cmap*.png", "*Deep*cmaps*.png"],
                 left_caption="Deep-MSA contact map panel",
                 right_patterns=[f"{pair_id}_best_clusters_cmap.png", "*best*clusters*cmap*.png", "*best*cmaps*.png"],
                 right_caption="Best clusters contact map panel",
                 fallback_left=[f"{pair_id}_all_clusters_cmap.png", "*all*clusters*cmap*.png"],
                 fallback_right=[f"{pair_id}_all_clusters_cmap.png", "*all*clusters*cmap*.png"],
                 mode=mode, publish_dir_for_pair=publish_dir_for_pair, html_out_dir=output_dir)

    # 3) All clusters mosaic (below)
    _figure_from_patterns(pb, fig_dir,
                          [f"{pair_id}_all_clusters_cmap.png", "*all*clusters*cmap*.png"],
                          "All clusters contact-map mosaic (small multiples).",
                          mode, publish_dir_for_pair, output_dir)

    # 4) Native structures (up to 2), self-alignment, fold1↔fold2 overlay
    shown_native = 0
    for pat in [f"{pair_id}_true_*.png", f"{pair_id}_native_*.png", "*true*.png", "*native*.png"]:
        p = _find_first(fig_dir, [pat])
        if p and p.is_file():
            src = _to_img_src(p, mode, publish_dir_for_pair, output_dir)
            if src:
                shown_native += 1
                pb.push(pb.fig_html(src, f"Native/true structure {shown_native}"))
                if shown_native >= 2: break
    if shown_native == 0:
        pb.push(pb.missing_html("Native/true structures"))

    p_self = _find_first(fig_dir, [f"{pair_id}_fold1_self.png", "*self*align*.png", "*self*.png"])
    if p_self and p_self.is_file():
        src = _to_img_src(p_self, mode, publish_dir_for_pair, output_dir)
        if src: pb.push(pb.fig_html(src, "Sanity check: self-alignment ~ TM≈1"))
    else:
        pb.push(pb.missing_html("Self-alignment (sanity)"))

    p_f12 = _find_first(fig_dir, [f"{pair_id}_3d_aligned.png", "*fold1*to*fold2*.png", "*f1*f2*align*.png", "*overlay*fold*1*2*.png"])
    if p_f12 and p_f12.is_file():
        src = _to_img_src(p_f12, mode, publish_dir_for_pair, output_dir)
        if src: pb.push(pb.fig_html(src, "Fold1 aligned to Fold2 (3D overlay)."))
    else:
        pb.push(pb.missing_html("Fold1↔Fold2 alignment image"))

    # 5) AF/ESM diagnostics
    for tag in ["AF2","AF3","ESM2","ESM3","AF","ESM"]:
        p = _find_first(fig_dir, [f"{pair_id}_fold_pair_scatter_plot_{tag}.png", f"fold_pair_scatter_plot_{tag}.png", f"*{tag}*scatter*plot*.png"])
        if p and p.is_file():
            src = _to_img_src(p, mode, publish_dir_for_pair, output_dir)
            if src: pb.push(pb.fig_html(src, f"{tag} diagnostic scatter vs. truth"))

    # 6) Per-pair cluster table
    if build_pair_cluster_table:
        try:
            df = build_pair_cluster_table(pair_id)
        except Exception as e:
            df = None
            pb.push(f"""<div class="warn">[warn] Could not build per-pair cluster table: {html.escape(str(e))}</div>""")
        if df is not None and not df.empty:
            # Round numerical columns
            for c in df.columns:
                try: df[c] = df[c].astype(float).round(2)
                except Exception: pass
            pb.push("<h2>Cluster metrics table</h2>")
            thead = "<tr>" + "".join(f"<th>{html.escape(str(c))}</th>" for c in df.columns) + "</tr>"
            rows = []
            for _, r in df.iterrows():
                tds = []
                for c in df.columns:
                    v = r[c]
                    s = "-" if v is None or (isinstance(v, float) and (v != v)) else str(v)
                    tds.append(f"<td>{html.escape(s)}</td>")
                rows.append("<tr>" + "".join(tds) + "</tr>")
            pb.push(f"<table>{thead}{''.join(rows)}</table>")
        else:
            pb.push("""<div class="warn">No per-pair cluster table available.</div>""")
    else:
        pb.push("""<div class="warn">build_pair_cluster_table not importable; table omitted.</div>""")

    out_html.write_text(pb.build(), encoding="utf-8")
    return out_html


def _discover_pairs() -> list[str]:
    """Find pair directories under DATA_DIR that match the configured regex and contain output_figs/."""
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

    # Warn if copying/linking while config says not to publish (not fatal)
    if args.mode in ("copy","link") and not PUBLISH_FIGURES_DEFAULT:
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
