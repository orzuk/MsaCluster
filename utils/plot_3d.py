"""3D structure rendering: PyMOL (API/CLI) and py3Dmol interactive HTML."""

from PIL import Image, ImageDraw, ImageFont

import pandas as pd
import numpy as np
import os, re, shutil
import subprocess, shlex, tempfile

PYMOL_AVAILABLE = False
try:
    import pymol
    from pymol import cmd
    PYMOL_AVAILABLE = True
except Exception:
    PYMOL_AVAILABLE = False

from config import PYMOL_EXE

from utils.energy_utils import read_energy_tuples, align_and_compare_residues


# ---------------------------------------------------------------------------
# PyMOL binary detection
# ---------------------------------------------------------------------------
PYMOL_BIN = os.environ.get("PYMOL_BIN")


def _detect_pymol_bin() -> str | None:
    # 1) explicit env
    cand = os.environ.get("PYMOL_BIN")
    if cand and os.path.isfile(cand):
        return cand
    # 2) known per-user install (adjust if you move it)
    known = PYMOL_EXE
    if os.path.isfile(known):
        return known
    # 3) in PATH
    for name in ("pymol", "pymol-open-source"):
        path = shutil.which(name)
        if path:
            return path
    return None


# cache at import time, but functions will re-check env so you can override per-run
_PYMOL_BIN_CACHED = _detect_pymol_bin()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _normalize_cluster_tag(s: str) -> str:
    """Canonical cluster tag: 'ShallowMsa_019', 'DeepMsa', etc."""
    s = s.strip()
    if s.lower() in ("deep", "deepmsa", "msa_deep", "deep_msa"):
        return "DeepMsa"
    m = re.fullmatch(r"(\d+)", s)
    if m:
        return f"ShallowMsa_{int(m.group(1)):03d}"
    m2 = re.fullmatch(r"ShallowMsa_(\d+)", s)
    if m2:
        return f"ShallowMsa_{int(m2.group(1)):03d}"
    return s


def _cluster_short_disp(tag: str) -> str:
    """ShallowMsa_019 -> C#19 ; DeepMsa/MSA_deep -> D; else pass-through."""
    t = _normalize_cluster_tag(tag or "")
    if "deep" in t.lower():
        return "D"
    m = re.search(r"(\d+)$", t)
    return f"C#{int(m.group(1))}" if m else t


def _add_png_legend(img_path: str, left_label: str, right_label: str, *,
                    title: str | None = None,
                    left_rgb=(255, 0, 0), right_rgb=(0, 80, 255)) -> None:
    """Append a white strip with a red/blue legend (and optional title) at the bottom of a PNG."""
    try:
        im = Image.open(img_path).convert("RGB")
    except Exception as e:
        print(f"[3D] legend skipped (open failed): {e}")
        return

    W, H = im.size
    extra = max(70, int(H * 0.12))
    out = Image.new("RGB", (W, H + extra), (255, 255, 255))
    out.paste(im, (0, 0))
    draw = ImageDraw.Draw(out)

    # font (fallback friendly)
    try:
        f_big = ImageFont.truetype("DejaVuSans-Bold.ttf", 20)
        f = ImageFont.truetype("DejaVuSans.ttf", 18)
    except Exception:
        f_big = f = ImageFont.load_default()

    y = H + (extra - 20) // 2
    pad = 18
    box = 16

    # left (red)
    x = pad
    draw.rectangle((x, y, x + box, y + box), fill=left_rgb)
    draw.text((x + box + 8, y - 2), f"Red: {left_label}", fill=(0, 0, 0), font=f)

    # right (blue)
    x2 = x + 360
    draw.rectangle((x2, y, x2 + box, y + box), fill=right_rgb)
    draw.text((x2 + box + 8, y - 2), f"Blue: {right_label}", fill=(0, 0, 0), font=f)

    # optional title (above legend, aligned left)
    if title:
        draw.text((pad, H - 26), title, fill=(0, 0, 0), font=f_big)

    try:
        out.save(img_path)
    except Exception as e:
        print(f"[3D] legend save failed: {e}")


def _pymol_cli_render(script_text: str, exe: str | None = None):
    exe = exe or os.environ.get("PYMOL_BIN") or _PYMOL_BIN_CACHED
    if not exe or not os.path.isfile(exe):
        raise RuntimeError("PyMOL CLI not found. Set PYMOL_BIN or install pymol in PATH.")
    with tempfile.NamedTemporaryFile("w", suffix=".pml", delete=False) as tf:
        tf.write(script_text)
        tf.flush()
        cmdline = f"{shlex.quote(exe)} -cq {shlex.quote(tf.name)}"
        subprocess.run(cmdline, shell=True, check=True)


def _pred_pdb_path(pair_id: str, family: str, ver: str, cluster: str, pdbid: str, chain: str) -> str:
    """family = 'AF' or 'ESM'; ver e.g. '2' or '3'."""
    tag = _normalize_cluster_tag(cluster)
    base = f"output_{family}/{family}{ver}/{tag}__{pdbid}{chain}.pdb"
    return os.path.join("Pipeline", "FoldPairs", pair_id, base)


def _best_cluster_from_table(pair_id: str, table: str, model_tag: str, fold_idx: int) -> tuple[str, float] | None:
    """
    Generic best-cluster picker from df_{table}.csv.
    Returns (canonical_cluster_tag, TMscore) for fold_idx {0,1}.
    """
    df_path = os.path.join('Pipeline', 'FoldPairs', pair_id, 'Analysis', f'df_{table}.csv')
    if not os.path.isfile(df_path):
        return None
    try:
        d = pd.read_csv(df_path)
    except Exception:
        return None

    if 'model' in d.columns:
        d = d[d['model'].astype(str).str.upper() == model_tag.upper()]
        if d.empty:
            return None

    score_col = 'TMscore_fold1' if fold_idx == 0 else 'TMscore_fold2'
    ccol = 'cluster_num' if 'cluster_num' in d.columns else ('cluster' if 'cluster' in d.columns else None)
    if (score_col not in d.columns) or (ccol is None):
        return None

    d = d.copy()
    d['_score'] = pd.to_numeric(d[score_col], errors='coerce')
    d = d.dropna(subset=['_score'])
    if d.empty:
        return None

    row = d.loc[d['_score'].idxmax()]
    tag = _normalize_cluster_tag(str(row[ccol]).strip())
    return tag, float(row['_score'])


# ---------------------------------------------------------------------------
# Rendering functions
# ---------------------------------------------------------------------------
def _render_true_structures(pair_dir, pdb1, pdb2, out_dir, out_prefix, make_interactive: bool = False):
    """Save (1) two_structures.png and (2) two_structures_aligned.png, with legends."""
    out_two = os.path.join(out_dir, f"{out_prefix}_two_structures.png")
    out_aln = os.path.join(out_dir, f"{out_prefix}_two_structures_aligned.png")
    os.makedirs(out_dir, exist_ok=True)

    label_red  = os.path.splitext(pdb1)[0]
    label_blue = os.path.splitext(pdb2)[0]

    if not PYMOL_AVAILABLE:
        exe = os.environ.get("PYMOL_BIN") or _PYMOL_BIN_CACHED
        if exe:
            p1 = os.path.join(pair_dir, pdb1)
            p2 = os.path.join(pair_dir, pdb2)
            script = f"""
reinitialize
load "{p1}", fold1
load "{p2}", fold2
color red, fold1
color blue, fold2
bg_color white
set orthoscopic, on
viewport 1200, 900
zoom all, 10
png {out_two}, dpi=220, width=1200, height=900
align fold2, fold1
png {out_aln}, dpi=220, width=1200, height=900
quit
""".lstrip()
            _pymol_cli_render(script)
            _add_png_legend(out_two,  label_red, label_blue, title=None)
            _add_png_legend(out_aln,  label_red, label_blue, title=None)
            return
        print("[PyMOL] skipped (neither import nor CLI available)")
        return

    pymol.finish_launching(['pymol', '-cq'])
    try:
        cmd.delete('all')
        cmd.load(os.path.join(pair_dir, pdb1), 'fold1')
        cmd.load(os.path.join(pair_dir, pdb2), 'fold2')
        cmd.color('red',  'fold1')
        cmd.color('blue', 'fold2')
        cmd.bg_color('white')
        cmd.set('orthoscopic', 1)
        cmd.zoom('all', buffer=10)
        cmd.png(out_two, dpi=220, width=1200, height=900)
        cmd.align('fold2', 'fold1')
        cmd.png(out_aln, dpi=220, width=1200, height=900)
    finally:
        cmd.quit()
    _add_png_legend(out_two,  label_red, label_blue, title=None)
    _add_png_legend(out_aln,  label_red, label_blue, title=None)

    # Optional interactive (3D)
    if make_interactive:
        try:
            html_unaln = os.path.join(out_dir, f"{out_prefix}_two_structures.html")
            write_py3dmol_overlay_html(
                os.path.join(pair_dir, pdb1),
                os.path.join(pair_dir, pdb2),
                html_unaln,
                left_label=label_red,
                right_label=label_blue,
                title="True vs True (unaligned)"
            )
            print(f"[3D] wrote {html_unaln}")
        except Exception as e:
            print(f"[3D] true-true interactive (unaligned) failed: {e}")

        try:
            html_aln = os.path.join(out_dir, f"{out_prefix}_two_structures_aligned.html")
            write_py3dmol_overlay_html(
                os.path.join(pair_dir, pdb1),
                os.path.join(pair_dir, pdb2),
                html_aln,
                left_label=label_red,
                right_label=label_blue,
                title="True vs True (aligned)"
            )
            print(f"[3D] wrote {html_aln}")
        except Exception as e:
            print(f"[3D] true-true interactive (aligned) failed: {e}")


def _best_cluster_from_df(pair_id: str, model_ver: str, fold_idx: int) -> str | None:
    """
    Return canonical cluster tag ('DeepMsa' or 'ShallowMsa_###') for the best TMscore
    in Pipeline/<pair>/Analysis/df_af.csv for AF{model_ver} and fold_idx in {0,1}.
    """
    df_path = os.path.join('Pipeline', 'FoldPairs', pair_id, 'Analysis', 'df_af.csv')
    if not os.path.isfile(df_path):
        return None
    d = pd.read_csv(df_path)

    if 'model' in d.columns:
        want = f"AF{model_ver}"
        d = d[d['model'].astype(str).str.upper() == want.upper()]
        if d.empty:
            return None

    score_col = 'TMscore_fold1' if fold_idx == 0 else 'TMscore_fold2'
    if score_col not in d.columns:
        return None

    ccol = 'cluster_num' if 'cluster_num' in d.columns else ('cluster' if 'cluster' in d.columns else None)
    if ccol is None:
        return None

    row = d.loc[d[score_col].astype(float).fillna(-np.inf).idxmax()]
    tag = str(row[ccol]).strip()
    if tag.lower().startswith('deep'):
        return 'DeepMsa'
    m = re.fullmatch(r'\d+', tag)
    if m:
        return f"ShallowMsa_{int(m.group(0)):03d}"
    return tag


def _best_esm_from_df(pair_id: str, model_ver: str, fold_idx: int):
    """
    Return (cluster_tag, sample_idx) for the best TMscore in df_esm.csv
    for ESM{model_ver} and fold_idx in {0,1}.
    """
    df_path = os.path.join('Pipeline', 'FoldPairs', pair_id, 'Analysis', 'df_esm.csv')
    if not os.path.isfile(df_path):
        return None
    d = pd.read_csv(df_path)

    want = f"ESM{model_ver}".upper()
    if 'model' in d.columns:
        d = d[d['model'].astype(str).str.upper() == want]
        if d.empty:
            return None

    score_col = 'TMscore_fold1' if fold_idx == 0 else 'TMscore_fold2'
    if score_col not in d.columns or d[score_col].isna().all():
        return None

    row = d.loc[d[score_col].astype(float).fillna(-np.inf).idxmax()]

    if 'cluster_num' in d.columns and pd.notna(row['cluster_num']):
        tag_raw = str(row['cluster_num']).strip()
    elif 'cluster' in d.columns and pd.notna(row['cluster']):
        tag_raw = str(row['cluster']).strip()
    else:
        return None
    cluster_tag = _normalize_cluster_tag(tag_raw)

    sample_col = next((c for c in ('sample', 'sample_id', 'sample_num') if c in d.columns), None)
    sample_idx = int(row[sample_col]) if sample_col and pd.notna(row[sample_col]) else 1

    return cluster_tag, sample_idx


def _render_true_vs_best_esm(pair_id: str, pdbids: list[str], pdbchains: list[str], fig_dir_root: str,
                             model_ver='2', make_interactive: bool = False):
    """For each fold, overlay true structure vs best ESM{ver} prediction."""
    for idx in (0, 1):
        picked = _best_esm_from_df(pair_id, model_ver, idx)
        if not picked:
            print(f"[3D] no best row in df_esm.csv for ESM{model_ver}, fold{idx+1}")
            continue
        cluster_tag, sample_idx = picked

        true_pdb = os.path.join('Pipeline', 'FoldPairs', pair_id, f"{pdbids[idx]}.pdb")
        pred_pdb = os.path.join(
            'Pipeline', 'FoldPairs', pair_id, 'output_esm_fold', f"esm{model_ver}",
            f"{cluster_tag}__sample_{sample_idx:03d}_esm{model_ver}.pdb"
        )
        if not (os.path.isfile(true_pdb) and os.path.isfile(pred_pdb)):
            print(f"[3D] missing PDB for ESM{model_ver}: {true_pdb} / {pred_pdb}")
            continue

        out = os.path.join(fig_dir_root, f"{pair_id}_fold{idx+1}_vs_best_ESM{model_ver}.png")
        try:
            align_and_visualize_proteins(true_pdb, pred_pdb, out, open_environment=True)
            print(f"[3D] wrote {out}")
        except Exception as e:
            print(f"[3D] ESM{model_ver} fold{idx+1} overlay failed: {e}")

        if make_interactive:
            title = f"ESM{model_ver}({_cluster_short_disp(cluster_tag)}) vs. {pdbids[idx]}{pdbchains[idx]}"
            html_out = os.path.join(fig_dir_root, f"{pair_id}_fold{idx+1}_vs_best_ESM{model_ver}.html")
            try:
                write_py3dmol_overlay_html(
                    true_pdb, pred_pdb, html_out,
                    left_label=f"{pdbids[idx]}{pdbchains[idx]}",
                    right_label=f"ESM{model_ver}",
                    title=title, color_left="red", color_right="blue"
                )
                print(f"[3D] wrote {html_out}")
            except Exception as e:
                print(f"[3D] ESM{model_ver} interactive failed: {e}")


def _render_true_vs_best_models_generic(pair_id: str, pdbids: list[str], pdbchains: list[str],
                                        fig_dir_root: str, *, family: str, ver: str, make_interactive: bool = False):
    """
    family='AF' or 'ESM'; ver='2' or '3'
    For each fold, overlay true structure vs best {family}{ver} prediction.
    """
    exe = os.environ.get("PYMOL_BIN") or _PYMOL_BIN_CACHED
    if not PYMOL_AVAILABLE and not exe:
        print(f"[3D] Skipping {family}{ver}: no PyMOL")
        return

    os.makedirs(fig_dir_root, exist_ok=True)
    table = 'af' if family == 'AF' else 'esm'
    model_tag = f"{family}{ver}"

    for idx in (0, 1):
        pick = _best_cluster_from_table(pair_id, table=table, model_tag=model_tag, fold_idx=idx)
        if not pick:
            print(f"[3D] no best cluster in df_{table}.csv for {model_tag}, fold{idx+1}")
            continue

        cluster, tm = pick
        true_pdb = os.path.join('Pipeline', 'FoldPairs', pair_id, f"{pdbids[idx]}.pdb")
        pred_pdb = _pred_pdb_path(pair_id, family, ver, cluster, pdbids[idx], pdbchains[idx])

        if not (os.path.isfile(true_pdb) and os.path.isfile(pred_pdb)):
            print(f"[3D] missing PDB for {model_tag}: {true_pdb} / {pred_pdb}")
            continue

        out = os.path.join(fig_dir_root, f"{pair_id}_fold{idx+1}_vs_best_{model_tag}.png")
        try:
            align_and_visualize_proteins(true_pdb, pred_pdb, out, open_environment=bool(PYMOL_AVAILABLE))
            left = f"{pdbids[idx]}{pdbchains[idx]}"
            right = f"{model_tag}"
            title = f"{model_tag}({_cluster_short_disp(cluster)}) vs. {left}: TM={tm:.2f}"
            _add_png_legend(out, left, right, title=title)
            print(f"[3D] wrote {out}")
        except Exception as e:
            print(f"[3D] {model_tag} fold{idx+1} overlay failed: {e}")

        if make_interactive:
            html_out = os.path.join(fig_dir_root, f"{pair_id}_fold{idx+1}_vs_best_{model_tag}.html")
            try:
                write_py3dmol_overlay_html(
                    true_pdb, pred_pdb, html_out,
                    left_label=f"{pdbids[idx]}{pdbchains[idx]}",
                    right_label=model_tag,
                    title=title, color_left="red", color_right="blue"
                )
                print(f"[3D] wrote {html_out}")
            except Exception as e:
                print(f"[3D] {model_tag} interactive failed: {e}")


def _render_ddg_aligned_image(pair_id: str, pdbids: list[str], pdbchains: list[str], fig_dir_root: str) -> None:
    """Build the aligned per-residue DDG image into output_figs/."""
    os.makedirs(fig_dir_root, exist_ok=True)
    pair_dir = os.path.join("Pipeline", "FoldPairs", pair_id)
    anal_dir = os.path.join(pair_dir, "Analysis")
    out_png = os.path.join(fig_dir_root, f"{pair_id}_ddg_aligned.png")

    print("fig_dir_root:", fig_dir_root)
    print("out_png:", out_png)

    # Load per-structure energy CSVs
    res_lists = []
    for i in range(2):
        tag4 = pdbids[i][:4]
        csv_path = os.path.join(anal_dir, f"df_energy_{tag4}.csv")
        txt_path = os.path.join(pair_dir, "output_deltaG", f"deltaG_{tag4}.txt")
        if os.path.isfile(csv_path):
            import pandas as _pd
            df_e = _pd.read_csv(csv_path)
            res_lists.append(df_e.to_dict("records"))
        elif os.path.isfile(txt_path):
            res_lists.append(read_energy_tuples(txt_path))
        else:
            print(f"[DDG] no energy data for {tag4} in {pair_id}; skipping DDG image")
            return

    try:
        align_and_compare_residues(
            res_lists[0], res_lists[1],
            f"{pdbids[0]}{pdbchains[0]}", f"{pdbids[1]}{pdbchains[1]}",
            output_file=out_png
        )
        print(f"[DDG] wrote {out_png}")
    except Exception as e:
        print(f"[DDG] align/plot failed for {pair_id}: {e}")


def write_py3dmol_overlay_html(
    pdb1: str,
    pdb2: str,
    out_html: str,
    *,
    left_label: str = "Fold1",
    right_label: str = "Fold2",
    title: str | None = None,
    color_left: str = "red",
    color_right: str = "blue",
    width: int = 700,
    height: int = 520,
):
    import py3Dmol, html as _html

    v = py3Dmol.view(width=width, height=height)
    with open(pdb1) as f:
        v.addModel(f.read(), "pdb")
    v.setStyle({"model": 0}, {"cartoon": {"color": color_left}})
    with open(pdb2) as f:
        v.addModel(f.read(), "pdb")
    v.setStyle({"model": 1}, {"cartoon": {"color": color_right}})
    v.setBackgroundColor("white")
    v.zoomTo()

    body = v._make_html()
    title_html = f"<h3 style='font-family:sans-serif;margin:6px 0 6px 0;'>{_html.escape(title) if title else ''}</h3>"
    legend_html = f"""
<div style="display:flex;justify-content:center;gap:18px;margin:4px 0 10px 0;font-family:sans-serif;">
  <span><span style="display:inline-block;width:12px;height:12px;background:{color_left};margin-right:6px;border:1px solid #444;"></span>{_html.escape(left_label)}</span>
  <span><span style="display:inline-block;width:12px;height:12px;background:{color_right};margin-right:6px;border:1px solid #444;"></span>{_html.escape(right_label)}</span>
</div>"""

    with open(out_html, "w") as f:
        f.write(title_html + legend_html + body)


def align_and_visualize_proteins(pdb_file1, pdb_file2, output_file, open_environment=True):
    """Align two protein structures and save the visualization as an image."""
    out_dir = os.path.dirname(os.path.abspath(output_file)) or "."
    os.makedirs(out_dir, exist_ok=True)

    if PYMOL_AVAILABLE:
        try:
            if open_environment:
                pymol.finish_launching(['pymol', '-cq'])
            cmd.delete('all')
            cmd.load(pdb_file1, 'protein1')
            cmd.load(pdb_file2, 'protein2')
            cmd.align('protein2', 'protein1')
            cmd.color('red', 'protein1')
            cmd.color('blue', 'protein2')
            cmd.bg_color('white')
            cmd.zoom('all', buffer=10)
            cmd.png(output_file, dpi=220, width=1200, height=900)
        finally:
            if open_environment:
                try: cmd.quit()
                except Exception: pass
        return

    # ---- CLI fallback ----
    exe = os.environ.get("PYMOL_BIN") or _PYMOL_BIN_CACHED
    if not exe:
        raise RuntimeError("PyMOL module not available and no CLI found.")

    p1 = os.path.abspath(pdb_file1)
    p2 = os.path.abspath(pdb_file2)
    out = os.path.abspath(output_file)

    script = f"""
reinitialize
load "{p1}", protein1
load "{p2}", protein2
align protein2, protein1
color red, protein1
color blue, protein2
bg_color white
viewport 1200, 900
zoom all, 10
set antialias, 2
png {out}, dpi=220, width=1200, height=900
quit
""".lstrip()

    if '_pymol_cli_render' in globals() and callable(globals()['_pymol_cli_render']):
        _pymol_cli_render(script)
    else:
        with tempfile.NamedTemporaryFile("w", suffix=".pml", delete=False) as tf:
            tf.write(script); tf.flush()
            cmdline = f"{shlex.quote(PYMOL_BIN)} -cq {shlex.quote(tf.name)}"
            subprocess.run(cmdline, shell=True, check=True)
