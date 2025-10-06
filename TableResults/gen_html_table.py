# File: TableResults/gen_html_table.py
import os, sys, re, html, shutil
import pandas as pd
import argparse

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, ROOT)
from config import *
from Analysis.postprocess_unified import build_unified_tables_from_cluster_dfs

# --- add near the top (after imports / _ensure_unified_csvs) ---

# Single source of truth for the main comparison table
DEFAULT_PREFERRED_COLS = [
    "pair_id", "#RES", "MSA: DEPTH; #RES; #Clusters", "PAIR_TM",
    "ΔG1", "ΔG2",
    "AF2Clust_TM1","AF2Clust_TM2","AF2Deep_TM1","AF2Deep_TM2",
    "AF3Clust_TM1","AF3Clust_TM2","AF3Deep_TM1","AF3Deep_TM2",
    "ESM2_TM1","ESM2_TM2","ESM3_TM1","ESM3_TM2",
    "MSATrans_CMAP_PR1","MSATrans_CMAP_PR2","MSATrans_CMAP_RE1","MSATrans_CMAP_RE2",
]

DEFAULT_EXPLANATIONS = {
    "#RES": "Number of residues in the longer chain of the pair.",
    "MSA: DEPTH; #RES; #Clusters": "DEPTH = number of sequences in DeepMsa.a3m; #RES = alignment width (columns, including gaps); #Clusters = count of ShallowMsa_* a3m files.",
    "PAIR_TM": "TM-score between the two ground-truth folds (max of the two directions).",
    "ΔG1": "Rosetta full-atom (Ref2015) energy of Fold1 (lower is more stable).",
    "ΔG2": "Rosetta full-atom (Ref2015) energy of Fold2 (lower is more stable).",
    "AF2Clust_TM1": "Best TM-score to Fold1 among AF2 predictions built from any shallow cluster (number in parentheses is the cluster id).",
    "AF2Clust_TM2": "Best TM-score to Fold2 among AF2 predictions from shallow clusters.",
    "AF2Deep_TM1":  "Best TM-score to Fold1 among AF2 predictions built from the DeepMsa alignment.",
    "AF2Deep_TM2":  "Best TM-score to Fold2 among AF2 predictions from DeepMsa.",
    "AF3Clust_TM1": "Best TM-score to Fold1 among AF3 predictions from shallow clusters.",
    "AF3Clust_TM2": "Best TM-score to Fold2 among AF3 predictions from shallow clusters.",
    "AF3Deep_TM1":  "Best TM-score to Fold1 among AF3 predictions from DeepMsa.",
    "AF3Deep_TM2":  "Best TM-score to Fold2 among AF3 predictions from DeepMsa.",
    "ESM2_TM1": "Best TM-score to Fold1 among all ESMFold(ESM2) predictions (across all sampled sequences from all clusters; no MSA used).",
    "ESM2_TM2": "Best TM-score to Fold2 among all ESMFold(ESM2) predictions.",
    "ESM3_TM1": "Best TM-score to Fold1 among all ESMFold(ESM3) predictions.",
    "ESM3_TM2": "Best TM-score to Fold2 among all ESMFold(ESM3) predictions.",
    "MSATrans_CMAP_PR1": "Maximum precision of MSA-Transformer contact map vs Fold1 truth (threshold=0.4; |i−j|≥6).",
    "MSATrans_CMAP_PR2": "Maximum precision vs Fold2 truth (same settings).",
    "MSATrans_CMAP_RE1": "Maximum recall vs Fold1 truth.",
    "MSATrans_CMAP_RE2": "Maximum recall vs Fold2 truth.",
}


# === NEW: ensure unified CSVs exist (or rebuild on demand) ===
def _ensure_unified_csvs(force_rerun: bool = False):
    need = force_rerun or (not os.path.exists(SUMMARY_RESULTS_TABLE)) or (not os.path.exists(DETAILED_RESULTS_TABLE))
    if need:
        # local import to avoid cycles when this module is imported elsewhere
        try:
            build_unified_tables_from_cluster_dfs(write_out=True)
            print("[gen_html] ensured unified CSVs (summary + detailed).")
        except Exception as e:
            print(f"[gen_html] WARN: could not build unified CSVs automatically: {e}")

NUM_RE = re.compile(r"^\s*([+-]?\d+(?:\.\d+)?)")

def _append_mean_std_row(df: pd.DataFrame, label_col: str, label_text: str = "Averages") -> pd.DataFrame:
    """Append a last row with per-column mean (std) for numeric columns.
       Non-numeric columns get '-' except the label_col which gets label_text."""
    out = df.copy()
    means = {}
    for c in out.columns:
        if c == label_col:
            means[c] = label_text
            continue
        # try parsing numeric (also when cells are strings like "0.78 (7)")
        s = out[c]
        if s.dtype.kind in "biufc":
            x = pd.to_numeric(s, errors="coerce")
        else:
            x = pd.to_numeric(s.apply(numeric_part), errors="coerce")
        if x.notna().any():
            mu = float(x.mean())
            sd = float(x.std(ddof=1)) if x.count() > 1 else 0.0
            means[c] = f"{mu:.3f} ({sd:.3f})"
        else:
            means[c] = "-"
    out = pd.concat([out, pd.DataFrame([means])], ignore_index=True)
    return out


def numeric_part(cell) -> str:
    if pd.isna(cell):
        return ""
    s = str(cell).strip()
    if s == "-":
        return ""
    m = NUM_RE.match(s)
    return m.group(1) if m else s


def write_global_tables(*, force_rerun_csv: bool = False,
                        fade_min_clusters: int = 2) -> tuple[str, str]:
    """
    Build both global HTML tables (summary + cluster-level) using the single
    source-of-truth column order and explanations. Returns (summary_html, clusters_html).
    """
    _ensure_unified_csvs(force_rerun=force_rerun_csv)

    out_summary = os.path.join(TABLES_RES, "protein_comparison_table.html")
    gen_html_from_summary_table(
        summary_csv=SUMMARY_RESULTS_TABLE,
        output_html=out_summary,
        preferred_column_order=DEFAULT_PREFERRED_COLS,
        column_explanations=DEFAULT_EXPLANATIONS,
        fade_min_clusters=fade_min_clusters,  # << NEW
    )

    out_clusters = os.path.join(TABLES_RES, "protein_clusters_table.html")
    gen_html_from_cluster_detailed_table(
        clusters_csv=DETAILED_RESULTS_TABLE,
        output_html=out_clusters,
    )
    return out_summary, out_clusters


def gen_html_from_summary_table(
    summary_csv: str | None = None,
    output_html: str | None = None,
    title: str = "Interactive Protein Comparison Table",
    base_pair_url: str | None = GITHUB_URL_HTML + "/{pair_id}.html",
    preferred_column_order: list[str] | None = None,
    column_explanations: dict[str, str] | None = None,
    fade_min_clusters: int | None = None,   # << NEW
) -> str:
    if summary_csv is None:
        summary_csv = SUMMARY_RESULTS_TABLE
    if output_html is None:
        output_html = os.path.join(TABLES_RES, "table.html")

    # Ensure CSVs exist if missing
    _ensure_unified_csvs(force_rerun=False)

    if not os.path.exists(summary_csv):
        raise FileNotFoundError(f"Summary CSV not found: {summary_csv}")

    df = pd.read_csv(summary_csv)

    if df.empty:
        os.makedirs(os.path.dirname(output_html), exist_ok=True)
        with open(output_html, "w", encoding="utf-8") as f:
            f.write(f"<!DOCTYPE html><html><head><meta charset='utf-8'><title>{html.escape(title)}</title></head>"
                    f"<body><h2>{html.escape(title)}</h2><p>No data available.</p></body></html>")
        return output_html

    # Identify the column that carries cluster count in parentheses
    cluster_col = None
    for c in df.columns:
        if "(#Clusters)" in c:
            cluster_col = c
            break
    cluster_col = "MSA: DEPTH; #RES; #Clusters" if "MSA: DEPTH; #RES; #Clusters" in df.columns else None

    # Helper to parse "(N)" from the “MSA DEPTH (#Clusters)” cell
    _paren = re.compile(r"\((\d+)\)")

    def _clusters_in_row(row) -> int | None:
        if cluster_col is None:
            return None
        val = str(row.get(cluster_col, ""))
        try:
            # split "DEPTH; WIDTH; CLUSTERS"
            parts = [p.strip() for p in val.split(";")]
            return int(parts[-1]) if len(parts) == 3 else None
        except Exception:
            return None


    # --- Normalize pair id column name to 'pair_id' (for display & linking) ---
    if "pair_id" in df.columns:
        pair_col = "pair_id"
    elif "fold_pair" in df.columns:
        df = df.rename(columns={"fold_pair": "pair_id"})
        pair_col = "pair_id"
    else:
        # Fallback: first column becomes pair id
        pair_col = df.columns[0]
        df = df.rename(columns={pair_col: "pair_id"})
        pair_col = "pair_id"

    # --- Build display column order: pair_id first, then preferred order (excluding pair_id), then the rest ---
    cols_all = list(df.columns)
    # Remove any duplicate/alias of pair id from the pool
    remaining = [c for c in cols_all if c != "pair_id"]

    ordered = ["pair_id"]
    if preferred_column_order:
        # Only keep preferred columns that actually exist and are not pair_id
        wanted = [c for c in preferred_column_order if c in remaining]
        ordered += wanted
        # Add everything else not already included
        the_rest = [c for c in remaining if c not in wanted]
        ordered += the_rest
    else:
        ordered += remaining

    df = df[ordered]

    # Round numeric columns to 3 decimals for display
    for c in df.columns:
        if c == "pair_id": continue
        s = pd.to_numeric(df[c].apply(numeric_part), errors="coerce")
        if s.notna().any():
            # replace display values with rounded (preserve cluster-id parentheses in AF2Clust_TM*)
            def _fmt(v):
                if pd.isna(v): return v
                v = str(v)
                # keep things like "0.812 (7)": round the first number and keep the parenthesis
                m = re.match(r"^\s*([+-]?\d+(?:\.\d+)?)\s*(.*)$", v)
                if not m: return v
                a, tail = m.group(1), m.group(2)
                return f"{float(a):.3f}{tail}"
            df[c] = df[c].map(_fmt)

    # Append "Averages (std)" row
    df = _append_mean_std_row(df, label_col="pair_id", label_text="Averages")



    # --- Table header (clickable for sorting) ---
    thead = "<tr>" + "".join(
        f'<th onclick="sortTable({i})">{html.escape(col)}</th>'
        for i, col in enumerate(df.columns)
    ) + "</tr>"

    # For sorting: pull numeric prefix from strings like "3219 (18)" -> "3219"

    # --- Rows. When building rows: choose row class based on cluster count ---
    rows = []
    for _, r in df.iterrows():
        pair = str(r["pair_id"])
        link = (base_pair_url or "{pair_id}.html").format(pair_id=html.escape(pair))
        tds = [f'<td><a href="{link}" target="_blank">{html.escape(pair)}</a></td>']
        for col in df.columns[1:]:
            val = r[col]
            if col == "#RES":
                # show integer in normal rows
                try:
                    ival = int(float(val))
                    disp = str(ival)
                    sortv = str(ival)
                except Exception:
                    disp = "-" if pd.isna(val) else str(val)
                    sortv = numeric_part(val)
                tds.append(f'<td data-sort-value="{html.escape(str(sortv))}">{html.escape(disp)}</td>')
            else:
                disp = "-" if pd.isna(val) else str(val)
                tds.append(f'<td data-sort-value="{html.escape(numeric_part(val))}">{html.escape(disp)}</td>')
        row_cls = ""
        if fade_min_clusters is not None:
            ncl = _clusters_in_row(r)
            if ncl is not None and ncl < fade_min_clusters:
                row_cls = ' class="lowclust"'
        rows.append(f"<tr{row_cls}>" + "".join(tds) + "</tr>")

    # ---- Averages row (2 decimals). Special handling for "MSA: DEPTH; #RES; #Clusters" ----
    avg_label = '<a href="https://orzuk.github.io/MsaCluster/pairs_global_analysis.html" target="_blank">Average</a>'
    avg_tds = [f'<td>{avg_label}</td>']

    # Pre-compute means for the triple column if present
    triple_col = "MSA: DEPTH; #RES; #Clusters"
    triple_disp = ""
    if triple_col in df.columns:
        depths, reses, clusts = [], [], []
        for s in df[triple_col].dropna():
            parts = [p.strip() for p in str(s).split(";")]
            def num(x):
                try: return float(numeric_part(x))
                except: return float("nan")
            vals = [num(p) for p in parts] + [float("nan")] * 3
            d, rnum, c = vals[0], vals[1], vals[2]
            if not pd.isna(d):    depths.append(d)
            if not pd.isna(rnum): reses.append(rnum)
            if not pd.isna(c):    clusts.append(c)
        def mean_or_nan(a):
            return (sum(a)/len(a)) if a else float("nan")
        md = mean_or_nan(depths)
        mr = mean_or_nan(reses)
        mc = mean_or_nan(clusts)
        triple_disp = f"{md:.2f} ; {mr:.2f} ; {mc:.2f}"

    # Build cells
    for col in df.columns[1:]:
        if col == triple_col:
            avg_tds.append(f'<td data-sort-value="">{html.escape(triple_disp)}</td>')
            continue
        s = pd.to_numeric(df[col].map(numeric_part), errors="coerce")
        if s.notna().any():
            mu = s.mean()
            sd = s.std(ddof=1) if s.count() > 1 else 0.0
            disp = f"{mu:.2f} ({sd:.2f})"
            avg_tds.append(f'<td data-sort-value="{mu:.6f}">{html.escape(disp)}</td>')
        else:
            avg_tds.append('<td data-sort-value=""></td>')
    rows.append("<tr>" + "".join(avg_tds) + "</tr>")


    # --- Explanations (only for existing columns) ---
    expl_lines = []
    if column_explanations:
        for c in df.columns:
            if c in column_explanations:
                expl_lines.append(f"<p><b>{html.escape(c)}</b>: {html.escape(column_explanations[c])}</p>")


    # --- update the legend block build (append a note if fading is active) ---
    if fade_min_clusters is not None:
        expl_lines.append(
            f"<p><i>Note:</i> rows shown dimmed have fewer than {fade_min_clusters} clusters in the shallow MSA and are de-emphasized.</p>"
        )
    expl_html = ("\n".join(expl_lines)) if expl_lines else ""

    # --- HTML doc ---
    html_doc = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8" />
<title>{html.escape(title)}</title>
<style>
  body {{
    font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
    background-color: #121212; color: #E0E0E0; margin: 20px;
  }}
  table {{
    width: 98%; margin: auto; border-collapse: collapse;
    box-shadow: 0 4px 8px rgba(0,0,0,0.5);
  }}
  th, td {{ border: 1px solid #333; padding: 10px 14px; text-align: left; white-space: nowrap; }}
  th {{ background-color: #b71c1c; color: #fff; font-size: 16px; cursor: pointer; position: sticky; top: 0; }}
  tr:nth-child(even) {{ background-color: #2b2b2b; }}
  tr:hover {{ background-color: #3a3a3a; }}
  a {{ color: #64B5F6; text-decoration: none; }}
  a:hover {{ text-decoration: underline; }}
  h2 {{ text-align: center; color: #E0E0E0; margin-top: 0; }}
  .legend {{ width: 98%; margin: 20px auto; line-height: 1.4; }}
    tr.lowclust{{opacity: .45; background - color:  # 1b1b1b; }}
    tr.lowclust:hover{{opacity: .65; background - color:  # 272727; }}
</style>
</head>
<body>
<h2>{html.escape(title)}</h2>
<table id="tbl">
  <thead>{thead}</thead>
  <tbody>
    {"".join(rows)}
  </tbody>
</table>
<div class="legend">
{expl_html}
</div>
<script>
let sortAsc = [];
function getCellSortValue(td) {{
  return td.getAttribute('data-sort-value') ?? td.textContent.trim();
}}
function sortTable(colIdx) {{
  const table = document.getElementById('tbl');
  const tbody = table.tBodies[0];
  const rows = Array.from(tbody.querySelectorAll('tr'));
  sortAsc[colIdx] = !sortAsc[colIdx];
  rows.sort((r1, r2) => {{
    const t1 = getCellSortValue(r1.cells[colIdx]);
    const t2 = getCellSortValue(r2.cells[colIdx]);
    const n1 = parseFloat(t1), n2 = parseFloat(t2);
    const bothNumeric = !isNaN(n1) && !isNaN(n2);
    if (bothNumeric) return sortAsc[colIdx] ? (n1 - n2) : (n2 - n1);
    const s1 = t1.toLowerCase(), s2 = t2.toLowerCase();
    if (s1 < s2) return sortAsc[colIdx] ? -1 : 1;
    if (s1 > s2) return sortAsc[colIdx] ? 1 : -1;
    return 0;
  }});
  rows.forEach(r => tbody.appendChild(r));
}}
</script>
</body>
</html>"""

    os.makedirs(os.path.dirname(output_html), exist_ok=True)
    with open(output_html, "w", encoding="utf-8") as f:
        f.write(html_doc)
    print(f"[html] wrote: {output_html}")
    return output_html


def gen_html_from_cluster_detailed_table(
    clusters_csv: str | None = None,
    output_html: str | None = None,
    title: str = "Cluster-level Results (one row per cluster)",
    base_pair_url: str | None = GITHUB_URL_HTML + "/{pair_id}.html",
) -> str:
    if clusters_csv is None:
        # next to summary table
        clusters_csv = os.path.join(os.path.dirname(SUMMARY_RESULTS_TABLE), "clusters_global_table.csv")
    if output_html is None:
        output_html = os.path.join(TABLES_RES, "protein_clusters_table.html")

    if not os.path.exists(clusters_csv):
        raise FileNotFoundError(f"Clusters CSV not found: {clusters_csv}")

    df = pd.read_csv(clusters_csv)
    if df.empty:
        os.makedirs(os.path.dirname(output_html), exist_ok=True)
        with open(output_html, "w", encoding="utf-8") as f:
            f.write(f"<!DOCTYPE html><html><head><meta charset='utf-8'><title>{html.escape(title)}</title></head>"
                    f"<body><h2>{html.escape(title)}</h2><p>No data available.</p></body></html>")
        return output_html

    # columns order
    preferred = ["pair_id","cluster","n","neff",
                 "AF2_TM1","AF2_TM2","AF3_TM1","AF3_TM2",
                 "CMAP_RE1","CMAP_RE2",
                 "ESM_SEQIDS","ESM_TM1_LIST","ESM_TM2_LIST"]
    cols = [c for c in preferred if c in df.columns] + [c for c in df.columns if c not in preferred]
    df = df[cols]

    # build header
    thead = "<tr>" + "".join(
        f'<th onclick="sortTable({i})">{html.escape(col)}</th>' for i, col in enumerate(df.columns)
    ) + "</tr>"

    rows = []
    for _, r in df.iterrows():
        pair = str(r[pair_col])
        link = (base_pair_url or "{pair_id}.html").format(pair_id=html.escape(pair))
        tds = [f'<td><a href="{link}" target="_blank">{html.escape(pair)}</a></td>']
        for col in df.columns[1:]:
            val = r[col]

            # special-case: clean cluster label if "unknown"
            if col == "cluster":
                disp = "" if (isinstance(val, str) and "unknown" in val.lower()) else ("" if pd.isna(val) else str(val))
                tds.append(f'<td data-sort-value="{html.escape(numeric_part(val))}">{html.escape(disp)}</td>')
                continue

            # special-case: n / neff as integers
            if col in ("n", "neff"):
                try:
                    ival = int(float(val))
                    tds.append(f'<td data-sort-value="{ival}">{ival}</td>')
                    continue
                except Exception:
                    pass

            disp = "-" if pd.isna(val) or val == "" else str(val)
            extra_cls = ' class="esmcol"' if col.startswith("ESM_") else ""
            tds.append(f'<td{extra_cls} data-sort-value="{html.escape(numeric_part(val))}">{html.escape(disp)}</td>')
        rows.append("<tr>" + "".join(tds) + "</tr>")


    # --- Append Average row (2 decimals). Leave cluster cell empty. Skip string ESM columns ---
    avg_tds = ['<td>Average</td>', '<td></td>']  # pair_id label, empty cluster id
    for col in df.columns[2:]:
        if col.startswith("ESM_"):
            avg_tds.append('<td data-sort-value=""></td>')
            continue
        s = pd.to_numeric(df[col], errors="coerce")
        if s.notna().any():
            mu = s.mean()
            sd = s.std(ddof=1) if s.count() > 1 else 0.0
            disp = f"{mu:.2f} ({sd:.2f})"
            avg_tds.append(f'<td data-sort-value="{mu:.6f}">{html.escape(disp)}</td>')
        else:
            avg_tds.append('<td data-sort-value=""></td>')
    rows.append("<tr>" + "".join(avg_tds) + "</tr>")

    html_doc = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8" />
<title>{html.escape(title)}</title>
<style>
  body {{ font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; background-color: #121212; color: #E0E0E0; margin: 20px; }}
  table {{ width: 98%; margin: auto; border-collapse: collapse; box-shadow: 0 4px 8px rgba(0,0,0,0.5); }}
  th, td {{ border: 1px solid #333; padding: 10px 14px; text-align: left; white-space: nowrap; }}
  th {{ background-color: #b71c1c; color: #fff; font-size: 16px; cursor: pointer; position: sticky; top: 0; }}
  tr:nth-child(even) {{ background-color: #2b2b2b; }}
  tr:hover {{ background-color: #3a3a3a; }}
  a {{ color: #64B5F6; text-decoration: none; }}
  a:hover {{ text-decoration: underline; }}
  h2 {{ text-align: center; color: #E0E0E0; margin-top: 0; }}
  td.esmcol {{ white-space: normal; word-break: break-word; max-width: 520px; }}
</style>
</head>
<body>
<h2>{html.escape(title)}</h2>
<table id="tbl">
  <thead>{thead}</thead>
  <tbody>
    {"".join(rows)}
  </tbody>
</table>
<script>
let sortAsc = [];
function getCellSortValue(td) {{ return td.getAttribute('data-sort-value') ?? td.textContent.trim(); }}
function sortTable(colIdx) {{
  const table = document.getElementById('tbl');
  const tbody = table.tBodies[0];
  const rows = Array.from(tbody.querySelectorAll('tr'));
  sortAsc[colIdx] = !sortAsc[colIdx];
  rows.sort((r1, r2) => {{
    const t1 = getCellSortValue(r1.cells[colIdx]);
    const t2 = getCellSortValue(r2.cells[colIdx]);
    const n1 = parseFloat(t1), n2 = parseFloat(t2);
    const bothNumeric = !isNaN(n1) && !isNaN(n2);
    if (bothNumeric) return sortAsc[colIdx] ? (n1 - n2) : (n2 - n1);
    const s1 = t1.toLowerCase(), s2 = t2.toLowerCase();
    if (s1 < s2) return sortAsc[colIdx] ? -1 : 1;
    if (s1 > s2) return sortAsc[colIdx] ? 1 : -1;
    return 0;
  }});
  rows.forEach(r => tbody.appendChild(r));
}}
</script>
</body>
</html>"""
    os.makedirs(os.path.dirname(output_html), exist_ok=True)
    with open(output_html, "w", encoding="utf-8") as f:
        f.write(html_doc)
    print(f"[html] wrote: {output_html}")
    return output_html


def gen_html_for_global_plots(
    images_dir: str,
    output_html: str = os.path.join("docs", "pairs_global_analysis.html"),
    title: str = "Global Comparison Plots"
) -> str:
    import os, html

    # Resolve absolute dir (for safe listing no matter CWD)
    images_dir_abs = os.path.abspath(images_dir)
    if not os.path.isdir(images_dir_abs):
        raise FileNotFoundError(f"Images dir not found: {images_dir_abs}")

    files = sorted([
        f for f in os.listdir(images_dir_abs)
        if f.lower().endswith((".png",".jpg",".jpeg",".webp",".svg"))
    ])

    # Ensure output dir exists and compute a RELATIVE path for <img src=...>
    out_dir = os.path.abspath(os.path.dirname(output_html))
    os.makedirs(out_dir, exist_ok=True)

    # Copy figures to output figs dir
    os.makedirs(FIGURE_PUBLISH_DIR, exist_ok=True)
    src_files = [f for f in os.listdir(images_dir_abs) if f.lower().endswith(('.png','.jpg','.jpeg','.webp','.svg'))]
    published = []
    for f in sorted(src_files):
        src = os.path.join(images_dir_abs, f)
        dst = os.path.join(FIGURE_PUBLISH_DIR, f)
        if PUBLISH_FIGURES:
            shutil.copy2(src, dst)
        if os.path.isfile(dst):  # only link files present in publish dir
            published.append(f)

    rel_prefix = os.path.relpath(os.path.abspath(FIGURE_PUBLISH_DIR), start=out_dir)


    if not files:
        rows = "<p>No figures were found in the images directory.</p>"
    else:
        rows = "\n".join(
            f'''<figure style="margin:18px 0;">
  <figcaption style="margin:4px 0 8px;font-family:Segoe UI;opacity:.85">{html.escape(f)}</figcaption>
  <img loading="lazy" src="{html.escape(os.path.join(rel_prefix, f))}" style="max-width:100%;border:1px solid #333;"/>
</figure>'''
            for f in files
        )

    html_doc = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<title>{html.escape(title)}</title>
<style>
  :root {{ color-scheme: dark; }}
  body {{
    background:#121212; color:#E0E0E0; font-family:'Segoe UI',Tahoma,Verdana,sans-serif; margin:20px;
  }}
  h2 {{ margin: 0 0 12px; }}
</style>
</head>
<body>
<h2>{html.escape(title)}</h2>
{rows}
</body>
</html>"""
    with open(output_html, "w", encoding="utf-8") as f:
        f.write(html_doc)
    print(f"[html] wrote: {output_html}")
    return output_html


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--force_rerun", action="store_true",
                    help="Recompute unified CSVs from per-pair cluster CSVs before building HTML.")
    ap.add_argument("--min_clusters_fade", type=int, default=2,
                    help="Fade rows whose (#Clusters) is below this number. Set <=0 to disable.")
    args = ap.parse_args()

    write_global_tables(force_rerun_csv=args.force_rerun,
                        fade_min_clusters=None if args.min_clusters_fade <= 0 else args.min_clusters_fade)

    print("[OK main html] done.")
