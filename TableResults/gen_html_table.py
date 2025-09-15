# File: TableResults/gen_html_table.py
import os, sys, re
import html
import pandas as pd


ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, ROOT)

from config import *

# >>> Use your config — no ad-hoc path strings
from config import SUMMARY_RESULTS_TABLE, TABLES_RES  # , OUTPUT_PATH_NOTEBOOKS, DATA_DIR, MAIN_DIR

def gen_html_from_summary_table(
    summary_csv: str | None = None,
    output_html: str | None = None,
    title: str = "Interactive Protein Comparison Table",
    # Keep the website link pattern you’re already using
    base_pair_url: str | None = GITHUB_URL_HTML + "/{pair_id}.html",
    preferred_column_order: list[str] | None = None,
) -> str:
    """
    Build a sortable HTML table from the *new* summary CSV (best-per-pair table).

    - If summary_csv/output_html are not provided, they default to config:
        summary_csv  -> config.SUMMARY_RESULTS_TABLE
        output_html  -> os.path.join(config.TABLES_RES, 'table.html')
    - First column links to per-pair page (base_pair_url).
    - Sorts numerically using the numeric part of "0.53 (6)".

    Returns:
        The path of the written HTML file.
    """
    # Resolve inputs from config if omitted
    if summary_csv is None:
        summary_csv = SUMMARY_RESULTS_TABLE
    if output_html is None:
        output_html = os.path.join(TABLES_RES, "table.html")

    if not os.path.exists(summary_csv):
        raise FileNotFoundError(f"Summary CSV not found: {summary_csv}")

    df = pd.read_csv(summary_csv)

    # Minimal empty-page fallback
    if df.empty:
        os.makedirs(os.path.dirname(output_html), exist_ok=True)
        with open(output_html, "w", encoding="utf-8") as f:
            f.write(f"<!DOCTYPE html><html><head><meta charset='utf-8'><title>{html.escape(title)}</title></head>"
                    f"<body><h2>{html.escape(title)}</h2><p>No data available.</p></body></html>")
        return output_html

    # Column ordering
    pair_col = "pair_id" if "pair_id" in df.columns else "fold_pair"
    if pair_col not in df.columns:
        raise KeyError("Expected a 'pair_id' (or 'fold_pair') column in the summary CSV.")

    cols = list(df.columns)
    if preferred_column_order is None:
        best_cols = [c for c in cols if c.startswith("BEST_")]
        other_cols = [c for c in cols if c not in best_cols + [pair_col]]
        ordered_cols = [pair_col] + sorted(best_cols) + other_cols
    else:
        wanted = [c for c in preferred_column_order if c in cols]
        rest = [c for c in cols if c not in wanted]
        ordered_cols = wanted + rest

    df = df[ordered_cols]

    # Build header
    thead = "<tr>" + "".join(
        f'<th onclick="sortTable({i})">{html.escape(col)}</th>'
        for i, col in enumerate(df.columns)
    ) + "</tr>"

    # Helper: numeric part from "0.53 (6)"
    num_re = re.compile(r"^\s*([+-]?\d+(?:\.\d+)?)")
    def numeric_part(cell) -> str:
        if pd.isna(cell):
            return ""
        s = str(cell).strip()
        if s == "-":
            return ""
        m = num_re.match(s)
        return m.group(1) if m else s

    # Build rows (first col = link to pair page)
    rows = []
    for _, r in df.iterrows():
        pair = str(r[pair_col])
        link = (base_pair_url or "{pair_id}.html").format(pair_id=html.escape(pair))
        tds = [f'<td><a href="{link}" target="_blank">{html.escape(pair)}</a></td>']
        for col in df.columns[1:]:
            val = r[col]
            disp = "-" if pd.isna(val) else str(val)
            tds.append(
                f'<td data-sort-value="{html.escape(numeric_part(val))}">{html.escape(disp)}</td>'
            )
        rows.append("<tr>" + "".join(tds) + "</tr>")

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
    width: 90%; margin: auto; border-collapse: collapse;
    box-shadow: 0 4px 8px rgba(0,0,0,0.5);
  }}
  th, td {{ border: 1px solid #333; padding: 10px 14px; text-align: left; white-space: nowrap; }}
  th {{ background-color: #004d40; color: #fff; font-size: 16px; cursor: pointer; position: sticky; top: 0; }}
  tr:nth-child(even) {{ background-color: #2b2b2b; }}
  tr:hover {{ background-color: #3a3a3a; }}
  a {{ color: #BB86FC; text-decoration: none; }}
  a:hover {{ text-decoration: underline; }}
  h2 {{ text-align: center; color: #E0E0E0; margin-top: 0; }}
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
    detailed_csv: str | None = None,
    output_html: str | None = None,
    title: str = "Cluster-level Results (one row per cluster)",
    base_pair_url: str | None = GITHUB_URL_HTML + "/{pair_id}.html",
) -> str:
    """
    Render the detailed (per-cluster) table to a sortable HTML page.

    - Input: the *detailed* CSV (one row per cluster with a `fold_pair` column)
      produced by postprocess_unified.py.
    - Output: an HTML table. First column links to the per-pair HTML page.
    """
    if detailed_csv is None:
        from config import DETAILED_RESULTS_TABLE
        detailed_csv = DETAILED_RESULTS_TABLE
    if output_html is None:
        from config import TABLES_RES
        output_html = os.path.join(TABLES_RES, "clusters_table.html")

    if not os.path.exists(detailed_csv):
        raise FileNotFoundError(f"Detailed CSV not found: {detailed_csv}")

    df = pd.read_csv(detailed_csv)
    if df.empty:
        os.makedirs(os.path.dirname(output_html), exist_ok=True)
        with open(output_html, "w", encoding="utf-8") as f:
            f.write(f"<!DOCTYPE html><html><head><meta charset='utf-8'><title>{html.escape(title)}</title></head>"
                    f"<body><h2>{html.escape(title)}</h2><p>No data available.</p></body></html>")
        return output_html

    # Ensure the pair column exists and is first
    pair_col = "fold_pair" if "fold_pair" in df.columns else ("pair_id" if "pair_id" in df.columns else None)
    if pair_col is None or pair_col not in df.columns:
        raise KeyError("Expected a 'fold_pair' (or 'pair_id') column in the detailed CSV.")
    # Reorder: put pair first; then the rest in original order
    cols = [pair_col] + [c for c in df.columns if c != pair_col]
    df = df[cols]

    # Header
    thead = "<tr>" + "".join(
        f'<th onclick="sortTable({i})">{html.escape(col)}</th>'
        for i, col in enumerate(df.columns)
    ) + "</tr>"

    num_re = re.compile(r"^\s*([+-]?\d+(?:\.\d+)?)")
    def numeric_part(cell) -> str:
        if pd.isna(cell):
            return ""
        s = str(cell).strip()
        if s == "-":
            return ""
        m = num_re.match(s)
        return m.group(1) if m else s

    # Rows (first col is a link to per-pair page)
    rows = []
    for _, r in df.iterrows():
        pair = str(r[pair_col])
        link = (base_pair_url or "{pair_id}.html").format(pair_id=html.escape(pair))
        tds = [f'<td><a href="{link}" target="_blank">{html.escape(pair)}</a></td>']
        for col in df.columns[1:]:
            val = r[col]
            disp = "-" if pd.isna(val) else str(val)
            tds.append(f'<td data-sort-value="{html.escape(numeric_part(val))}">{html.escape(disp)}</td>')
        rows.append("<tr>" + "".join(tds) + "</tr>")

    # Same CSS/JS as the other table
    html_doc = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8" />
<title>{html.escape(title)}</title>
<style>
  body {{ font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; background-color: #121212; color: #E0E0E0; margin: 20px; }}
  table {{ width: 95%; margin: auto; border-collapse: collapse; box-shadow: 0 4px 8px rgba(0,0,0,0.5); }}
  th, td {{ border: 1px solid #333; padding: 10px 14px; text-align: left; white-space: nowrap; }}
  th {{ background-color: #004d40; color: #fff; font-size: 16px; cursor: pointer; position: sticky; top: 0; }}
  tr:nth-child(even) {{ background-color: #2b2b2b; }}
  tr:hover {{ background-color: #3a3a3a; }}
  a {{ color: #BB86FC; text-decoration: none; }}
  a:hover {{ text-decoration: underline; }}
  h2 {{ text-align: center; color: #E0E0E0; margin-top: 0; }}
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


if __name__ == "__main__":
    # 1) Best-per-pair summary table
    out1 = os.path.join(TABLES_RES, "protein_comparison_table.html")
    gen_html_from_summary_table(
        preferred_column_order=[
            "pair_id",
            "BEST_AF_TM_FOLD1",
            "BEST_CMAP_T1_F1",
            "BEST_CMAP_T1_PRECISION",
            "BEST_CMAP_T1_RECALL",
            "BEST_CMAP_T1_JACCARD",
            "BEST_CMAP_T1_MCC",
        ],
        output_html=out1,
    )

    # 2) Cluster-level (detailed) table
    out2 = os.path.join(TABLES_RES, "protein_clusters_table.html")
    gen_html_from_cluster_detailed_table(output_html=out2)

    print("OK:\n ", out1, "\n ", out2)

