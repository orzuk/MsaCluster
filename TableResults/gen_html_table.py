# File: TableResults/gen_html_table.py
import os, sys, re, html
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, ROOT)
from config import TABLES_RES, SUMMARY_RESULTS_TABLE, DETAILED_RESULTS_TABLE, GITHUB_URL_HTML

def gen_html_from_summary_table(
    summary_csv: str | None = None,
    output_html: str | None = None,
    title: str = "Interactive Protein Comparison Table",
    base_pair_url: str | None = GITHUB_URL_HTML + "/{pair_id}.html",
    preferred_column_order: list[str] | None = None,
    column_explanations: dict[str, str] | None = None,
) -> str:
    if summary_csv is None:
        summary_csv = SUMMARY_RESULTS_TABLE
    if output_html is None:
        output_html = os.path.join(TABLES_RES, "table.html")

    if not os.path.exists(summary_csv):
        raise FileNotFoundError(f"Summary CSV not found: {summary_csv}")

    df = pd.read_csv(summary_csv)

    if df.empty:
        os.makedirs(os.path.dirname(output_html), exist_ok=True)
        with open(output_html, "w", encoding="utf-8") as f:
            f.write(f"<!DOCTYPE html><html><head><meta charset='utf-8'><title>{html.escape(title)}</title></head>"
                    f"<body><h2>{html.escape(title)}</h2><p>No data available.</p></body></html>")
        return output_html

    pair_col = "pair_id" if "pair_id" in df.columns else "fold_pair"
    cols = list(df.columns)

    if preferred_column_order is not None:
        wanted = [c for c in preferred_column_order if c in cols]
        rest = [c for c in cols if c not in wanted]
        df = df[wanted + rest]
    else:
        df = df[[pair_col] + [c for c in cols if c != pair_col]]

    # Build header
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

    # Build explanations block (only for columns we actually have)
    expl_lines = []
    if column_explanations:
        for c in df.columns:
            if c in column_explanations:
                expl_lines.append(f"<p><b>{html.escape(c)}</b>: {html.escape(column_explanations[c])}</p>")
    expl_html = ("\n".join(expl_lines)) if expl_lines else ""

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
    detailed_csv: str | None = None,
    output_html: str | None = None,
    title: str = "Cluster-level Results (one row per cluster)",
    base_pair_url: str | None = GITHUB_URL_HTML + "/{pair_id}.html",
) -> str:
    if detailed_csv is None:
        detailed_csv = DETAILED_RESULTS_TABLE
    if output_html is None:
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

    pair_col = "fold_pair" if "fold_pair" in df.columns else ("pair_id" if "pair_id" in df.columns else None)
    if pair_col is None:
        raise KeyError("Expected a 'fold_pair' (or 'pair_id') column in the detailed CSV.")
    cols = [pair_col] + [c for c in df.columns if c != pair_col]
    df = df[cols]

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
    # Your main comparison table (to docs/protein_comparison_table.html)
    out1 = os.path.join(TABLES_RES, "protein_comparison_table.html")
    preferred = [
    "pair_id", "#RES",
    "AF2Clust_TM1","AF2Clust_TM2","AF2Deep_TM1","AF2Deep_TM2",
    "AF3Clust_TM1","AF3Clust_TM2","AF3Deep_TM1","AF3Deep_TM2",
    "ESM2_TM1","ESM2_TM2","ESM3_TM1","ESM3_TM2",
    "MSATrans_CMAP_PR1","MSATrans_CMAP_PR2","MSATrans_CMAP_RE1","MSATrans_CMAP_RE2",
    ]
    explanations = {
        "#RES": "Number of residues in the longer chain of the pair.",
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
    gen_html_from_summary_table(
        preferred_column_order=preferred,
        output_html=out1,
        column_explanations=explanations,
    )

    # Cluster-level (detailed) table unchanged
    out2 = os.path.join(TABLES_RES, "protein_clusters_table.html")
    gen_html_from_cluster_detailed_table(output_html=out2)

    print("OK:\n ", out1, "\n ", out2)
