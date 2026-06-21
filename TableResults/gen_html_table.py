# File: TableResults/gen_html_table.py
import os, sys, re, html, shutil
import pandas as pd
import argparse

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, ROOT)
from config import *
# NOTE: build_unified_tables_from_cluster_dfs is imported lazily inside
# _ensure_unified_csvs() so that the lightweight HTML/table builders here can be
# imported and run without pulling in torch (only needed to REBUILD the CSVs).

# --- add near the top (after imports / _ensure_unified_csvs) ---
NUM_RE = re.compile(r"^\s*([+-]?\d+(?:\.\d+)?)")
TRIPLET_COL = "MSA: DEPTH; #RES; #Clusters"


# Single source of truth for the main comparison table
DEFAULT_PREFERRED_COLS = [
    "pair_id", "NAME", "ORGANISM", "trigger_class", "#RES", "MSA: DEPTH; #RES; #Clusters",
    "SEQ_ID%", "PAIR_TM", "SS_SWITCH%", "DDG_bias",
    "AF2Clust_TM1","AF2Clust_TM2","AF2Deep_TM1","AF2Deep_TM2",
    "AF3Clust_TM1","AF3Clust_TM2","AF3Deep_TM1","AF3Deep_TM2",
    "ESM2_TM1","ESM2_TM2","ESM3_TM1","ESM3_TM2",
    "MSATrans_CMAP_PR1","MSATrans_CMAP_PR2","MSATrans_CMAP_RE1","MSATrans_CMAP_RE2",
]

DEFAULT_EXPLANATIONS = {
    "NAME": "Protein name(s) of the fold-switch pair (short common alias where known, e.g. XCL1, otherwise the RCSB molecule description). Source: docs/pair_protein_names.csv (scripts/fetch_pair_protein_names.py).",
    "ORGANISM": "Source organism(s) of the two PDB chains. Source: docs/pair_protein_names.csv.",
    "trigger_class": "Auto-classified fold-switching trigger from PDB metadata: ligand / oligomerization / protein_binding / mutation / equilibrium_or_unknown. Source: docs/triggers_from_pdb.csv (scripts/classify_triggers_from_pdb.py).",
    "SEQ_ID%": "Pairwise sequence identity between the two PDB chains of the pair, computed as best-matching chain-pair local-alignment identity (score over the local alignment divided by the length of the shorter sequence). Robust to terminal truncations. The 'mutation' trigger class is defined as SEQ_ID% < 95%. Source: docs/triggers_from_pdb.csv (scripts/classify_triggers_from_pdb.py).",
    "#RES": "Number of residues in the longer chain of the pair.",
    "MSA: DEPTH; #RES; #Clusters": "DEPTH = number of sequences in DeepMsa.a3m; #RES = alignment width (columns, including gaps); #Clusters = count of ShallowMsa_* a3m files.",
    "PAIR_TM": "TM-score between the two ground-truth folds (max of the two directions).",
    "SS_SWITCH%": "Porter-style secondary-structure switch fraction: fraction of aligned residues where DSSP class transitions between H (helix) and E (strand). H<->C and E<->C transitions are NOT counted. Source: docs/ss_switch_per_pair.csv (scripts/compute_ss_switch_for_pairs.py).",
    "DDG_bias": "ThermoMPNN conformation-biasing energetic fold preference: per-pair mean over MSA clusters of ΣΔΔG(impose cluster sequence on Fold2) − ΣΔΔG(on Fold1), in kcal/mol. >0 = the family's sequences are energetically more compatible with Fold1, <0 = Fold2. Relative measure (each fold referenced to its own native), summed over positions. The fold-switching signal is the cluster-to-cluster SPREAD of this quantity (analyzed separately), not the mean. Source: df_ddg.csv (run_DDG.py). Replaces the earlier raw-Rosetta ΔG1/ΔG2 (unrelaxed whole-assembly REU, not meaningful).",
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
    "CONCORD": "AGGREGATE: cross-method per-cluster sign concordance, sequence-divergence-corrected, averaged over the 28 method-pairs of the 8 point-estimate methods (shuffle null ~0.50; higher = orthogonal methods agree more on each cluster's fold preference). Source: docs/fold_diversity_concordance_corrected_perm_8method_bh.csv.",
    "CONCORD_pBH": "AGGREGATE: Benjamini-Hochberg-adjusted permutation p-value for CONCORD over the full pool of pairs (<=0.05 = significant cross-method concordance). Same source.",
    "MORANS_I": "AGGREGATE (clade-preference): strongest positive Moran's I of per-cluster fold preference across the methods (on the cluster tree, inverse-patristic weights). Effect size of phylogenetic clade structure; ~0 = no clade structure, larger = fold-preferring clusters form clades. Source: docs/phylo_placement_corrected.csv.",
    "MORANS_Inorm": "AGGREGATE (clade-preference): Moran's I normalised by the maximum achievable on each pair's tree (0 = null, 1 = maximal clade segregation). Rank clade structure by this, not by the BH p-value (the test is over-powered at ~100 clusters). Same source.",
    "MORANS_pBH": "AGGREGATE (clade-preference): BH-adjusted permutation p-value for the strongest-method Moran's I. NOTE: at ~100 clusters this test is over-powered, so many pairs pass BH at negligible effect size; interpret with MORANS_I / MORANS_Inorm. Same source.",
    "MORANS_method": "AGGREGATE: which method gives the strongest positive Moran's I for this pair.",
}


# === NEW: ensure unified CSVs exist (or rebuild on demand) ===
def _summary_looks_malformed() -> bool:
    """True if the on-disk summary CSV is not the wide per-pair comparison table.
    Guards against a stale/partial summary (e.g. a single-pair per-cluster file
    left by `postprocess --pairs <one>`): such a file has a 'cluster' column,
    lacks 'PAIR_TM', or covers very few pairs. In any of those cases we must
    rebuild rather than render a broken one-pair table."""
    try:
        if not os.path.exists(SUMMARY_RESULTS_TABLE):
            return True
        head = pd.read_csv(SUMMARY_RESULTS_TABLE, nrows=5)
        cols = set(head.columns)
        if "cluster" in cols or "PAIR_TM" not in cols:
            return True
        n_pairs = pd.read_csv(SUMMARY_RESULTS_TABLE, usecols=[0]).iloc[:, 0].nunique()
        return n_pairs < 50
    except Exception:
        return True


def _ensure_unified_csvs(force_rerun: bool = False):
    need = (force_rerun or (not os.path.exists(SUMMARY_RESULTS_TABLE))
            or (not os.path.exists(DETAILED_RESULTS_TABLE)) or _summary_looks_malformed())
    if need:
        # local import to avoid cycles AND to keep torch out of the import path
        # for the (torch-free) HTML/table builders below.
        try:
            from Analysis.postprocess_unified import build_unified_tables_from_cluster_dfs
            build_unified_tables_from_cluster_dfs(write_out=True)
            print("[gen_html] ensured unified CSVs (summary + detailed).")
        except Exception as e:
            print(f"[gen_html] WARN: could not build unified CSVs automatically: {e}")


def _append_mean_std_row(df: pd.DataFrame, label_col: str, label_text: str = "Averages") -> pd.DataFrame:
    out = df.copy()
    means = {}
    for c in out.columns:
        if c == label_col:
            means[c] = label_text
            continue
        if c == TRIPLET_COL:
            triples = []
            for v in out[c]:
                if pd.isna(v): continue
                parts = [p.strip() for p in str(v).split(";")]
                if len(parts) != 3: continue
                nums = []
                for p in parts:
                    m = NUM_RE.match(p)
                    nums.append(float(m.group(1)) if m else float("nan"))
                if all(pd.notna(nums)):
                    triples.append(nums)
            if triples:
                import numpy as _np
                mu = _np.nanmean(_np.array(triples, dtype=float), axis=0)
                means[c] = f"{mu[0]:.2f}; {mu[1]:.2f}; {mu[2]:.2f}"
            else:
                means[c] = "-"
            continue

        s = out[c]
        if s.dtype.kind in "biufc":
            x = pd.to_numeric(s, errors="coerce")
        else:
            x = pd.to_numeric(s.apply(numeric_part), errors="coerce")
        if x.notna().any():
            mu = float(x.mean())
            sd = float(x.std(ddof=1)) if x.count() > 1 else 0.0
            means[c] = f"{mu:.2f} ({sd:.2f})"
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
    gen_html_from_cluster_detailed_table(clusters_csv=DETAILED_RESULTS_TABLE, output_html=out_clusters)
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

    # Merge in trigger_class + seq_identity from docs/triggers_from_pdb.csv
    # (auto-classified by scripts/classify_triggers_from_pdb.py; seq_identity
    # is the best-matching chain-pair local-alignment identity in [0,1]).
    _trig_csv = os.path.join(os.path.dirname(summary_csv), "triggers_from_pdb.csv")
    if os.path.isfile(_trig_csv):
        try:
            _trig_raw = pd.read_csv(_trig_csv)
            _trig_cols = ["pair_id", "trigger_class"]
            if "seq_identity" in _trig_raw.columns:
                _trig_cols.append("seq_identity")
            _trig = _trig_raw[_trig_cols].drop_duplicates("pair_id")
            _id_col = "pair_id" if "pair_id" in df.columns else ("fold_pair" if "fold_pair" in df.columns else df.columns[0])
            if _id_col != "pair_id":
                df = df.rename(columns={_id_col: "pair_id"})
            df = df.merge(_trig, on="pair_id", how="left")
            # Format seq_identity as percent for display, e.g. "0.983" -> "98.3%"
            if "seq_identity" in df.columns:
                df["SEQ_ID%"] = df["seq_identity"].map(
                    lambda x: f"{float(x) * 100:.1f}%" if pd.notnull(x) else ""
                )
                df = df.drop(columns=["seq_identity"])
        except Exception as e:
            print(f"[gen_html] WARN: failed to merge trigger_class/seq_identity: {e}")

    # Merge in SS-switch fraction from docs/ss_switch_per_pair.csv
    # (computed by scripts/compute_ss_switch_for_pairs.py via DSSP).
    _ss_csv = os.path.join(os.path.dirname(summary_csv), "ss_switch_per_pair.csv")
    if os.path.isfile(_ss_csv):
        try:
            _ss = pd.read_csv(_ss_csv)
            if "pair_id" in _ss.columns and "ss_switch_fraction" in _ss.columns:
                _ss = _ss[["pair_id", "ss_switch_fraction"]].drop_duplicates("pair_id")
                df = df.merge(_ss, on="pair_id", how="left")
                # Format as percent, e.g. "0.27" -> "27.0%"
                df["SS_SWITCH%"] = df["ss_switch_fraction"].map(
                    lambda x: f"{float(x) * 100:.1f}%" if pd.notnull(x) else ""
                )
                df = df.drop(columns=["ss_switch_fraction"])
        except Exception as e:
            print(f"[gen_html] WARN: failed to merge ss_switch_fraction: {e}")

    # --- Merge AGGREGATE cross-method results: cross-method concordance and
    #     phylogenetic clade-preference (Moran's I). These are the per-pair
    #     summary statistics, kept in their own column group on the right. ---
    _docs_dir = os.path.dirname(summary_csv)
    # (1) Corrected full-pool concordance perm-BH (8 point-estimate methods)
    _conc_csv = os.path.join(_docs_dir, "fold_diversity_concordance_corrected_perm_8method_bh.csv")
    if os.path.isfile(_conc_csv):
        try:
            _c = pd.read_csv(_conc_csv)
            _ccol = "observed_mean_concordance_corrected"
            if "pair_id" in _c.columns and _ccol in _c.columns:
                _c = _c[["pair_id", _ccol, "p_bh"]].drop_duplicates("pair_id").rename(
                    columns={_ccol: "CONCORD", "p_bh": "CONCORD_pBH"})
                df = df.merge(_c, on="pair_id", how="left")
        except Exception as e:
            print(f"[gen_html] WARN: failed to merge concordance: {e}")
    # (2) Clade-preference: per-pair strongest positive Moran's I across methods
    _phylo_csv = os.path.join(_docs_dir, "phylo_placement_corrected.csv")
    if os.path.isfile(_phylo_csv):
        try:
            _m = pd.read_csv(_phylo_csv)
            if {"pair_id", "morans_I"}.issubset(_m.columns):
                _pcol = "morans_p_bh" if "morans_p_bh" in _m.columns else None
                _keep = ["pair_id", "morans_I", "morans_I_norm", "method"]
                _keep = [c for c in _keep if c in _m.columns] + ([_pcol] if _pcol else [])
                _mb = (_m[_m["morans_I"].notna()]
                       .sort_values("morans_I", ascending=False)
                       .drop_duplicates("pair_id"))[_keep]
                _ren = {"morans_I": "MORANS_I", "morans_I_norm": "MORANS_Inorm",
                        "method": "MORANS_method"}
                if _pcol: _ren[_pcol] = "MORANS_pBH"
                _mb = _mb.rename(columns=_ren)
                df = df.merge(_mb, on="pair_id", how="left")
        except Exception as e:
            print(f"[gen_html] WARN: failed to merge Moran's I: {e}")


    # Now also add RMSD for structure


    # If SEQ_ID exists from another source, add a formatted column for display
    # (keep raw numeric in df). seq_identity merged from triggers_from_pdb.csv
    # is already formatted as SEQ_ID% above, so we skip if already present.
    if "SEQ_ID" in df.columns and "SEQ_ID%" not in df.columns:
        df["SEQ_ID%"] = df["SEQ_ID"].map(lambda x: f"{x:.1f}%" if pd.notnull(x) else "")

    # Column ordering: insert "SEQ_ID%" right after "trigger_class" (it is a
    # property of the PDB pair, so it lives next to the trigger-class column);
    # fall back to "PAIR_TM" if trigger_class is not in the frame.
    cols = list(df.columns)
    if "SEQ_ID%" in cols:
        cols.remove("SEQ_ID%")
        if "trigger_class" in cols:
            anchor = "trigger_class"
        elif "PAIR_TM" in cols:
            anchor = "PAIR_TM"
        else:
            anchor = None
        if anchor is not None:
            pos = cols.index(anchor)
            cols.insert(pos + 1, "SEQ_ID%")
            df = df[cols]

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

    # Energetic fold preference: surface DDG_Bias_mean as "DDG_bias" and drop
    # the broken raw-Rosetta ΔG1/ΔG2 columns (unrelaxed whole-assembly REU).
    if "DDG_Bias_mean" in df.columns and "DDG_bias" not in df.columns:
        df = df.rename(columns={"DDG_Bias_mean": "DDG_bias"})
    df = df.drop(columns=[c for c in ("ΔG1", "ΔG2") if c in df.columns])

    # Merge human-readable protein NAME + ORGANISM from the canonical table
    # docs/pair_protein_names.csv (scripts/fetch_pair_protein_names.py).
    _names_csv = os.path.join(os.path.dirname(summary_csv), "pair_protein_names.csv")
    if os.path.isfile(_names_csv):
        try:
            _nm = pd.read_csv(_names_csv)[["pair_id", "display_name", "organism"]]
            _nm = _nm.drop_duplicates("pair_id").rename(
                columns={"display_name": "NAME", "organism": "ORGANISM"})
            df = df.merge(_nm, on="pair_id", how="left")
        except Exception as e:
            print(f"[gen_html] WARN: failed to merge protein names: {e}")

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

    # Force the aggregate (cross-method) columns to the very end so they form a
    # single right-hand column group.
    _AGG_COLS = ["CONCORD", "CONCORD_pBH", "MORANS_I", "MORANS_Inorm",
                 "MORANS_pBH", "MORANS_method"]
    _agg_present = [c for c in _AGG_COLS if c in ordered]
    ordered = [c for c in ordered if c not in _agg_present] + _agg_present

    df = df[ordered]

    # --- Column-group bookkeeping (for coloured headers + thick separators) ---
    # Three groups: description / single-method / aggregate. The first column of
    # each group gets a thick left border; NAME/ORGANISM wrap to two lines.
    _DESC_COLS = {"pair_id", "NAME", "ORGANISM", "trigger_class", "#RES",
                  "MSA: DEPTH; #RES; #Clusters", "SEQ_ID%", "PAIR_TM", "SS_SWITCH%"}
    _AGG_SET = set(_AGG_COLS)
    _WRAP_COLS = {"NAME", "ORGANISM"}

    def _col_group(col: str) -> str:
        if col in _AGG_SET:
            return "agg"
        if col in _DESC_COLS:
            return "desc"
        return "method"

    _group_start = set()
    _prev_g = None
    for _c in df.columns:
        _g = _col_group(_c)
        if _g != _prev_g:
            _group_start.add(_c)
        _prev_g = _g

    def _th_cls(col: str) -> str:
        cl = [f"g-{_col_group(col)}"]
        if col in _group_start:
            cl.append("gstart")
        if col in _WRAP_COLS:
            cl.append("wrapcol")
        return ' class="' + " ".join(cl) + '"'

    def _td_cls(col: str) -> str:
        cl = []
        if col in _group_start:
            cl.append("gstart")
        if col in _WRAP_COLS:
            cl.append("wrapcol")
        return (' class="' + " ".join(cl) + '"') if cl else ""

    # Round numeric columns to 3 decimals for display. Skip pre-formatted
    # percent columns (already 1-decimal, e.g. SEQ_ID%, SS_SWITCH%) and text
    # columns (NAME, ORGANISM) so they are not re-rounded to 3 decimals.
    _no_round = {"pair_id", "NAME", "ORGANISM", "trigger_class", "SEQ_ID%", "SS_SWITCH%"}
    for c in df.columns:
        if c in _no_round: continue
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


    # --- Table header (clickable for sorting) ---
    # Display-name overrides: underlying CSV column names are kept for
    # backward compatibility, but the header rendered to HTML uses the
    # Greek symbols for ΔG and ΔΔG so the table matches the paper's
    # notation. Anywhere "DDG" appears in a CSV column name is rewritten
    # to "ΔΔG" in the display header.
    def _display_header(col: str) -> str:
        return col.replace("DDG", "ΔΔG")

    thead = "<tr>" + "".join(
        f'<th{_th_cls(col)} onclick="sortTable({i})">{html.escape(_display_header(col))}</th>'
        for i, col in enumerate(df.columns)
    ) + "</tr>"

    def _format_triplet_cell(val, is_avg: bool) -> str:
        if pd.isna(val):
            return "-"
        s = str(val)
        parts = [p.strip() for p in s.split(";")]
        out = []
        for p in parts:
            m = NUM_RE.match(p)
            if not m:
                out.append(p)
                continue
            x = float(m.group(1))
            out.append(f"{x:.2f}" if is_avg else str(int(round(x))))
        return "; ".join(out)

    # For sorting: pull numeric prefix from strings like "3219 (18)" -> "3219"

    # --- Rows. When building rows: choose row class based on cluster count ---
    rows = []
    for _, r in df.iterrows():
        pair = str(r["pair_id"])
        link = (base_pair_url or "{pair_id}.html").format(pair_id=html.escape(pair))
        tds = [f'<td{_td_cls("pair_id")}><a href="{link}" target="_blank">{html.escape(pair)}</a></td>']
        is_avg = pair.lower().startswith("average")  # "Average" or "Averages"

        for col in df.columns[1:]:
            val = r[col]

            if col == TRIPLET_COL:
                disp = _format_triplet_cell(val, is_avg)
            elif col == "#RES" and not is_avg and pd.notna(val):
                # show integers for #RES in normal rows
                try:
                    disp = str(int(round(float(val))))
                except Exception:
                    disp = "-" if pd.isna(val) else str(val)
            else:
                disp = "-" if pd.isna(val) else str(val)

            tds.append(f'<td{_td_cls(col)} data-sort-value="{html.escape(numeric_part(val))}">{html.escape(disp)}</td>')

        row_cls = ""
        if fade_min_clusters is not None:
            ncl = _clusters_in_row(r)
            if ncl is not None and ncl < fade_min_clusters:
                row_cls = ' class="lowclust"'
        rows.append(f"<tr{row_cls}>" + "".join(tds) + "</tr>")


    # ---- Averages row (2 decimals). Special handling for "MSA: DEPTH; #RES; #Clusters" ----
    avg_label = '<a href="https://orzuk.github.io/MsaCluster/pairs_global_analysis.html" target="_blank">Averages</a>'
    avg_tds = [f'<td{_td_cls("pair_id")}>{avg_label}</td>']

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
        triple_disp = f"{int(md)} ; {mr:.2f} ; {mc:.2f}"

    # Build cells
    for col in df.columns[1:]:
        if col == triple_col:
            avg_tds.append(f'<td{_td_cls(col)} data-sort-value="">{html.escape(triple_disp)}</td>')
            continue
        s = pd.to_numeric(df[col].map(numeric_part), errors="coerce")
        if s.notna().any():
            mu = s.mean()
            sd = s.std(ddof=1) if s.count() > 1 else 0.0
            disp = f"{mu:.2f} ({sd:.2f})"
            avg_tds.append(f'<td{_td_cls(col)} data-sort-value="{mu:.6f}">{html.escape(disp)}</td>')
        else:
            avg_tds.append(f'<td{_td_cls(col)} data-sort-value=""></td>')
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
  th, td {{ border: 1px solid #333; padding: 10px 14px; text-align: left; white-space: nowrap; vertical-align: top; }}
  th {{ background-color: #b71c1c; color: #fff; font-size: 16px; cursor: pointer; position: sticky; top: 0; }}
  /* column-group header colours: description (red) / single-method (blue) / aggregate (green) */
  th.g-desc {{ background-color: #b71c1c; }}
  th.g-method {{ background-color: #1565c0; }}
  th.g-agg {{ background-color: #2e7d32; }}
  /* thick separator at the first column of each group; thin (#333) within a group */
  th.gstart, td.gstart {{ border-left: 3px solid #cfcfcf; }}
  /* let NAME / ORGANISM wrap to ~two lines instead of forcing one very wide column */
  th.wrapcol, td.wrapcol {{ white-space: normal; max-width: 16ch; }}
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

    # Put the pair id column first for display & links
    pair_col = "fold_pair" if "fold_pair" in df.columns else ("pair_id" if "pair_id" in df.columns else None)
    if pair_col is None:
        raise KeyError("Expected a 'fold_pair' (or 'pair_id') column in the detailed CSV.")
    cols = [pair_col] + [c for c in df.columns if c != pair_col]
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
    images_dir: str = None,
    output_html: str = os.path.join("docs", "pairs_global_analysis.html"),
    title: str = "Global Analysis — Fold-Switch Evolution"
) -> str:
    """Render the global analysis page from the ordered, captioned figure spec.

    The spec (utils.plot_global.GLOBAL_FIGURE_SPEC) is the single place that
    defines which global figures appear, in what order, and the explanatory
    paragraph shown before each. `images_dir` is accepted for backward
    compatibility but ignored — image paths come from the spec (relative to the
    page's own directory). Spec entries whose file is missing are skipped.
    """
    from utils.plot_global import write_global_html
    return write_global_html(output_html=output_html, title=title)


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
