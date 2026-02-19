"""Cross-pair global statistics scatter plots."""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from config import DATA_DIR, FIGURE_RES_DIR, DETAILED_RESULTS_TABLE, SUMMARY_RESULTS_TABLE
from utils.utils import list_protein_pairs, numify


def global_pairs_statistics_plots(output_dir: str | None = None) -> None:
    """
    Build ALL global scatter plots from the unified CSV(s) and save them under output_dir.
    """
    if output_dir is None:
        output_dir = FIGURE_RES_DIR
    os.makedirs(output_dir, exist_ok=True)

    # -------- AF/ESM from unified detailed table --------
    if (not os.path.isfile(DETAILED_RESULTS_TABLE)) or os.path.getsize(DETAILED_RESULTS_TABLE) == 0:
        print(f"[warn] No detailed CSV at {DETAILED_RESULTS_TABLE}. Run postprocess first.")
        df_det = pd.DataFrame(columns=["fold_pair","model","TMscore_fold1","TMscore_fold2"])
    else:
        df_det = pd.read_csv(DETAILED_RESULTS_TABLE)

    if "TMscore_fold1" not in df_det.columns and "score_pdb1" in df_det.columns:
        df_det = df_det.rename(columns={"score_pdb1": "TMscore_fold1", "score_pdb2": "TMscore_fold2"})
    if "fold_pair" not in df_det.columns:
        if "pair_id" in df_det.columns:
            df_det = df_det.rename(columns={"pair_id": "fold_pair"})
        else:
            df_det["fold_pair"] = "UNK"

    if "model" not in df_det.columns:
        df_det["model"] = ""
    df_det["model"] = df_det["model"].astype(str).str.upper()

    def _save_tm_scatter(sub_df: pd.DataFrame, model_tag: str, out_path: str) -> None:
        if sub_df.empty:
            print(f"[info] No rows for {model_tag} in {DETAILED_RESULTS_TABLE}")
            return

        g = sub_df.groupby("fold_pair", as_index=False).agg(
            mean_TM1=("TMscore_fold1","mean"),
            mean_TM2=("TMscore_fold2","mean"),
            max_TM1 =("TMscore_fold1","max"),
            max_TM2 =("TMscore_fold2","max"),
        )

        plt.figure(figsize=(8,7))
        first = True
        for _, r in g.iterrows():
            plt.scatter(r["mean_TM1"], r["mean_TM2"], marker="+", label="Mean" if first else "")
            plt.scatter(r["max_TM1"],  r["max_TM2"],  marker="*", label="Max"  if first else "")
            plt.plot([r["mean_TM1"], r["max_TM1"]], [r["mean_TM2"], r["max_TM2"]],
                     linestyle="--", color="gray", alpha=0.5)
            plt.text(r["max_TM1"], r["max_TM2"], r["fold_pair"], fontsize=7, ha="right")
            first = False

        plt.xlabel("TMscore fold1")
        plt.ylabel("TMscore fold2")
        plt.title(f"Per-pair TM (mean vs max): {model_tag}")
        plt.legend(loc="upper left")
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(out_path, dpi=160)
        plt.close()
        print(f"[plot] wrote {out_path}")

    for tag in ["AF2", "AF3", "ESM2", "ESM3"]:
        _save_tm_scatter(
            df_det[df_det["model"] == tag],
            model_tag=tag,
            out_path=os.path.join(output_dir, f"fold_pair_scatter_plot_{tag}.png"),
        )

    # -------- MSA-Transformer (recall) from per-pair df_cmap.csv --------
    pairs = list_protein_pairs(parsed=False)
    rows = []
    for pid in pairs:
        p_csv = os.path.join(DATA_DIR, pid, "Analysis", "df_cmap.csv")
        if not (os.path.isfile(p_csv) and os.path.getsize(p_csv) > 0):
            continue
        try:
            d = pd.read_csv(p_csv)
        except Exception:
            continue

        c1 = "t1_recall" if "t1_recall" in d.columns else None
        c2 = "t2_recall" if "t2_recall" in d.columns else None
        if not c1 or not c2:
            continue

        rows.append({
            "fold_pair": pid,
            "max_R1": pd.to_numeric(d[c1], errors="coerce").max(),
            "max_R2": pd.to_numeric(d[c2], errors="coerce").max(),
        })

    if rows:
        df_msat = pd.DataFrame(rows)
        out = os.path.join(output_dir, "fold_pair_scatter_plot_MSA_TRANS.png")
        plt.figure(figsize=(8,7))
        first = True
        for _, r in df_msat.iterrows():
            plt.scatter(r["max_R1"], r["max_R2"], marker="*", label="Max" if first else "")
            plt.text(r["max_R1"], r["max_R2"], r["fold_pair"], fontsize=7, ha="right")
            first = False
        plt.xlabel("MSA-Transformer Recall (fold1)")
        plt.ylabel("MSA-Transformer Recall (fold2)")
        plt.title("Per-pair max recall: MSA-Transformer")
        plt.legend(loc="upper left")
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(out, dpi=160)
        plt.close()
        print(f"[plot] wrote {out}")
    else:
        print("[info] No df_cmap.csv files found; skipping MSA-Transformer plot.")

    # === Cross-method pairwise comparisons ===
    if (not os.path.isfile(SUMMARY_RESULTS_TABLE)) or os.path.getsize(SUMMARY_RESULTS_TABLE) == 0:
        print(f"[warn] No summary CSV at {SUMMARY_RESULTS_TABLE}. Skipping cross-method plots.")
        return

    df_sum = pd.read_csv(SUMMARY_RESULTS_TABLE)

    cols_deep = dict(
        AF2_TM1="AF2Deep_TM1", AF2_TM2="AF2Deep_TM2",
        AF3_TM1="AF3Deep_TM1", AF3_TM2="AF3Deep_TM2",
    )
    cols_bestcl = dict(
        AF2_TM1="AF2Clust_TM1", AF2_TM2="AF2Clust_TM2",
        AF3_TM1="AF3Clust_TM1", AF3_TM2="AF3Clust_TM2",
    )
    cols_esm = dict(
        ESM2_TM1="ESM2_TM1", ESM2_TM2="ESM2_TM2",
        ESM3_TM1="ESM3_TM1", ESM3_TM2="ESM3_TM2",
    )

    def _split_pair_id(pid: str):
        s = str(pid)
        return tuple(s.split("_", 1)) if "_" in s else (s, "fold2")

    def _plot_xy_pairs(
        df: pd.DataFrame,
        x1_col: str, x2_col: str,
        y1_col: str, y2_col: str,
        title: str,
        outfile: str
    ):
        sub = df.copy()
        for c in [x1_col, x2_col, y1_col, y2_col]:
            if c not in sub.columns:
                print(f"[skip] Missing column {c} for {title}")
                return
            sub[c] = sub[c].map(numify)
        sub = sub.dropna(subset=[x1_col, x2_col, y1_col, y2_col])
        if sub.empty:
            print(f"[info] No rows for plot: {title}")
            return

        plt.figure(figsize=(6.4, 6.0))
        plt.plot([0,1],[0,1], 'k-', linewidth=1, alpha=0.8)

        for _, r in sub.iterrows():
            pid = r.get("pair_id") or r.get("fold_pair") or "pair"
            f1, f2 = _split_pair_id(pid)

            x1, y1 = r[x1_col], r[y1_col]
            x2, y2 = r[x2_col], r[y2_col]

            plt.plot([x1, x2], [y1, y2], linestyle='--', linewidth=0.8, alpha=0.5)
            plt.scatter(x1, y1, s=22, marker='o')
            plt.scatter(x2, y2, s=22, marker='o')

            try:
                plt.text(x1, y1, f1, fontsize=6, ha='right', va='bottom', alpha=0.9)
                plt.text(x2, y2, f2, fontsize=6, ha='left', va='top', alpha=0.9)
            except Exception:
                pass

        plt.xlim(0, 1); plt.ylim(0, 1)
        plt.xlabel("X method TM-score")
        plt.ylabel("Y method TM-score")
        plt.title(title)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        out_path = os.path.join(output_dir, outfile)
        plt.savefig(out_path, dpi=180)
        plt.close()
        print(f"[plot] wrote {out_path}")

    _plot_xy_pairs(
        df_sum,
        x1_col=cols_deep["AF2_TM1"], x2_col=cols_deep["AF2_TM2"],
        y1_col=cols_deep["AF3_TM1"], y2_col=cols_deep["AF3_TM2"],
        title="AF2 (Deep) vs AF3 (Deep) — TM-scores (per fold, pairs linked)",
        outfile="compare_AF2_AF3_DEEP.png"
    )

    _plot_xy_pairs(
        df_sum,
        x1_col=cols_bestcl["AF2_TM1"], x2_col=cols_bestcl["AF2_TM2"],
        y1_col=cols_bestcl["AF3_TM1"], y2_col=cols_bestcl["AF3_TM2"],
        title="AF2 (Best Shallow Cluster) vs AF3 (Best Shallow Cluster)",
        outfile="compare_AF2_AF3_BESTCLUST.png"
    )

    _plot_xy_pairs(
        df_sum,
        x1_col=cols_esm["ESM2_TM1"], x2_col=cols_esm["ESM2_TM2"],
        y1_col=cols_esm["ESM3_TM1"], y2_col=cols_esm["ESM3_TM2"],
        title="ESM2 (Best Prediction) vs ESM3 (Best Prediction)",
        outfile="compare_ESM2_ESM3_BESTPRED.png"
    )
