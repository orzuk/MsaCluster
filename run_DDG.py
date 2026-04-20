#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_DDG.py — per-cluster conformation-biasing (ThermoMPNN) scoring.

Scores cluster homolog sequences against both fold backbones using
ThermoMPNN's DDG prediction, following the contrastive scoring idea of
Cavanagh et al. (Science 2026, conformation-biasing / "CB"):

    Bias(s) = sum_i [DDG_F2(i, s_i)] - sum_i [DDG_F1(i, s_i)]

    Positive  -> sequence s prefers fold 1 (destabilized more on F2).
    Negative  -> sequence s prefers fold 2.

We only call ThermoMPNN twice per pair (once per backbone), caching the
full L x 19 DDG matrix; per-sequence scoring is a cheap lookup.

Inputs
  Pipeline/<pair>/output_get_msa/DeepMsa.a3m        (for seed residue indices)
  Pipeline/<pair>/output_msa_cluster/ShallowMsa_*.a3m (clusters)
  Pipeline/<pair>/chain_pdb_files/{tag}.pdb          (truth backbones; falls back to whole PDB)

Outputs
  Pipeline/<pair>/output_ddg/ddg_matrix_{tag1}.npz   (cached L x 20 DDG, columns = AA_LIST)
  Pipeline/<pair>/output_ddg/ddg_matrix_{tag2}.npz
  Pipeline/<pair>/output_ddg/cluster_{NNN}.json      (per-sequence + per-residue detail)
  Pipeline/<pair>/Analysis/df_ddg.csv                (per-cluster summary)

Usage
  python3 run_DDG.py -input 4dxrA_4dxtA --thermompnn_dir $THERMOMPNN_DIR
"""
from __future__ import annotations
import os
import sys
import json
import time
import argparse
from pathlib import Path
from typing import Dict, List, Optional, Tuple  # noqa: F401

import numpy as np
import pandas as pd

# Project config & utils
from config import DATA_DIR, MAIN_DIR  # noqa: F401
from utils.protein_utils import read_msa, truth_pdbs, greedy_select
from utils.utils import pair_str_to_tuple
from Analysis.cmap_analysis import seed_residue_indices


# ----------------------------------------------------------------------
# Amino acid constants (must match ThermoMPNN / ProteinMPNN ordering)
# ----------------------------------------------------------------------
AA_LIST = list("ACDEFGHIKLMNPQRSTVWY")
AA_TO_IDX = {aa: i for i, aa in enumerate(AA_LIST)}
N_AA = len(AA_LIST)


# ======================================================================
# ThermoMPNN scorer (adapted from robustness-dynamics/compute_robustness.py)
# ======================================================================

class ThermoMPNNScorer:
    """Compute a full L x 20 DDG matrix for a backbone.

    Entry at (i, a) is the predicted DDG (kcal/mol) of mutating position i
    from the PDB wild-type residue to amino acid AA_LIST[a]. Wild-type
    entries (a == wt_idx) are 0.0 by construction.
    """

    def __init__(self, thermompnn_dir: Optional[str] = None,
                 checkpoint: str = "thermoMPNN_default.pt",
                 batch_size: int = 1000):
        self._model = None
        self._device = "cuda"
        self._thermompnn_dir = thermompnn_dir or os.environ.get(
            "THERMOMPNN_DIR", None)
        self._checkpoint = checkpoint
        self._batch_size = batch_size
        self._alt_parse_PDB = None

    def load(self, device: str = "cuda") -> None:
        self._device = device
        if self._thermompnn_dir is None:
            raise ValueError(
                "ThermoMPNN directory not set. Either set THERMOMPNN_DIR env "
                "var or pass --thermompnn_dir.")
        thermompnn_dir = Path(self._thermompnn_dir)
        repo_str = str(thermompnn_dir)
        if repo_str not in sys.path:
            sys.path.insert(0, repo_str)

        import torch
        from omegaconf import OmegaConf
        from train_thermompnn import TransferModelPL
        from protein_mpnn_utils import alt_parse_PDB
        self._alt_parse_PDB = alt_parse_PDB

        cfg_dict = {
            "training": {"num_workers": 1, "learn_rate": 1e-3, "epochs": 100, "lr_schedule": True},
            "model": {"hidden_dims": [64, 32], "subtract_mut": True,
                      "num_final_layers": 2, "freeze_weights": True,
                      "load_pretrained": True, "lightattn": True},
            "platform": {"thermompnn_dir": str(thermompnn_dir)},
        }
        local_yaml = thermompnn_dir / "local.yaml"
        if local_yaml.exists():
            cfg = OmegaConf.merge(OmegaConf.load(str(local_yaml)), cfg_dict)
        else:
            cfg = OmegaConf.create(cfg_dict)

        ckpt_path = thermompnn_dir / "models" / self._checkpoint
        if not ckpt_path.exists():
            raise FileNotFoundError(f"ThermoMPNN checkpoint not found: {ckpt_path}")

        self._model = TransferModelPL.load_from_checkpoint(str(ckpt_path), cfg=cfg).model
        self._model = self._model.eval().to(device)
        print(f"[ddg] ThermoMPNN loaded from {ckpt_path}", flush=True)

    def compute_matrix(self, pdb_path: str, chain_id: str) -> Tuple[np.ndarray, str]:
        """Return (L x N_AA float32 matrix, parsed_sequence).

        Uses the PDB's own sequence as wild-type; column order = AA_LIST.
        """
        import torch
        from datasets import Mutation

        pdb_data = self._alt_parse_PDB(pdb_path, input_chain_list=chain_id)
        if not pdb_data or not pdb_data[0].get("seq", ""):
            pdb_data = self._alt_parse_PDB(pdb_path)
        if not pdb_data or not pdb_data[0].get("seq", ""):
            raise RuntimeError(f"alt_parse_PDB returned empty sequence for {pdb_path} chain={chain_id}")

        seq = pdb_data[0]["seq"]
        L = len(seq)
        ddg = np.zeros((L, N_AA), dtype=np.float32)

        all_muts, all_idx = [], []
        for i in range(L):
            wt_aa = seq[i]
            if wt_aa not in AA_TO_IDX:
                continue
            for j, aa in enumerate(AA_LIST):
                if aa == wt_aa:
                    continue  # WT stays 0.0
                all_muts.append(Mutation(position=i, wildtype=wt_aa, mutation=aa,
                                         ddG=None, pdb=pdb_data[0].get("name", "protein")))
                all_idx.append((i, j))

        with torch.no_grad():
            for start in range(0, len(all_muts), self._batch_size):
                end = min(start + self._batch_size, len(all_muts))
                preds, _ = self._model(pdb_data, all_muts[start:end])
                for (i, j), pred in zip(all_idx[start:end], preds):
                    if pred is not None and "ddG" in pred:
                        ddg[i, j] = float(pred["ddG"].cpu().item())
        return ddg, seq


# ======================================================================
# A3M / alignment helpers
# ======================================================================

def _strip_insertions(a3m_seq: str) -> str:
    """Drop lowercase insertions, keep uppercase residues and '-' gaps.
    Returns the match-state projection of an A3M sequence.
    """
    return "".join(c for c in a3m_seq if c.isupper() or c == "-")


# ======================================================================
# Main per-pair driver
# ======================================================================

def _cluster_a3ms(pair_id: str) -> List[Path]:
    d = Path(DATA_DIR) / pair_id / "output_msa_cluster"
    return sorted(d.glob("ShallowMsa_*.a3m"))


def _deep_a3m_path(pair_id: str) -> Path:
    return Path(DATA_DIR) / pair_id / "output_get_msa" / "DeepMsa.a3m"


def _ensure_outdir(pair_id: str) -> Path:
    out = Path(DATA_DIR) / pair_id / "output_ddg"
    out.mkdir(parents=True, exist_ok=True)
    return out


def _analysis_dir(pair_id: str) -> Path:
    out = Path(DATA_DIR) / pair_id / "Analysis"
    out.mkdir(parents=True, exist_ok=True)
    return out


def _score_sequence(match_seq: str,
                    ddg_F1: np.ndarray, pdb_seq_F1: str, idx_F1: np.ndarray,
                    ddg_F2: np.ndarray, pdb_seq_F2: str, idx_F2: np.ndarray
                    ) -> Dict:
    """Score one match-state sequence against both folds.

    One global alignment: idx_F1 maps chain-F1's k-th PDB residue to a
    match-state column; same for F2. Per-residue contribution is DDG at
    the observed homolog residue; '-' or non-standard residues contribute 0.
    """
    def score_one(ddg: np.ndarray, cols: np.ndarray):
        L = len(cols)
        contribs = np.zeros(L, dtype=np.float32)
        for k in range(L):
            aa = match_seq[cols[k]]
            if aa == "-" or aa not in AA_TO_IDX:
                continue
            contribs[k] = ddg[k, AA_TO_IDX[aa]]
        return contribs

    if ddg_F1.shape[0] != len(idx_F1):
        raise ValueError(f"F1 length mismatch: PDB={ddg_F1.shape[0]} vs seed idx={len(idx_F1)}")
    if ddg_F2.shape[0] != len(idx_F2):
        raise ValueError(f"F2 length mismatch: PDB={ddg_F2.shape[0]} vs seed idx={len(idx_F2)}")

    c1 = score_one(ddg_F1, idx_F1)
    c2 = score_one(ddg_F2, idx_F2)
    return {
        "contribs_F1": c1.tolist(),
        "contribs_F2": c2.tolist(),
        "sum_F1": float(c1.sum()),
        "sum_F2": float(c2.sum()),
        "bias": float(c2.sum() - c1.sum()),  # >0 => prefers F1 (more destabilized on F2)
    }


def run_pair(pair_id: str, thermompnn_dir: str, device: str = "cuda",
             cluster_sample_n: int = 10, batch_size: int = 1000,
             force: bool = False) -> None:
    t_start = time.time()
    outdir = _ensure_outdir(pair_id)
    analysis_dir = _analysis_dir(pair_id)

    # Pair tags & PDBs
    tagA, tagB = pair_str_to_tuple(pair_id)       # e.g. "4dxrA", "4dxtA"
    pdb1, c1, pdb2, c2 = truth_pdbs(pair_id)
    deep_a3m = _deep_a3m_path(pair_id)
    if not deep_a3m.is_file():
        raise FileNotFoundError(f"DeepMsa.a3m not found: {deep_a3m}")

    # Residue indices per chain in the match-state frame (uses DeepMsa seeds)
    idx_F1 = seed_residue_indices(str(deep_a3m), tagA)
    idx_F2 = seed_residue_indices(str(deep_a3m), tagB)

    # -------- ThermoMPNN: two matrix computations (cached on disk) -------
    scorer: Optional[ThermoMPNNScorer] = None

    def _get_matrix(tag: str, pdb: str, chain: str) -> Tuple[np.ndarray, str]:
        cache = outdir / f"ddg_matrix_{tag}.npz"
        if cache.exists() and not force:
            z = np.load(cache, allow_pickle=False)
            print(f"[ddg] cached matrix for {tag}: shape={z['ddg'].shape}", flush=True)
            return z["ddg"], str(z["seq"])
        nonlocal scorer
        if scorer is None:
            scorer = ThermoMPNNScorer(thermompnn_dir=thermompnn_dir, batch_size=batch_size)
            scorer.load(device=device)
        t0 = time.time()
        m, seq = scorer.compute_matrix(pdb, chain)
        print(f"[ddg] {tag}: L={m.shape[0]}  ({time.time()-t0:.1f}s)", flush=True)
        np.savez(cache, ddg=m, seq=seq, aa_list="".join(AA_LIST))
        return m, seq

    ddg1, seq1 = _get_matrix(tagA, pdb1, c1)
    ddg2, seq2 = _get_matrix(tagB, pdb2, c2)

    # -------- Per-cluster sampling + scoring ---------
    rows = []
    for a3m in _cluster_a3ms(pair_id):
        cl_tag = a3m.stem.split("_")[-1]   # "003"
        entries = read_msa(str(a3m))
        if not entries:
            continue
        picked = greedy_select(entries, min(cluster_sample_n, len(entries)))
        if not picked:
            continue

        per_seq = []
        for name, a3m_seq in picked:
            match = _strip_insertions(a3m_seq)
            scored = _score_sequence(match, ddg1, seq1, idx_F1,
                                     ddg2, seq2, idx_F2)
            scored["name"] = name
            per_seq.append(scored)

        if not per_seq:
            continue

        biases = np.asarray([x["bias"] for x in per_seq], dtype=np.float32)
        bias_mean = float(biases.mean())
        bias_std = float(biases.std(ddof=0))
        frac_F1 = float((biases > 0).mean())

        # Per-cluster JSON (includes per-residue contributions for later deep dives)
        (outdir / f"cluster_{cl_tag}.json").write_text(json.dumps({
            "pair_id": pair_id,
            "cluster": cl_tag,
            "tagA": tagA, "tagB": tagB,
            "n_seqs": len(per_seq),
            "bias_mean": bias_mean,
            "bias_std": bias_std,
            "frac_F1_preferring": frac_F1,
            "sequences": per_seq,
        }, indent=2))

        rows.append({
            "pair_id": pair_id,
            "cluster": cl_tag,
            "n_seqs": len(per_seq),
            "bias_mean": bias_mean,
            "bias_std": bias_std,
            "frac_F1_preferring": frac_F1,
            "sum_F1_mean": float(np.mean([x["sum_F1"] for x in per_seq])),
            "sum_F2_mean": float(np.mean([x["sum_F2"] for x in per_seq])),
        })

    # Per-pair summary CSV
    df = pd.DataFrame(rows)
    df_path = analysis_dir / "df_ddg.csv"
    df.to_csv(df_path, index=False)
    print(f"[ddg] wrote {df_path} ({len(df)} clusters) in {time.time()-t_start:.1f}s",
          flush=True)


# ======================================================================
# CLI
# ======================================================================

def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(description="Per-cluster ThermoMPNN conformation-bias scoring.")
    p.add_argument("-input", dest="pair_id", required=True,
                   help="Pair id (e.g., 4dxrA_4dxtA).")
    p.add_argument("--thermompnn_dir", default=os.environ.get("THERMOMPNN_DIR",
                   "/sci/labs/orzuk/orzuk/github/ThermoMPNN"),
                   help="Path to ThermoMPNN cloned repo (default: $THERMOMPNN_DIR).")
    p.add_argument("--device", default="cuda", choices=["cuda", "cpu"])
    p.add_argument("--cluster_sample_n", type=int, default=10)
    p.add_argument("--batch_size", type=int, default=1000,
                   help="ThermoMPNN mutations per forward pass.")
    p.add_argument("--force", action="store_true",
                   help="Recompute ThermoMPNN matrices even if cached.")
    args = p.parse_args(argv)

    run_pair(args.pair_id,
             thermompnn_dir=args.thermompnn_dir,
             device=args.device,
             cluster_sample_n=args.cluster_sample_n,
             batch_size=args.batch_size,
             force=args.force)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except SystemExit:
        raise
    except Exception:
        import traceback
        traceback.print_exc()
        raise SystemExit(1)
