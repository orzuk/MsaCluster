#!/usr/bin/env python3
"""Run S4PRED per-cluster on sampled representative sequences.

Strategy
--------
For each pair x cluster, sample N (default 10) representative sequences from
the cluster's a3m alignment, write them to a single FASTA, invoke S4PRED in
--save-files mode, parse the per-sequence ss2 outputs, and aggregate to a
single per-cluster Q3 SS prediction by taking the per-residue MAJORITY VOTE
across the 10 representatives. Also writes per-cluster %H, %E, %C and median
H/E/C confidences.

Outputs (per pair, written under Pipeline/FoldPairs/<pair>/Analysis/):
  s4pred_<cluster>.ss2          : aggregated SS prediction in PSIPRED VFORMAT-like
                                  format, plus per-residue %H / %E / %C across
                                  the sampled reps.
  df_s4pred.csv                 : one row per cluster (incl. DeepMsa) with
                                  pair_id, cluster, n_reps, n_residues,
                                  frac_H, frac_E, frac_C, ss_q3

Per-cluster ss_q3 is the H/E/C string after majority vote.
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
import random
import subprocess
import sys
import tempfile
from collections import Counter

GAP_CHARS = set("-.")


def read_a3m(path: str) -> list[tuple[str, str]]:
    """Return list of (header, aligned_seq) from an a3m (insert columns stripped)."""
    out: list[tuple[str, str]] = []
    cur_hdr = None
    cur: list[str] = []
    try:
        with open(path) as f:
            for line in f:
                line = line.rstrip()
                if not line:
                    continue
                if line.startswith(">"):
                    if cur_hdr is not None:
                        out.append((cur_hdr, "".join(cur)))
                    cur_hdr = line[1:].split()[0]
                    cur = []
                else:
                    cur.append("".join(c for c in line if not c.islower()))
            if cur_hdr is not None:
                out.append((cur_hdr, "".join(cur)))
    except OSError:
        return []
    return out


def sample_representatives(seqs: list[tuple[str, str]], n: int,
                            rng: random.Random) -> list[tuple[str, str]]:
    """Sample n diverse representatives. Strategy: take the FIRST sequence (the
    query / cluster reference) plus n-1 random others (without replacement).
    For small clusters return all seqs."""
    if len(seqs) <= n:
        return seqs
    head = [seqs[0]]
    tail = rng.sample(seqs[1:], k=n - 1)
    return head + tail


def write_fasta(reps: list[tuple[str, str]], path: str) -> None:
    """Write a FASTA with ungapped sequences (S4PRED is single-sequence)."""
    with open(path, "w") as f:
        for i, (hdr, seq) in enumerate(reps):
            ungapped = "".join(c for c in seq if c not in GAP_CHARS).upper()
            # Replace any non-canonical letter with X for safety
            ungapped = "".join(c if c.isalpha() and c in "ACDEFGHIKLMNPQRSTVWY" else "X"
                                for c in ungapped)
            if not ungapped:
                continue
            f.write(f">rep_{i}\n{ungapped}\n")


def run_s4pred(fasta: str, s4pred_dir: str, out_dir: str, device: str) -> bool:
    """Invoke S4PRED with --save-files. Returns True on success."""
    os.makedirs(out_dir, exist_ok=True)
    script = os.path.join(s4pred_dir, "run_model.py")
    cmd = ["python", script,
           "--device", device,
           "--save-files",
           "--outdir", out_dir,
           "--save-by-idx",
           "--silent",
           fasta]
    try:
        subprocess.run(cmd, capture_output=True, text=True, check=True,
                       cwd=s4pred_dir, timeout=600)
        return True
    except subprocess.CalledProcessError as e:
        print(f"[warn] S4PRED failed on {fasta}: {e.stderr[-300:]}",
              file=sys.stderr)
        return False


def parse_ss2(path: str) -> list[tuple[str, str, float, float, float]]:
    """Parse a PSIPRED VFORMAT .ss2 file.
    Returns list of (aa, ss_class, p_C, p_H, p_E) per residue, in order.
    """
    out = []
    try:
        with open(path) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.split()
                if len(parts) < 6:
                    continue
                try:
                    idx = int(parts[0])
                    aa = parts[1]
                    ss = parts[2]
                    p_c = float(parts[3])
                    p_h = float(parts[4])
                    p_e = float(parts[5])
                    out.append((aa, ss, p_c, p_h, p_e))
                except (ValueError, IndexError):
                    continue
    except OSError:
        return []
    return out


def majority_q3(positions: list[list[str]]) -> str:
    """positions[i] = list of Q3 calls (H/E/C) across reps for residue i.
    Return majority-vote string."""
    out = []
    for pos_calls in positions:
        if not pos_calls:
            out.append("-")
            continue
        c = Counter(pos_calls).most_common(1)[0][0]
        out.append(c)
    return "".join(out)


def process_cluster(pair_id: str, cluster_tag: str, a3m_path: str,
                     s4pred_dir: str, analysis_dir: str, work_dir: str,
                     n_reps: int, device: str, rng: random.Random) -> dict | None:
    seqs = read_a3m(a3m_path)
    if not seqs:
        return None
    reps = sample_representatives(seqs, n_reps, rng)
    fasta_path = os.path.join(work_dir, f"{pair_id}__{cluster_tag}.fas")
    pred_dir = os.path.join(work_dir, f"{pair_id}__{cluster_tag}_preds")
    write_fasta(reps, fasta_path)
    ok = run_s4pred(fasta_path, s4pred_dir, pred_dir, device)
    if not ok:
        return None
    # Collect per-rep predictions
    rep_preds: list[list[tuple]] = []
    for f in sorted(glob.glob(os.path.join(pred_dir, "*.ss2"))):
        parsed = parse_ss2(f)
        if parsed:
            rep_preds.append(parsed)
    if not rep_preds:
        return None
    # All reps may have different sequence lengths if the reference seq frame
    # differs; align by taking the MIN length and truncating.
    L = min(len(p) for p in rep_preds)
    pos_q3 = [[p[i][1] for p in rep_preds] for i in range(L)]
    ss_string = majority_q3(pos_q3)
    aa_string = "".join(rep_preds[0][i][0] for i in range(L))
    n_h = ss_string.count("H"); n_e = ss_string.count("E"); n_c = ss_string.count("C")
    n = n_h + n_e + n_c if (n_h + n_e + n_c) else 1
    # Write aggregated ss2-style file
    os.makedirs(analysis_dir, exist_ok=True)
    out_ss2 = os.path.join(analysis_dir, f"s4pred_{cluster_tag}.ss2")
    with open(out_ss2, "w") as f:
        f.write(f"# S4PRED majority-vote across {len(rep_preds)} reps\n")
        for i in range(L):
            # average confidences across reps
            pCs = [p[i][2] for p in rep_preds]
            pHs = [p[i][3] for p in rep_preds]
            pEs = [p[i][4] for p in rep_preds]
            mean_pC = sum(pCs) / len(pCs)
            mean_pH = sum(pHs) / len(pHs)
            mean_pE = sum(pEs) / len(pEs)
            f.write(f"{i+1:4d} {aa_string[i]} {ss_string[i]}  "
                    f"{mean_pC:.3f}  {mean_pH:.3f}  {mean_pE:.3f}\n")
    # Optional cleanup of per-rep predictions
    for f in glob.glob(os.path.join(pred_dir, "*.ss2")):
        try: os.remove(f)
        except OSError: pass
    try:
        os.rmdir(pred_dir)
    except OSError:
        pass
    try:
        os.remove(fasta_path)
    except OSError:
        pass
    return {
        "pair_id": pair_id,
        "cluster": cluster_tag,
        "n_reps": len(rep_preds),
        "n_residues": L,
        "frac_H": n_h / n,
        "frac_E": n_e / n,
        "frac_C": n_c / n,
        "ss_q3": ss_string,
    }


def process_pair(pair_id: str, fpairs_root: str, s4pred_dir: str,
                  work_root: str, n_reps: int, device: str,
                  rng: random.Random) -> list[dict]:
    rows: list[dict] = []
    shallow_dir = os.path.join(fpairs_root, pair_id, "output_msa_cluster")
    deep_file = os.path.join(fpairs_root, pair_id, "output_get_msa",
                              "DeepMsa.a3m")
    analysis_dir = os.path.join(fpairs_root, pair_id, "Analysis")
    pair_work = os.path.join(work_root, pair_id)
    os.makedirs(pair_work, exist_ok=True)
    cluster_files = sorted(glob.glob(os.path.join(shallow_dir, "ShallowMsa_*.a3m")))
    if os.path.isfile(deep_file):
        cluster_files.append(deep_file)
    for cf in cluster_files:
        ctag = "DeepMsa" if "DeepMsa" in cf else os.path.splitext(os.path.basename(cf))[0]
        r = process_cluster(pair_id, ctag, cf, s4pred_dir,
                             analysis_dir, pair_work, n_reps, device, rng)
        if r:
            rows.append(r)
    # Save per-pair df
    if rows:
        out_csv = os.path.join(analysis_dir, "df_s4pred.csv")
        with open(out_csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader(); w.writerows(rows)
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pairs", default="ALL")
    ap.add_argument("--fpairs_root", default="Pipeline/FoldPairs")
    ap.add_argument("--s4pred_dir", required=True,
                    help="Path to s4pred repo (containing run_model.py).")
    ap.add_argument("--work_root", default="Pipeline/FoldPairs/_s4pred_work")
    ap.add_argument("--n_reps", type=int, default=10)
    ap.add_argument("--device", default="cpu", choices=["cpu", "gpu"])
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", default="docs/s4pred_per_cluster_summary.csv")
    args = ap.parse_args()

    rng = random.Random(args.seed)

    if args.pairs.upper() == "ALL":
        pair_dirs = sorted(d for d in os.listdir(args.fpairs_root)
                            if os.path.isdir(os.path.join(args.fpairs_root, d)))
    else:
        pair_dirs = [s.strip() for s in args.pairs.replace(",", " ").split()
                      if s.strip()]

    all_rows: list[dict] = []
    done = 0
    for pid in pair_dirs:
        rows = process_pair(pid, args.fpairs_root, args.s4pred_dir,
                             args.work_root, args.n_reps, args.device, rng)
        if rows:
            all_rows.extend(rows); done += 1
            print(f"  [{done}] {pid}: {len(rows)} clusters", flush=True)
    if not all_rows:
        print("[error] no rows collected", file=sys.stderr); sys.exit(2)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    fields_full = list(all_rows[0].keys())
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields_full)
        w.writeheader(); w.writerows(all_rows)
    print(f"\n[ok] wrote {len(all_rows)} cluster rows -> {args.out}")


if __name__ == "__main__":
    main()
