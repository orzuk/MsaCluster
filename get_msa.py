#!/usr/bin/env python3
import os
from pathlib import Path
from argparse import ArgumentParser
from typing import List, Tuple, Iterable, Dict, Set

# ColabFold (MMseqs2) for MSA building
from colabfold.batch import get_msa_and_templates, msa_to_str


def _ungap_upper(s: str) -> str:
    """Remove gaps/lowercase; keep only uppercase amino acids."""
    return "".join(ch for ch in (s or "") if ch.isupper())


def _clean_a3m_seq(s: str) -> str:
    """Keep uppercase AA and '-' only (strip lowercase insertions)."""
    return "".join(ch for ch in (s or "") if (ch == "-" or ch.isupper()))


def _read_seed_a3m(seed_path: Path) -> Tuple[str, str, str, str]:
    """
    Read a 2-sequence A3M seed. Return:
      (seed_aln_A, seed_aln_B, raw_seq_A, raw_seq_B)
    seed_aln_* are cleaned (uppercase + '-') and have identical length (L_seed).
    """
    headers: List[str] = []
    seqs: List[str] = []
    with open(seed_path, "r") as fh:
        cur_h = None
        cur_s: List[str] = []
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if cur_h is not None:
                    seqs.append("".join(cur_s))
                headers.append(line[1:].strip())
                cur_h = line[1:].strip()
                cur_s = []
            else:
                cur_s.append(line.strip())
        if headers and cur_s:
            seqs.append("".join(cur_s))

    if len(seqs) < 2:
        raise ValueError(f"Seed A3M must contain 2 sequences: {seed_path}")

    # Clean to get the aligned seed rows in the pairwise frame
    seed_aln_A = _clean_a3m_seq(seqs[0])
    seed_aln_B = _clean_a3m_seq(seqs[1])
    if len(seed_aln_A) != len(seed_aln_B):
        raise ValueError("Seed alignment rows must have the same length.")

    raw_seq_A = _ungap_upper(seqs[0])
    raw_seq_B = _ungap_upper(seqs[1])
    return seed_aln_A, seed_aln_B, raw_seq_A, raw_seq_B


def _parse_a3m_text(a3m_text: str) -> List[Tuple[str, str]]:
    """
    Parse an A3M (string) into [(header, seq)] with cleaned sequences (uppercase + '-').
    """
    out: List[Tuple[str, str]] = []
    h = None
    cur: List[str] = []
    for line in a3m_text.splitlines():
        line = line.rstrip("\n")
        if not line:
            continue
        if line.startswith(">"):
            if h is not None:
                out.append((h, _clean_a3m_seq("".join(cur))))
            h = line[1:].strip()
            cur = []
        else:
            cur.append(line.strip())
    if h is not None:
        out.append((h, _clean_a3m_seq("".join(cur))))
    return out


def _expand_to_seed_frame(row_seq: str, seed_aln_for_chain: str) -> str:
    """
    Expand a monomer MSA row (length == len(raw_query)) to the seed alignment frame:
      - For each seed column: if the seed has '-', emit '-';
        else emit next character from row_seq.
    Assumes row_seq length equals number of non-gap positions in seed_aln_for_chain.
    """
    res: List[str] = []
    j = 0
    for c in seed_aln_for_chain:
        if c == "-":
            res.append("-")
        else:
            if j >= len(row_seq):
                # Defensive: pad with '-' if something is off
                res.append("-")
            else:
                res.append(row_seq[j])
                j += 1
    return "".join(res)


def build_seedlocked_deepmsa(seed_a3m: Path, result_dir: Path) -> Path:
    """
    Build a seed-locked Deep MSA:
      1) Read the 2-seq seed alignment (pairwise frame).
      2) Build monomer MSAs for seqA and seqB (ColabFold/MMseqs2).
      3) Project all rows onto the seed frame (no extra columns).
      4) Merge and dedupe, add the two seed rows at the top.
      5) Write DeepMsa.a3m in result_dir.
    """
    result_dir.mkdir(parents=True, exist_ok=True)
    seed_aln_A, seed_aln_B, raw_A, raw_B = _read_seed_a3m(seed_a3m)
    L_seed = len(seed_aln_A)
    L_A = len(raw_A)
    L_B = len(raw_B)

    # --- Build monomer MSA for A ---
    tupA = get_msa_and_templates(
        query_sequences=[raw_A],
        jobname="DeepMsa_S1",
        result_dir=str(result_dir),
        use_templates=False,
        custom_template_path="",
        pair_mode="unpaired",
        msa_mode="mmseqs2_uniref_env",
    )
    unpairedA, pairedA, quniqA, qcardA, _ = tupA
    a3mA = msa_to_str(unpairedA, None, quniqA, qcardA)
    rowsA = _parse_a3m_text(a3mA)  # cleaned, len == L_A

    # --- Build monomer MSA for B ---
    tupB = get_msa_and_templates(
        query_sequences=[raw_B],
        jobname="DeepMsa_S2",
        result_dir=str(result_dir),
        use_templates=False,
        custom_template_path="",
        pair_mode="unpaired",
        msa_mode="mmseqs2_uniref_env",
    )
    unpairedB, pairedB, quniqB, qcardB, _ = tupB
    a3mB = msa_to_str(unpairedB, None, quniqB, qcardB)
    rowsB = _parse_a3m_text(a3mB)  # cleaned, len == L_B

    # --- Project onto the seed frame ---
    expA: List[Tuple[str, str]] = [(f"A|{h}", _expand_to_seed_frame(s, seed_aln_A)) for h, s in rowsA]
    expB: List[Tuple[str, str]] = [(f"B|{h}", _expand_to_seed_frame(s, seed_aln_B)) for h, s in rowsB]

    # --- Merge, dedupe, and add seed rows (first) ---
    out_rows: List[Tuple[str, str]] = [
        (">S1", seed_aln_A),
        (">S2", seed_aln_B),
    ]
    out_rows.extend([(">" + h, s) for (h, s) in expA + expB])

    # Deduplicate by expanded (seed-frame) sequence; keep first occurrence
    seen: Set[str] = set()
    final_rows: List[Tuple[str, str]] = []
    for h, s in out_rows:
        if s in seen:
            continue
        seen.add(s)
        final_rows.append((h, s))

    # --- Write DeepMsa.a3m (seed-locked width) ---
    out_path = result_dir / "DeepMsa.a3m"
    with open(out_path, "w") as fh:
        for h, s in final_rows:
            fh.write(f"{h}\n")
            fh.write(f"{s}\n")

    print(f"[get_msa] Seed-locked DeepMsa.a3m written: {out_path} (width={L_seed}, rows={len(final_rows)})")

    # Also write a small metadata file for downstream mapping if needed
    (result_dir / "DeepMsa.meta.txt").write_text(f"L_seed={L_seed}\nL_A={L_A}\nL_B={L_B}\n")
    return out_path


def main():
    ap = ArgumentParser(description="Build seed-locked Deep MSA (width equals the pairwise seed alignment).")
    ap.add_argument("seed_a3m", type=str, help="Path to 2-sequence pairwise alignment (A3M).")
    ap.add_argument("out_dir", type=str, help="Output directory (e.g., Pipeline/<pair>/output_get_msa).")
    ap.add_argument("--pair", type=str, default="", help="Pair id for logging (optional).")
    args = ap.parse_args()

    seed_a3m = Path(args.seed_a3m).resolve()
    out_dir  = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    # Hint the ColabFold backend (optional; harmless if unset)
    os.environ.setdefault("COLABFOLD_USE_LOCAL_MMSEQS", "0")
    os.environ.setdefault("COLABFOLD_USER_AGENT", "MsaCluster/seedlocked")

    print(f"[get_msa] pair={args.pair} seed={seed_a3m} out_dir={out_dir}")
    build_seedlocked_deepmsa(seed_a3m, out_dir)


if __name__ == "__main__":
    main()
