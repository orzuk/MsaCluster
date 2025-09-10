#!/usr/bin/env python3
import os
from pathlib import Path
from argparse import ArgumentParser

from colabfold.batch import get_msa_and_templates, msa_to_str

def _ungap(s: str) -> str:
    return (s or "").replace("-", "").replace(".", "").strip()

def _read_fasta_one(path: str) -> tuple[str, str]:
    name, seq = None, []
    with open(path, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                name = line[1:].strip()
            else:
                seq.append(line.strip())
    return (name or Path(path).stem, _ungap("".join(seq)))

def main():
    p = ArgumentParser()
    p.add_argument("input", help="Seed A3M path (e.g., Pipeline/<pair>/_seed_both.a3m)")
    p.add_argument("results", help="Output dir (e.g., Pipeline/<pair>/output_get_msa)")
    p.add_argument("--pair", required=True, help="Pair id, e.g. 1fzpD_2frhA")
    args = p.parse_args()

    seed_a3m   = Path(args.input)
    result_dir = Path(args.results)
    pair_id    = args.pair

    result_dir.mkdir(parents=True, exist_ok=True)
    pair_dir   = seed_a3m.parent

    # figure chain tags from pair id
    foldA, foldB = pair_id.split("_")  # e.g., 1fzpD, 2frhA

    # read the two chain sequences
    faA = pair_dir / "fasta_chain_files" / f"{foldA}.fasta"
    faB = pair_dir / "fasta_chain_files" / f"{foldB}.fasta"
    nameA, seqA = _read_fasta_one(str(faA))
    nameB, seqB = _read_fasta_one(str(faB))

    # read seed a3m lines (some builds prefer having it)
    with open(seed_a3m, "r") as fh:
        a3m_lines = fh.read().splitlines()

    def _make_deep(jobname: str, seq: str):
        try:
            tup = get_msa_and_templates(
                query_sequences=seq,
                jobname=jobname,
                result_dir=result_dir,
                use_templates=False,
                custom_template_path="",
                pair_mode="",
                msa_mode="mmseqs2_uniref_env",
                a3m_lines=a3m_lines,  # works on builds that accept it
            )
        except TypeError:
            tup = get_msa_and_templates(
                query_sequences=seq,
                jobname=jobname,
                result_dir=result_dir,
                use_templates=False,
                custom_template_path="",
                pair_mode="",
                msa_mode="mmseqs2_uniref_env",
            )
        unpaired_msa, paired_msa, quniq, qcard, _ = tup
        msa = msa_to_str(unpaired_msa, paired_msa, quniq, qcard)
        (result_dir / f"{jobname}.a3m").write_text(msa)

    _make_deep(f"DeepMsa_{foldA}", seqA)
    _make_deep(f"DeepMsa_{foldB}", seqB)

if __name__ == "__main__":
    main()
