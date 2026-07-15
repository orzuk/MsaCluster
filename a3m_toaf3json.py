# a3m_to_af3json.py
import sys, json, pathlib, re


def read_a3m(p):
    import re
    aln = []
    seq = None
    with open(p) as f:
        for raw in f:
            line = raw.rstrip("\r\n")
            if not line:
                continue
            # skip comment / metadata lines commonly present in ColabFold A3M
            if line.startswith(("#", ";")):
                continue
            if line.startswith(">"):
                if seq is not None:
                    aln.append(seq)
                seq = ""
                continue
            # if we haven't seen a header yet, ignore stray lines
            if seq is None:
                continue
            # strip lowercase insertions and '.' gaps typical for A3M
            seq += re.sub(r"[a-z.]", "", line)
    if seq is not None:
        aln.append(seq)
    return aln

# Build a single-string A3M for AF3 (target first). Use the *raw* A3M lines if you prefer,
# but here we reconstruct a minimal A3M with synthetic headers to guarantee target-first.
def to_a3m_string(target_seq, msa_seqs):
    lines = []
    lines.append(">query")
    lines.extend([target_seq])  # target first, uppercase (no lowercase inserts)
    for i, s in enumerate(msa_seqs):
        if s == target_seq:
            continue
        lines.append(f">aln{i+1}")
        lines.extend([s])
    return "\n".join(lines) + "\n"


def _chain_ids(n):
    """Return n AF3 chain IDs: A, B, ... Z, AA, AB, ... (AF3 accepts multi-char)."""
    import string
    out, i = [], 0
    while len(out) < n:
        q, r = divmod(i, 26)
        out.append((string.ascii_uppercase[q - 1] if q else "") + string.ascii_uppercase[r])
        i += 1
    return out


def build_payload(seq, a3m_text, name, ligand_ccds=None, n_copies=1):
    """Build the AF3 input JSON payload.

    apo (default): a single protein chain 'A' with the per-cluster MSA.
    holo:
      - n_copies > 1 : homo-oligomer -- repeat the protein chain (A, B, ...),
        each carrying the same unpaired MSA (oligomerization trigger class).
      - ligand_ccds  : append a ligand entry with the given PDB CCD codes,
        on the next free chain ID (ligand trigger class).
    The protein-partner (protein_binding) case is not built here because the
    partner sequence is not available from the per-cluster inputs; supply it
    upstream if needed.
    """
    n_copies = max(1, int(n_copies))
    ids = _chain_ids(n_copies + (1 if ligand_ccds else 0))
    sequences = []
    for cid in ids[:n_copies]:
        sequences.append({
            "protein": {
                "id": [cid],
                "sequence": seq,
                "unpairedMsa": a3m_text,
                "pairedMsa": "",
                "templates": [],
            }
        })
    if ligand_ccds:
        sequences.append({
            "ligand": {
                "id": [ids[n_copies]],
                "ccdCodes": list(ligand_ccds),
            }
        })
    return {
        "name": name,
        "sequences": sequences,
        "modelSeeds": [1],
        "dialect": "alphafold3",
        "version": 1,
    }


def _main():
    import argparse
    ap = argparse.ArgumentParser(description="A3M (+FASTA) -> AlphaFold3 JSON input")
    ap.add_argument("fasta")
    ap.add_argument("a3m")
    ap.add_argument("out_json")
    ap.add_argument("--ligands", default="",
                    help="Comma/space-separated PDB CCD codes for a holo "
                         "prediction (ligand trigger class), e.g. 'ATP,MG'.")
    ap.add_argument("--copies", type=int, default=1,
                    help="Homo-oligomer protein-chain copies (oligomerization "
                         "trigger class). Default 1 (monomer / apo).")
    args = ap.parse_args()

    inp_fa = pathlib.Path(args.fasta)
    inp_a3m = pathlib.Path(args.a3m)
    out_json = pathlib.Path(args.out_json)

    fa = inp_fa.read_text().splitlines()
    seq = "".join(l.strip() for l in fa if not l.startswith(">"))

    msa = read_a3m(inp_a3m)
    if not msa or msa[0] != seq:
        msa = [seq] + [s for s in msa if s != seq]
    a3m_text = to_a3m_string(seq, msa)

    ligand_ccds = [t for t in args.ligands.replace(",", " ").split() if t]
    payload = build_payload(seq, a3m_text, inp_fa.stem,
                            ligand_ccds=ligand_ccds, n_copies=args.copies)

    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(payload, indent=2))
    print("Wrote", out_json,
          f"(copies={max(1, args.copies)}, ligands={ligand_ccds})")


if __name__ == "__main__":
    _main()
