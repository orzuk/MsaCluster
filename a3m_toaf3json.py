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


inp_fa = pathlib.Path(sys.argv[1])
inp_a3m = pathlib.Path(sys.argv[2])
out_json = pathlib.Path(sys.argv[3])

# read target sequence from FASTA
fa = inp_fa.read_text().splitlines()
seq = "".join(l.strip() for l in fa if not l.startswith(">"))

# read MSA; ensure target is first row of the alignment
msa = read_a3m(inp_a3m)
if not msa or msa[0] != seq:
    # put target on top if missing / order differs
    msa = [seq] + [s for s in msa if s != seq]

# convert MSA to A3M string
a3m_text = to_a3m_string(seq, msa)

payload = {
  "name": inp_fa.stem,
  "sequences": [
    {
      "protein": {
        "id": ["A"],
        "sequence": seq,
        "unpairedMsa": a3m_text,   # <- correct per-chain field
        "pairedMsa": "",           # empty for monomer
        "templates": []            # none
      }
    }
  ],
  "modelSeeds": [1],
  "dialect": "alphafold3",
  "version": 1
}


out_json.parent.mkdir(parents=True, exist_ok=True)
out_json.write_text(json.dumps(payload, indent=2))
print("Wrote", out_json)
