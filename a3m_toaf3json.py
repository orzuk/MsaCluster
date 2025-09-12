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

# AF3 expects per-chain unpaired MSA(s). We provide one MSA block for chain A.
# deletionMatrix: put zeros (length = sequence length per row) if you don’t have gaps per column handy.
delmat = [[0]*len(s) for s in msa]

payload = {
  "name": inp_fa.stem,
  "sequences": [{"protein": {"id": ["A"], "sequence": seq}}],
  "unpairedMSAs": [
    [  # list per chain; here only chain A
      {
        "descriptions": ["colabfold_a3m"]*len(msa),
        "sequences": msa,
        "deletionMatrix": delmat
      }
    ]
  ],
  "modelSeeds": [1],
  "dialect": "alphafold3",
  "version": 1
}

out_json.parent.mkdir(parents=True, exist_ok=True)
out_json.write_text(json.dumps(payload, indent=2))
print("Wrote", out_json)
