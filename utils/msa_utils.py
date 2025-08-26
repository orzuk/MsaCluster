from Bio import SeqIO, AlignIO, pairwise2
from Bio.Align import substitution_matrices
import parasail  # pip/conda package: parasail
import requests
from io import StringIO
# from __future__ import annotations
from typing import List, Tuple, Dict
import re



# Identity matrix for "standard" (globalxx-like but with positive gap costs)
_ID_ALPHABET = "ARNDCQEGHILKMFPSTWYVBZXUO"
_ID_MATRIX = parasail.matrix_create(_ID_ALPHABET, 1, 0)  # match=1, mismatch=0
_BLOSUM62   = parasail.blosum62
_CIGAR_RE   = re.compile(r"(\d+)([=XMDI])")
_BLOSUM_ALPHABET = set("ARNDCQEGHILKMFPSTWYVBZX")  # no J/U/O in common blosum tables

def _sanitize_blosum(seq: str) -> str:
    # map unsupported letters to reasonable substitutes; other unknowns -> X
    s = (seq or "").upper().translate(str.maketrans({"U":"C","O":"K","J":"L"}))
    return "".join(ch if ch in _BLOSUM_ALPHABET else "X" for ch in s)

def _aligned_strings_from_cigar(res, sA: str, sB: str) -> Tuple[str, str, str, Dict[str,int]]:
    cig = getattr(res, "cigar", None)
    if cig is None:
        tb = getattr(res, "traceback", None)
        if tb is None:
            raise RuntimeError("parasail result has neither cigar nor traceback")
        a, b = tb.query, tb.ref
        return a, b, "<no-cigar>", {}
    cigar_str = str(cig)
    i = getattr(cig, "beg_query", 0)
    j = getattr(cig, "beg_ref", 0)
    a_parts, b_parts = [], []
    counts = {"=":0,"X":0,"M":0,"I":0,"D":0}
    for length, op in _CIGAR_RE.findall(cigar_str):
        n = int(length); counts[op] = counts.get(op,0) + n
        if op in ("=","X","M"):
            a_parts.append(sA[i:i+n]); b_parts.append(sB[j:j+n]); i += n; j += n
        elif op == "I":  # insertion in A (gap in B)
            a_parts.append(sA[i:i+n]); b_parts.append("-"*n); i += n
        elif op == "D":  # deletion from A (gap in A)
            a_parts.append("-"*n); b_parts.append(sB[j:j+n]); j += n
    return "".join(a_parts), "".join(b_parts), cigar_str, counts

def get_align_indexes(
    seqA: str,
    seqB: str,
    mode: str = "blosum",     # "blosum" or "standard"
    gap_open: int = 11,       # POSITIVE costs (parasail requirement)
    gap_extend: int = 1,
    debug: bool = False,
) -> Tuple[List[int], List[int]]:
    """
    Global alignment via parasail; returns index lists (i,j) where both sequences are aligned (non-gaps).
    Robust across Windows builds (parses CIGAR), with helpful debug output.
    """
    if mode not in ("blosum","standard"):
        raise ValueError("mode must be 'standard' or 'blosum'")

    # Prep sequences
    if mode == "blosum":
        sA = _sanitize_blosum(seqA); sB = _sanitize_blosum(seqB)
        mat = _BLOSUM62
        go, ge = int(abs(gap_open)), int(abs(gap_extend))
    else:
        sA = (seqA or "").upper(); sB = (seqB or "").upper()
        mat = _ID_MATRIX
        # IMPORTANT: use POSITIVE gap costs to avoid all-gap alignments
        go, ge = 5, 1  # conservative defaults for "standard"

    # Align (trace-enabled, SIMD)
    res = parasail.nw_trace_striped_16(sA, sB, go, ge, mat)
    a, b, cigar_str, op_counts = _aligned_strings_from_cigar(res, sA, sB)

    # Build index mapping (both non-gaps)
    i = j = 0
    idxA: List[int] = []
    idxB: List[int] = []
    for ca, cb in zip(a, b):
        if ca != "-" and cb != "-":
            idxA.append(i); idxB.append(j)
        if ca != "-": i += 1
        if cb != "-": j += 1

    if debug:
        aligned_len = len(a)
        both_non_gap = len(idxA)
        matches = sum(1 for x,y in zip(a,b) if x != "-" and y != "-" and x == y)
        mism    = both_non_gap - matches
        print(f"[align] mode={mode} | aligned_len={aligned_len} both_non_gap={both_non_gap} "
              f"matches={matches} mismatches={mism} gapsA={a.count('-')} gapsB={b.count('-')}")
        if cigar_str != "<no-cigar>":
            print(f"[align] cigar(head): {cigar_str[:80]} ...")
            print(f"[align] ops: {op_counts}")

    # Safety fallback: if something went wrong (shouldn’t, with positive gaps), try a stricter identity scoring
    if not idxA:
        if debug:
            print("[align][warn] zero aligned columns; retrying with stricter identity scoring")
        # Stricter: identity match=2, mismatch=-1 (forces aligned pairs)
        mat2 = parasail.matrix_create(_ID_ALPHABET, 2, -1)
        res2 = parasail.nw_trace_striped_16(sA, sB, 5, 1, mat2)
        a2, b2, *_ = _aligned_strings_from_cigar(res2, sA, sB)
        i = j = 0; idxA=[]; idxB=[]
        for ca, cb in zip(a2, b2):
            if ca != "-" and cb != "-":
                idxA.append(i); idxB.append(j)
            if ca != "-": i += 1
            if cb != "-": j += 1

    # Final guard
    if not idxA:
        raise ValueError("No matched residues found from sequence alignment.")

    return idxA, idxB


def lprint(string, f):
    print(string)
    f.write(string + '\n')


def load_fasta(fil):
    seqs, IDs = [], []
    with open(fil) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq = ''.join([x for x in record.seq])
            IDs.append(record.id)
            seqs.append(seq)
    return IDs, seqs


def write_fasta(names, seqs, outfile='tmp.fasta'):
    with open(outfile, 'w') as f:
        for nm, seq in list(zip(names, seqs)):
            f.write(">%s\n%s\n" % (nm, seq))


# Do pairwise sequence alignment and keep the matching indices
def get_align_indexes_old(seqA, seqB):
    alignments = pairwise2.align.globalxx(seqA, seqB, one_alignment_only=True)
    best_align = alignments[0]
    seqA = best_align.seqA
    seqB = best_align.seqB
    cursA, cursB = 0, 0
    seqA_idxs, seqB_idxs = [], []
    for aa in range(len(seqA)):
        if (seqA[aa] != '-') & (seqB[aa] != '-'):
            seqA_idxs.append(cursA)
            seqB_idxs.append(cursB)
            cursA += 1
            cursB += 1
        if (seqA[aa] == '-') & (seqB[aa] != '-'):
            cursB += 1
        if (seqA[aa] != '-') & (seqB[aa] == '-'):
            cursA += 1
    return seqA_idxs, seqB_idxs



def download_and_parse_pfam_msa(pfam_id, alignment_type="seed"):
    """
    Downloads and parses a Pfam MSA using the REST API (SEED only).

    Parameters:
        pfam_id (str): Pfam family ID (e.g., 'PF00069')
        alignment_type (str): Only 'seed' is supported currently

    Returns:
        MultipleSeqAlignment: Parsed MSA object from Biopython
    """
    if alignment_type != "seed":
        raise NotImplementedError("Only 'seed' alignment is available via the Pfam REST API")

    url = f"https://www.ebi.ac.uk/interpro/api/pfam/entry/pfam/{pfam_id}/alignment/seed/stockholm/"

    response = requests.get(url)
    if response.status_code != 200:
        raise Exception(f"Download failed with status code {response.status_code}: {url}")

    sto_text = response.text
    msa_io = StringIO(sto_text)
    alignment = AlignIO.read(msa_io, "stockholm")

    return alignment
