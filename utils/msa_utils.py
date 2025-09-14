from Bio import SeqIO, AlignIO, pairwise2
from Bio.Align import substitution_matrices
import parasail  # pip/conda package: parasail
import requests
from io import StringIO
from typing import List, Tuple, Dict
import re, os

_HAS_PARASAIL = True # Assume import parasail will work

# Identity matrix for "standard" (globalxx-like but with positive gap costs)
_ID_ALPHABET = "ARNDCQEGHILKMFPSTWYVBZXUO"
_ID_MATRIX = parasail.matrix_create(_ID_ALPHABET, 1, 0)  # match=1, mismatch=0
_BLOSUM62   = parasail.blosum62
_CIGAR_RE   = re.compile(r"(\d+)([=XMDI])")
_BLOSUM_ALPHABET = set("ARNDCQEGHILKMFPSTWYVBZX")  # no J/U/O in common blosum tables


def _sanitize_blosum(seq: str) -> str:
    # U->C, O->K, J->L ; others outside table -> X
    s = (seq or "").upper().translate(str.maketrans({"U": "C", "O": "K", "J": "L"}))
    return "".join(ch if ch in _BLOSUM_ALPHABET else "X" for ch in s)

_CIGAR_RE = re.compile(r"(\d+)([=XMDISH])")  # include S/H just in case

def _aligned_strings_from_parasail(res, sA: str, sB: str) -> Tuple[str, str, str, Dict[str,int]]:
    """
    Robustly reconstruct (aligned_A, aligned_B) by using traceback if available,
    else parse the CIGAR string. Returns (a, b, cigar_str, op_counts).
    """
    cig = getattr(res, "cigar", None)
    tb  = getattr(res, "traceback", None)

    if tb is not None:
        a, b = tb.query, tb.ref
        return a, b, "<traceback>", {}

    if cig is None:
        raise RuntimeError("parasail result has neither cigar nor traceback")

    cigar_str = str(cig)
    i = getattr(cig, "beg_query", 0)
    j = getattr(cig, "beg_ref", 0)
    a_parts: List[str] = []
    b_parts: List[str] = []
    counts: Dict[str, int] = {"=":0,"X":0,"M":0,"I":0,"D":0,"S":0,"H":0}

    for length, op in _CIGAR_RE.findall(cigar_str):
        n = int(length)
        counts[op] = counts.get(op, 0) + n
        if op in ("=", "X", "M"):
            a_parts.append(sA[i:i+n]); b_parts.append(sB[j:j+n]); i += n; j += n
        elif op == "I" or op == "S":    # insertion in A (gap in B) ; treat 'S' like 'I' for query soft-clip
            a_parts.append(sA[i:i+n]); b_parts.append("-" * n); i += n
        elif op == "D" or op == "H":    # deletion from A (gap in A) ; treat 'H' like 'D' for ref-clip
            a_parts.append("-" * n); b_parts.append(sB[j:j+n]); j += n

    return "".join(a_parts), "".join(b_parts), cigar_str, counts


def _biopython_align(seqA: str, seqB: str, mode: str) -> Tuple[str, str]:
    from Bio import pairwise2
    if mode == "standard":
        aln = pairwise2.align.globalxx(seqA, seqB, one_alignment_only=True)[0]
    else:
        from Bio.Align import substitution_matrices
        blosum = substitution_matrices.load("BLOSUM62")
        aln = pairwise2.align.globalds(seqA, seqB, blosum, -11, -1, one_alignment_only=True)[0]
    return aln.seqA, aln.seqB


def get_align_indexes(
    seqA: str,
    seqB: str,
    mode: str = "blosum",     # "blosum" or "standard"
    gap_open: int = 11,       # POSITIVE for parasail
    gap_extend: int = 1,
    debug: bool = False,
    debug_print_alignment: bool = True,
    debug_max_chars: int = 240,   # max chars to print from aligned strings (head/tail shown)
) -> Tuple[List[int], List[int]]:
    """
    Global alignment (Needleman–Wunsch). Returns lists of indices (i,j) where both seqs have residues (non-gaps).
    Prints detailed debug info if debug=True, including aligned strings (with gaps) and first/last index pairs.
    """
    if mode not in ("blosum", "standard"):
        raise ValueError("mode must be 'standard' or 'blosum'")

    # --- prep sequences
    if mode == "blosum":
        sA = _sanitize_blosum(seqA)
        sB = _sanitize_blosum(seqB)
    else:
        sA = (seqA or "").upper()
        sB = (seqB or "").upper()

    # --- attempt parasail first
    a = b = ""
    cigar_str = ""
    op_counts: Dict[str, int] = {}
    used = "parasail"
    try:
        if not _HAS_PARASAIL:
            raise RuntimeError("parasail not available")

        # choose scoring
        if mode == "blosum":
            matrix = parasail.blosum62
            go, ge = int(abs(gap_open)), int(abs(gap_extend))
        else:
            # identity-like, but with mismatch penalty to avoid all-gap oddities
            matrix = parasail.matrix_create(_ID_ALPHABET, 2, -1)  # match=2, mismatch=-1
            go, ge = 5, 1

        res = parasail.nw_trace_striped_16(sA, sB, go, ge, matrix)
        a, b, cigar_str, op_counts = _aligned_strings_from_parasail(res, sA, sB)

    except Exception as e:
        used = f"biopython (fallback due to: {e})"
        a, b = _biopython_align(sA, sB, mode=mode)
        cigar_str = "<pairwise2>"
        op_counts = {}

    # --- build indices
    i = j = 0
    idxA: List[int] = []
    idxB: List[int] = []
    for ca, cb in zip(a, b):
        if ca != "-" and cb != "-":
            idxA.append(i); idxB.append(j)
        if ca != "-": i += 1
        if cb != "-": j += 1

    # --- debug prints
    if debug:
        aligned_len = len(a)
        both_non_gap = len(idxA)
        matches = sum(1 for x, y in zip(a, b) if x != "-" and y != "-" and x == y)
        mism    = both_non_gap - matches
        gapsA   = a.count("-")
        gapsB   = b.count("-")
        print(f"[align] engine={used}  mode={mode}  aligned_len={aligned_len}  "
              f"both_non_gap={both_non_gap}  matches={matches}  mismatches={mism}  "
              f"gapsA={gapsA}  gapsB={gapsB}")

        if cigar_str:
            show = cigar_str if len(cigar_str) <= 200 else cigar_str[:200] + " ..."
            print(f"[align] CIGAR: {show}")
            if op_counts:
                print(f"[align] op_counts: {op_counts}")

        if debug_print_alignment:
            def _clip(s: str) -> str:
                if len(s) <= debug_max_chars:
                    return s
                half = debug_max_chars // 2
                return s[:half] + " ... " + s[-half:]
            print(f"[align] A_aln: {_clip(a)}")
            print(f"[align] B_aln: {_clip(b)}")

        if both_non_gap:
            head = list(zip(idxA, idxB))[:20]
            tail = list(zip(idxA, idxB))[-10:]
            print(f"[align] first idx pairs (i,j): {head}")
            print(f"[align] last  idx pairs (i,j): {tail}")

    # final guard
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


# --- NEW: write a simple A3M and build a 2-seq seed alignment from a pair ---

def write_a3m(names, aligned_seqs, outfile='seed.a3m'):
    """
    Write a minimal A3M (FASTA-compatible) with gaps as '-'.
    Uppercase = match states (fine for a seed). No inserts needed.
    """
    os.makedirs(os.path.dirname(outfile), exist_ok=True)
    with open(outfile, 'w') as f:
        for nm, seq in zip(names, aligned_seqs):
            f.write(f">{nm}\n{seq}\n")


def build_pair_seed_a3m_from_pair(
    pair_id: str,
    data_dir: str = "Pipeline",
    prefer_chain_files: bool = True,
    out_name: str = "_seed_both.a3m",
    align_mode: str = "blosum"  # "blosum" or "standard"
) -> str:
    """
    Build a 2-sequence aligned seed (A3M) for the given pair using the two chains.

    We look for the sequences in this order:
    1) Pipeline/<pair>/fasta_chain_files/<pdbid><chain>.fasta  (if prefer_chain_files=True)
    2) Pipeline/<pair>/<pdbid>.fasta                           (your existing 'main' FASTA)

    Returns the path to the created A3M file.
    """
    pair_dir = os.path.join(data_dir, pair_id)
    pA, pB = pair_id.split("_", 1)  # e.g. "1fzpD", "2frhA"

    def _pick_fasta(p_with_chain: str) -> Tuple[str, str]:
        """Return (name, seq) from the best available fasta for this PDB+chain."""
        base = p_with_chain[:-1]  # "1fzpD" -> "1fzp"
        candidates = []
        if prefer_chain_files:
            candidates.append(os.path.join(pair_dir, "fasta_chain_files", f"{p_with_chain}.fasta"))
            candidates.append(os.path.join(pair_dir, f"{p_with_chain}.fasta"))
        candidates.append(os.path.join(pair_dir, f"{base}.fasta"))

        for c in candidates:
            if os.path.exists(c):
                ids, seqs = load_fasta(c)
                if len(seqs) >= 1 and len(seqs[0]) > 0:
                    return ids[0], seqs[0]

        raise FileNotFoundError(
            f"Missing FASTA for {p_with_chain}. Tried: {', '.join(candidates)}"
        )

    nameA, seqA = _pick_fasta(pA)
    nameB, seqB = _pick_fasta(pB)

    # Align the two chains once (we keep it simple and robust)
    if align_mode == "blosum":
        blosum = substitution_matrices.load("BLOSUM62")
        aln = pairwise2.align.globalds(seqA, seqB, blosum, -11, -1, one_alignment_only=True)[0]
    else:
        aln = pairwise2.align.globalxx(seqA, seqB, one_alignment_only=True)[0]

    alignedA, alignedB = aln.seqA, aln.seqB
    out_a3m = os.path.join(pair_dir, out_name)
    write_a3m([nameA, nameB], [alignedA, alignedB], out_a3m)
    return out_a3m
