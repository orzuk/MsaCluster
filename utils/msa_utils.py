from Bio import SeqIO, AlignIO, pairwise2
from Bio.Align import substitution_matrices, PairwiseAligner


try:
    import parasail  # pip/conda package: parasail
    _HAS_PARASAIL = True
except ImportError:
    parasail = None  # type: ignore[assignment]
    _HAS_PARASAIL = False

import requests
from io import StringIO
from typing import List, Tuple, Dict
import re, os
import numpy as np
from collections import Counter
from pathlib import Path
import math

# Identity matrix for "standard" (globalxx-like but with positive gap costs)
_ID_ALPHABET = "ARNDCQEGHILKMFPSTWYVBZXUO"
_ID_MATRIX = parasail.matrix_create(_ID_ALPHABET, 1, 0) if _HAS_PARASAIL else None
_BLOSUM62   = parasail.blosum62 if _HAS_PARASAIL else None
_BLOSUM_ALPHABET = set("ARNDCQEGHILKMFPSTWYVBZX")  # no J/U/O in common blosum tables
_CIGAR_RE = re.compile(r"(\d+)([=XMDISH])")  # include S/H just in case


def _sanitize_blosum(seq: str) -> str:
    # U->C, O->K, J->L ; others outside table -> X
    s = (seq or "").upper().translate(str.maketrans({"U": "C", "O": "K", "J": "L"}))
    return "".join(ch if ch in _BLOSUM_ALPHABET else "X" for ch in s)


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
    if mode == "standard":
        aln = pairwise2.align.globalxx(seqA, seqB, one_alignment_only=True)[0]
    else:
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


# utils/seq_utils.py (or utils/protein_utils.py)

def sequence_identity_percent(seq1: str, seq2: str) -> float:
    """
    Global alignment identity in percent, counting exact matches over aligned (non-gap) positions.
    Returns a float in [0, 100].
    """
    aligner = PairwiseAligner()
    # Needleman–Wunsch style global alignment
    aligner.mode = "global"
    # Simple DNA/protein-ish scoring that measures identity, not similarity
    aligner.match_score = 1.0
    aligner.mismatch_score = 0.0
    # Penalize gaps so we don't collapse; modest affine-like penalties
    aligner.open_gap_score = -1.0
    aligner.extend_gap_score = -0.5

    aln = aligner.align(seq1, seq2)[0]
    a = aln.aligned  # tuple (list of blocks for seq1, list of blocks for seq2)
    # Build aligned strings to count matches cleanly
    s1_aln, s2_aln = aln.format().splitlines()[1], aln.format().splitlines()[3]  # formatted 3-line block
    # If format() changes in your Biopython version, fallback:
    # s1_aln = str(aln).splitlines()[1]
    # s2_aln = str(aln).splitlines()[3]

    matches = 0
    aligned = 0
    for c1, c2 in zip(s1_aln, s2_aln):
        if c1 == '-' or c2 == '-':
            continue
        aligned += 1
        if c1 == c2:
            matches += 1
    if aligned == 0:
        return 0.0
    return 100.0 * matches / aligned


def trim_a3m_for_pair_union(a3m_path, out_path, s1_tokens, s2_tokens, max_len=None):
    """
    Keep columns where S1 OR S2 has a residue; optionally cap to max_len (e.g., 1024).
    Uses a robust downsampler + final hard cap to guarantee length <= max_len.
    """

    lines = Path(a3m_path).read_text().splitlines(True)
    headers, seqs = parse_a3m_records(lines, strip_inserts=True)

    i1 = find_row_idx(headers, s1_tokens, default=0)
    i2 = find_row_idx(headers, s2_tokens, default=1 if len(headers) > 1 else 0)
    s1_aln, s2_aln = seqs[i1], seqs[i2]

    # Columns where at least one of S1/S2 is not a gap
    keep = [j for j, (a, b) in enumerate(zip(s1_aln, s2_aln)) if (a != '-' or b != '-')]

    # Uniformly downsample if too long
    if max_len is not None and len(keep) > max_len:
        # indices via linspace (inclusive of endpoints), floor to ints
        idx = np.floor(np.linspace(0, len(keep) - 1, max_len)).astype(int).tolist()

        # deduplicate while preserving order
        seen = set()
        dedup = []
        for k in idx:
            if k not in seen:
                dedup.append(k); seen.add(k)

        # If dedup became shorter than max_len due to rounding collisions,
        # fill by scanning forward for missing indices.
        ptr = 0
        while len(dedup) < max_len and ptr < len(keep):
            if ptr not in seen:
                dedup.append(ptr); seen.add(ptr)
            ptr += max(1, len(keep) // max_len)

        # Final monotone, bounded indices
        dedup = sorted(set(dedup))[:max_len]
        keep = [keep[k] for k in dedup]

    # Project
    trimmed = project_columns(seqs, keep)

    # Final hard cap (safety against any weirdness): guarantee <= max_len
    if max_len is not None and len(trimmed[0]) > max_len:
        trimmed = [s[:max_len] for s in trimmed]
        keep = keep[:max_len]

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    write_a3m(headers, trimmed, out_path)

    # Save mapping for downstream (optional)
    try:
        Path(out_path).with_suffix(".colmap.txt").write_text(",".join(map(str, keep)))
    except Exception:
        pass

    return keep

def clean_a3m_line(s: str) -> str:
    """
    Make A3M lines alignment-consistent:
    - remove lowercase inserts
    - treat '.' and '*' as gaps
    - keep A–Z, '-', '.','*' semantics
    """
    if not s:
        return ""
    # normalize known gap/stop symbols to '-'
    s = s.replace('.', '-').replace('*', '-')
    # remove all lowercase letters (inserts)
    s = re.sub(r'[a-z]', '', s)
    # keep uppercase letters and '-' only
    s = ''.join(ch for ch in s if (ch == '-') or ('A' <= ch <= 'Z'))
    return s


# --- ADD to msa_utils.py ---

def read_first_two_a3m_rows(a3m_path: str) -> tuple[str, str]:
    """
    Return the first two aligned (gapped) sequences from an A3M/FASTA file,
    cleaned with clean_a3m_line(). Assumes your two query sequences are first.
    """
    rows = []
    with open(a3m_path, "r") as fh:
        name = None
        buf = []
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    rows.append(clean_a3m_line("".join(buf)))
                    if len(rows) == 2:
                        break
                name = line[1:].strip()
                buf = []
            else:
                buf.append(line)
        if name is not None and len(rows) < 2:
            rows.append(clean_a3m_line("".join(buf)))

    if len(rows) < 2:
        raise ValueError(f"{a3m_path}: expected ≥2 sequences, found {len(rows)}")
    if len(rows[0]) != len(rows[1]):
        raise ValueError("First two A3M rows differ in aligned length.")
    return rows[0], rows[1]


def msa_row_to_residx(aln_seq: str) -> np.ndarray:
    """
    Map a single aligned (gapped) A3M row to residue indices:
    returns an array of length = MSA columns, with 0-based residue index at each column,
    or -1 where the row has a gap.
    """
    idx = np.full(len(aln_seq), -1, dtype=int)
    r = 0
    for j, ch in enumerate(aln_seq):
        if ch != '-':
            idx[j] = r
            r += 1
    return idx


# --- add to msa_utils.py ---

def find_two_query_rows_in_a3m(a3m_path: str, keyA: str | None, keyB: str | None) -> tuple[str, str]:
    """
    Return the aligned (gapped) strings of the two *query* sequences from a full A3M/FASTA.
    If keyA/keyB (e.g. '2qqjA','4qdsA') are provided, prefer rows whose headers contain them.
    Fallback: the first two sequences in the file.
    """
    rows = []          # (header, cleaned_seq)
    with open(a3m_path, "r") as fh:
        name, buf = None, []
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    rows.append((name, clean_a3m_line("".join(buf))))
                name, buf = line[1:].strip(), []
            else:
                buf.append(line)
        if name is not None:
            rows.append((name, clean_a3m_line("".join(buf))))

    if len(rows) < 2:
        raise ValueError(f"{a3m_path}: expected ≥2 sequences, found {len(rows)}")

    def _pick(key: str | None) -> int | None:
        if not key:
            return None
        for i, (hdr, _) in enumerate(rows):
            if key in hdr:
                return i
        return None

    iA = _pick(keyA)
    iB = _pick(keyB)
    if iA is not None and iB is not None:
        alnA = rows[iA][1]; alnB = rows[iB][1]
    else:
        # Fallback to first two rows
        alnA = rows[0][1]; alnB = rows[1][1]

    if len(alnA) != len(alnB):
        raise ValueError("Chosen A3M rows differ in aligned length.")
    return alnA, alnB


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


# === A3M parsing/projection helpers ===
def parse_a3m_records(a3m_lines, strip_inserts=True):
    """Return (headers, seqs) lists; optionally strip lowercase inserts."""
    headers, seqs = [], []
    i, n = 0, len(a3m_lines)
    while i < n:
        if a3m_lines[i].startswith(">"):
            h = a3m_lines[i].rstrip("\n")
            i += 1
            s_chunks = []
            while i < n and not a3m_lines[i].startswith(">"):
                s_chunks.append(a3m_lines[i].strip())
                i += 1
            s = "".join(s_chunks)
            if strip_inserts:
                s = "".join(ch for ch in s if (ch == '-' or ch.isupper()))
            headers.append(h)
            seqs.append(s)
        else:
            i += 1
    if not headers:
        raise ValueError("Empty A3M")
    L = len(seqs[0])
    assert all(len(s) == L for s in seqs), "A3M rows must be equal length after stripping."
    return headers, seqs


def project_columns(seqs, keep_idx):
    return ["".join(s[j] for j in keep_idx) for s in seqs]

def keep_mask_union(seq_a, seq_b):
    """Columns where either S1 or S2 has a residue."""
    return [j for j, (a, b) in enumerate(zip(seq_a, seq_b)) if (a != '-' or b != '-')]

def keep_mask_query(seq_q):
    """Columns where the query has a residue (for AF)."""
    return [j for j, ch in enumerate(seq_q) if ch != '-']

def find_row_idx(headers, tokens, default=None):
    """Find a header containing all tokens (case-insensitive)."""
    toks = [t.lower() for t in tokens if t]
    for i, h in enumerate(headers):
        hl = h.lower()
        if all(t in hl for t in toks):
            return i
    return default


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




# Helper function for neff computation
def _pid_row_from_char_array(A: np.ndarray, B_block: np.ndarray) -> np.ndarray:
    """Vectorized pairwise PID(i, block) ignoring gaps. A: (L,), B_block: (B,L)."""
    both = (A != '-') & (B_block != '-')
    denom = both.sum(axis=1)
    matches = (A == B_block) & both
    num = matches.sum(axis=1)
    with np.errstate(divide='ignore', invalid='ignore'):
        pid = np.where(denom > 0, num / denom, 0.0)
    return pid

def compute_neff_from_a3m(
    a3m_path: str,
    id_thresh: float = 0.8,
    use_query_columns: bool = True,
    strip_inserts: bool = True,
    mode: str = "exact",                 # "exact" or "approx"
    approx_n_hashes: int = 3,            # number of random signatures for blocking
    approx_sig_len: int = 32,            # columns per signature
    approx_bucket_cap: int = 5000,       # if a bucket would generate >cap pairs, sample within
    approx_seed: int = 12345,            # RNG seed for reproducibility
    approx_subsample_cap: int = 0,       # optional: if n > cap, subsample sequences first (0=disabled)
) -> int:
    """
    Compute N_effective (Meff) from an A3M file.

    Definition (pairwise, redundancy-aware):
      Neff = sum_i 1 / (|{ j : PID(i,j) >= id_thresh }|)

    - Identity is computed over positions with residues in BOTH sequences (ignoring gaps).
    - If use_query_columns=True, restrict columns to those where the query (first row) has residues.
    - strip_inserts=True removes lowercase insert chars before computing identity.

    Modes:
      * mode="exact":   O(n^2) pairwise (blockwise) – precise.
      * mode="approx":  LSH-ish blocking with several random signatures; do exact PID only within
                        candidate buckets. Non-candidate pairs are treated as below-threshold.
                        This usually **overestimates** Neff slightly (missed neighbors → larger weights).
      You can optionally reduce n first with approx_subsample_cap to bound runtime.

    Returns:
      Integer Meff (rounded) for convenience.
    """
    # --- load & clean ---
    with open(a3m_path, "r") as f:
        lines = f.readlines()

    headers, seqs = parse_a3m_records(lines, strip_inserts=strip_inserts)
    n = len(seqs)
    if n == 0:
        return 0

    # Optional column restriction to query residue positions
    if use_query_columns:
        keep_idx = keep_mask_query(seqs[0])
        seqs = project_columns(seqs, keep_idx)

    # Convert to char array
    arr = np.array([list(s) for s in seqs], dtype='<U1')  # (n, L)
    n, L = arr.shape

    # Optional top-level subsample to bound runtime
    if mode == "approx" and approx_subsample_cap and n > approx_subsample_cap:
        rng = np.random.default_rng(approx_seed)
        pick = rng.choice(n, size=approx_subsample_cap, replace=False)
        arr = arr[pick]
        n = arr.shape[0]

    if mode == "exact":
        # ---- EXACT: O(n^2) blockwise PID counts ----
        counts = np.zeros(n, dtype=np.int32)
        block = 1024
        for i in range(n):
            Ai = arr[i]
            c = 0
            for start in range(0, n, block):
                end = min(n, start + block)
                pid = _pid_row_from_char_array(Ai, arr[start:end])
                c += int((pid >= id_thresh).sum())
            counts[i] = max(c, 1)
        weights = 1.0 / counts.astype(np.float64)
        return int(round(weights.sum()))

    # ---- APPROX: blocking with multiple random signatures ----
    rng = np.random.default_rng(approx_seed)

    # choose candidate columns for signatures – prefer non-gap heavy columns
    # (simple heuristic: any column; feel free to bias toward columns with fewer gaps)
    if L == 0:
        return 0
    col_pool = np.arange(L)

    # For each signature k, pick sig_len columns and make a key for each sequence.
    # A bucket holds indices that match exactly on those columns.
    # Union candidates across multiple signatures reduces false negatives.
    sig_cols: List[np.ndarray] = [
        rng.choice(col_pool, size=min(approx_sig_len, L), replace=False)
        for _ in range(max(1, approx_n_hashes))
    ]

    # Build buckets for each signature
    buckets: List[dict] = []
    for cols in sig_cols:
        # Build keys: join chars at the signature columns; keep '-' as-is.
        # Use bytes to avoid Python string overhead.
        key_arr = arr[:, cols]  # (n, s)
        # Turn into bytes keys
        keys = [bytes("".join(row.tolist()), "ascii", "ignore") for row in key_arr]
        table: dict = {}
        for idx, k in enumerate(keys):
            table.setdefault(k, []).append(idx)
        buckets.append(table)

    # Candidate pairs: we will count neighbors only inside each bucket (with exact PID),
    # and we update counts symmetrically. To avoid O(n^2), we skip massive buckets.
    counts = np.ones(n, dtype=np.int32)  # count self to avoid div-zero

    # helper: exact PID between i and a block of j's
    def pid_row_i_to_js(i: int, js: List[int]) -> np.ndarray:
        A = arr[i]
        B = arr[np.array(js, dtype=int)]
        return _pid_row_from_char_array(A, B)

    # Iterate buckets
    for bidx, table in enumerate(buckets):
        for key, members in table.items():
            m = len(members)
            if m <= 1:
                continue
            # If the bucket is huge, subsample pairs inside to cap runtime
            if approx_bucket_cap > 0 and (m * (m - 1) // 2) > approx_bucket_cap:
                # sample up to approx_bucket_cap pairs uniformly without replacement
                # Generate a small random set of centers and their neighbors
                kcent = max(1, int(np.sqrt(approx_bucket_cap)))
                centers = rng.choice(m, size=min(kcent, m), replace=False)
                for ci in centers:
                    i = members[ci]
                    # sample neighbors for i
                    neigh = rng.choice(m, size=min(kcent, m-1), replace=False)
                    # ensure no self
                    neigh = [members[j] for j in neigh if members[j] != i]
                    if not neigh:
                        continue
                    pid = pid_row_i_to_js(i, neigh)
                    for jj, j in enumerate(neigh):
                        if pid[jj] >= id_thresh:
                            counts[i] += 1
                            counts[j] += 1  # symmetric update
                continue

            # Small/medium bucket: evaluate all pairs inside
            # Do it by rows to keep vectorization benefits
            members_sorted = sorted(members)
            for idx_pos, i in enumerate(members_sorted):
                js = members_sorted[idx_pos + 1 :]
                if not js:
                    continue
                pid = pid_row_i_to_js(i, js)
                hit_mask = (pid >= id_thresh)
                if np.any(hit_mask):
                    hits = [js[k] for k, h in enumerate(hit_mask) if h]
                    counts[i] += len(hits)
                    for j in hits:
                        counts[j] += 1

    # Ensure minimum 1
    counts = np.maximum(counts, 1)
    weights = 1.0 / counts.astype(np.float64)
    meff = weights.sum()
    return int(round(meff))

