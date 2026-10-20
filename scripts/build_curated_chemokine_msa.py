#!/usr/bin/env python3
"""Curated (Dishman-like) vertebrate chemokine sequence set for the XCL1 pilot.

Replaces the ~8000-sequence MMseqs2/UniRef DeepMsa of pair 2n54B_2hdmA (XCL1;
F1 = 2n54B alternative fold, F2 = 2hdmA chemokine fold) with a curated set of
vertebrate chemokines (cf. Dishman et al., Science 2021: 457 chemokines from 30
vertebrates), so the per-cluster fold-preference + clade-shift analysis can be
re-run on true chemokine clades.

Subcommands
-----------
fetch     (network: UniProt REST + RCSB) species TSV -> curated FASTA + manifest
          + reference sequences (2N54 / 2HDM entity sequences, human XCL1).
align     (offline) curated FASTA -> pipeline DeepMsa.a3m. Row 1 = query
          (2n54B chain), row 2 = partner chain (2hdmA), both headed by their
          chain tags so the pipeline's chain lookup works. Columns = query
          residues (uppercase / '-'), insertions lowercase. Uses MAFFT when
          `mafft` is on PATH, else Biopython pairwise-to-query projection.
build     fetch + align.
annotate  (offline) per-cluster family composition (XCL/CCL/CXCL/CX3CL),
          species groups and identity to human XCL1 for a pilot pair dir, and
          (given a phylo_placement CSV with clade-shift columns) whether each
          method's best clade is the XCL clade or an XCL sub-clade.

FASTA headers: >{accession}|{gene}|{species}|{family}

Examples
--------
  python3 scripts/build_curated_chemokine_msa.py build \\
      --species data/pilots/xcl1_dishman_species.tsv --outdir data/pilots/xcl1_dishman

  python3 scripts/build_curated_chemokine_msa.py align \\
      --fasta data/pilots/xcl1_dishman/chemokines_curated.fasta \\
      --query_fasta Pipeline/FoldPairs/2n54B_2hdmA/fasta_chain_files/2n54B.fasta \\
      --partner_fasta Pipeline/FoldPairs/2n54B_2hdmA/fasta_chain_files/2hdmA.fasta \\
      --out Pipeline/Pilots/xcl1_dishman/2n54B_2hdmA/output_get_msa/DeepMsa.a3m

  python3 scripts/build_curated_chemokine_msa.py annotate \\
      --pair_dir Pipeline/Pilots/xcl1_dishman/2n54B_2hdmA \\
      --manifest data/pilots/xcl1_dishman/manifest.tsv \\
      --phylo_csv Pipeline/Pilots/xcl1_dishman/results/phylo_placement.csv
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from typing import Dict, List, Sequence, Tuple

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
RCSB_FASTA = "https://www.rcsb.org/fasta/entry/{}"
USER_AGENT = "MsaCluster/0.1 (fold-switching pilot; or.zuk@mail.huji.ac.il)"

# UniProt "Protein families" names of the four chemokine subfamilies.
UNIPROT_FAMILIES = {
    "intercrine alpha (chemokine CxC) family": "CXCL",
    "intercrine beta (chemokine CC) family": "CCL",
    "intercrine gamma family": "XCL",
    "intercrine delta family": "CX3CL",
}
# Pfam IL8 / chemokine domain: catches unreviewed entries with no family
# annotation; those are then classified by gene name or cysteine motif.
PFAM_CHEMOKINE = "PF00048"
UNIPROT_FIELDS = ("accession,id,gene_primary,organism_name,organism_id,"
                  "protein_families,reviewed,annotation_score,ft_signal,length,sequence")

GENE_PREFIX_FAMILY = [
    (r"^CX3CL", "CX3CL"), (r"^SCYD", "CX3CL"),
    (r"^CXCL", "CXCL"), (r"^SCYB", "CXCL"), (r"^IL8", "CXCL"),
    (r"^PF4", "CXCL"), (r"^PPBP", "CXCL"),
    (r"^XCL", "XCL"), (r"^SCYC", "XCL"), (r"^LTN", "XCL"),
    (r"^CCL", "CCL"), (r"^SCYA", "CCL"),
]
FAMILIES = ("XCL", "CCL", "CXCL", "CX3CL", "unclassified")
HUMAN_XCL1_ACC = "P47992"
MAMMAL_GROUPS = {"placental", "marsupial", "monotreme"}
MANIFEST_COLS = ["seq_id", "accession", "entry_name", "gene", "species", "taxon_id",
                 "group", "family", "family_source", "uniprot_family", "reviewed",
                 "annotation_score", "signal_end", "precursor_len", "mature_len"]


# ---------------------------------------------------------------------------
# small I/O helpers
# ---------------------------------------------------------------------------
def _http_get(url: str, retries: int = 4, timeout: int = 90):
    last = None
    for k in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read().decode("utf-8"), dict(r.headers)
        except Exception as e:          # transient network / 5xx: back off
            last = e
            time.sleep(2 * (k + 1))
    raise RuntimeError(f"GET failed after {retries} tries: {url} ({last})")


def read_fasta(path: str) -> List[Tuple[str, str]]:
    out, h, cur = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line:
                continue
            if line.startswith(">"):
                if h is not None:
                    out.append((h, "".join(cur)))
                h, cur = line[1:].strip(), []
            else:
                cur.append(line.strip())
    if h is not None:
        out.append((h, "".join(cur)))
    return out


def write_fasta(records: Sequence[Tuple[str, str]], path: str) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", newline="\n") as fh:
        for h, s in records:
            fh.write(f">{h}\n{s}\n")


def read_species(path: str) -> List[Dict[str, str]]:
    with open(path) as fh:
        lines = [l for l in fh if l.strip() and not l.lstrip().startswith("#")]
    return [{"taxon_id": r["taxon_id"].strip(), "name": r["name"].strip(),
             "group": (r.get("group") or "").strip()}
            for r in csv.DictReader(lines, delimiter="\t")]


def _safe(s: str) -> str:
    """Header-safe token: no whitespace, no newick-special characters."""
    s = re.sub(r"\s+", "_", (s or "").strip())
    return re.sub(r"[|(),:;\[\]'\"]", "-", s) or "NA"


def _make_aligner(mode: str = "global"):
    """BLOSUM62 aligner with free end gaps (the chemokine domain sits inside a
    longer precursor / stalk, and the PDB query is short)."""
    from Bio import Align
    from Bio.Align import substitution_matrices
    al = Align.PairwiseAligner()
    al.mode = mode
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score = -10.0
    al.extend_gap_score = -0.5
    al.end_open_gap_score = 0.0
    al.end_extend_gap_score = 0.0
    return al


def _identity(a: str, b: str, aligner) -> float:
    if not a or not b:
        return 0.0
    aln = aligner.align(a, b)[0]
    return aln.counts().identities / float(min(len(a), len(b)))


# ---------------------------------------------------------------------------
# fetch
# ---------------------------------------------------------------------------
def uniprot_query_for_taxon(taxon_id: str) -> str:
    fam = " OR ".join(f'family:"{f}"' for f in UNIPROT_FAMILIES)
    return (f"({fam} OR xref:pfam-{PFAM_CHEMOKINE}) "
            f"AND organism_id:{taxon_id} AND fragment:false")


def fetch_uniprot_tsv(query: str) -> List[Dict[str, str]]:
    url = UNIPROT_SEARCH + "?" + urllib.parse.urlencode(
        {"query": query, "fields": UNIPROT_FIELDS, "format": "tsv", "size": 500})
    rows: List[Dict[str, str]] = []
    while url:
        txt, headers = _http_get(url)
        lines = txt.splitlines()
        if lines:
            rows.extend(csv.DictReader(lines, delimiter="\t"))
        link = headers.get("Link") or headers.get("link") or ""
        m = re.search(r'<([^>]+)>;\s*rel="next"', link)
        url = m.group(1) if m else None
    return rows


def fetch_rcsb_entity_seq(pdb_id: str) -> Tuple[str, str]:
    txt, _ = _http_get(RCSB_FASTA.format(pdb_id))
    recs, h, cur = [], None, []
    for line in txt.splitlines():
        if line.startswith(">"):
            if h is not None:
                recs.append((h, "".join(cur)))
            h, cur = line[1:], []
        elif line.strip():
            cur.append(line.strip())
    if h is not None:
        recs.append((h, "".join(cur)))
    if not recs:
        raise RuntimeError(f"no FASTA records for PDB {pdb_id}")
    return recs[0]


def _signal_end(ft_signal: str) -> int:
    """'SIGNAL 1..21; /evidence=...' -> 21 (0 if absent / unparsable)."""
    m = re.search(r"SIGNAL\s+<?(\d+)\.\.>?(\d+)", ft_signal or "")
    return int(m.group(2)) if m else 0


def classify_family(uniprot_family: str, gene: str, mature: str) -> Tuple[str, str]:
    """(family, source). Priority: UniProt family > gene name > cysteine motif."""
    fam_l = (uniprot_family or "").lower()
    for name, fam in UNIPROT_FAMILIES.items():
        if name in fam_l:
            return fam, "uniprot_family"
    g = (gene or "").upper()
    for pat, fam in GENE_PREFIX_FAMILY:
        if re.match(pat, g):
            return fam, "gene_name"
    # N-terminal cysteine motif of the mature chemokine: CC / CXC / CX3C,
    # or XC (a single N-terminal Cys, two Cys in total) for the XCL family.
    head = mature[:45]
    i = head.find("C")
    if i >= 0:
        nxt = mature[i + 1:i + 5]
        if nxt[:1] == "C":
            return "CCL", "motif"
        if nxt[1:2] == "C":
            return "CXCL", "motif"
        if nxt[3:4] == "C":
            return "CX3CL", "motif"
        if mature.count("C") <= 2:
            return "XCL", "motif"
    return "unclassified", "none"


def record_header(d: Dict[str, str]) -> str:
    return "|".join([_safe(d["accession"]), _safe(d["gene"] or "NA"),
                     _safe(d["species"]), d["family"]])


def curate(raw_by_species: Dict[str, List[Dict[str, str]]],
           species: List[Dict[str, str]], min_len: int, max_len: int,
           max_within_species_id: float) -> Tuple[List[Dict[str, str]], Counter]:
    """Trim signal peptides, classify, filter and deduplicate.

    Deduplication: one entry per (species, gene); within a species, entries
    >= max_within_species_id identical (alleles / isoforms / redundant TrEMBL)
    collapse to one; globally identical sequences collapse to one. Reviewed
    entries (and higher UniProt annotation scores) win.
    Returns (kept records, drop-reason counter).
    """
    aligner = _make_aligner("global")
    drops: Counter = Counter()
    kept_all: List[Dict[str, str]] = []
    seen_seq: Dict[str, str] = {}
    sp_meta = {s["taxon_id"]: s for s in species}
    for taxon, rows in raw_by_species.items():
        meta = sp_meta[taxon]
        recs = []
        for r in rows:
            seq = (r.get("Sequence") or "").strip().upper()
            if not seq:
                drops["no_sequence"] += 1
                continue
            if not (min_len <= len(seq) <= max_len):
                drops["precursor_length"] += 1
                continue
            if sum(ch not in "ACDEFGHIKLMNPQRSTVWY" for ch in seq) > 0.05 * len(seq):
                drops["nonstandard_aa"] += 1
                continue
            sig = _signal_end(r.get("Signal peptide", ""))
            mature = seq[sig:] if 0 < sig < len(seq) - 30 else seq
            mature = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "X", mature)
            gene = (r.get("Gene Names (primary)") or "").strip()
            fam, src = classify_family(r.get("Protein families", ""), gene, mature)
            try:
                score = float(r.get("Annotation") or r.get("Annotation score") or 0)
            except ValueError:
                score = 0.0
            recs.append({
                "accession": r["Entry"], "entry_name": r.get("Entry Name", ""),
                "gene": gene, "species": meta["name"], "taxon_id": taxon,
                "group": meta["group"], "family": fam, "family_source": src,
                "uniprot_family": r.get("Protein families", ""),
                "reviewed": ("reviewed" if (r.get("Reviewed", "").lower() == "reviewed")
                             else "unreviewed"),
                "annotation_score": score, "signal_end": sig,
                "precursor_len": len(seq), "mature_len": len(mature), "mature": mature,
            })
        recs.sort(key=lambda d: (d["reviewed"] != "reviewed", -d["annotation_score"],
                                 d["family"] == "unclassified", -d["mature_len"]))
        by_gene: Dict[str, Dict[str, str]] = {}
        stage1 = []
        for d in recs:
            g = d["gene"].upper()
            if g:
                if g in by_gene:
                    drops["same_gene_in_species"] += 1
                    continue
                by_gene[g] = d
            stage1.append(d)
        stage2: List[Dict[str, str]] = []
        for d in stage1:
            dup = False
            if max_within_species_id < 1.0:
                for k in stage2:
                    la, lb = len(d["mature"]), len(k["mature"])
                    if min(la, lb) / max(la, lb) < 0.8:
                        continue
                    if _identity(d["mature"], k["mature"], aligner) >= max_within_species_id:
                        dup = True
                        break
            if dup:
                drops["near_identical_in_species"] += 1
                continue
            stage2.append(d)
        for d in stage2:
            if d["mature"] in seen_seq:
                drops["identical_sequence"] += 1
                continue
            seen_seq[d["mature"]] = d["accession"]
            kept_all.append(d)
    return kept_all, drops


def summarize(kept: List[Dict[str, str]], species: List[Dict[str, str]],
              drops: Counter) -> str:
    L = []
    fam = Counter(d["family"] for d in kept)
    src = Counter(d["family_source"] for d in kept)
    L.append(f"curated sequences: {len(kept)}  (species: {len(species)})")
    L.append("per family: " + ", ".join(f"{f}={fam.get(f, 0)}" for f in FAMILIES))
    L.append("family source: " + ", ".join(f"{k}={v}" for k, v in src.most_common()))
    L.append("reviewed: " + str(sum(d["reviewed"] == "reviewed" for d in kept)))
    L.append("dropped: " + (", ".join(f"{k}={v}" for k, v in drops.most_common()) or "none"))
    L.append("")
    L.append(f"{'species':<30} {'group':<18} {'n':>4} {'XCL':>4} {'CCL':>4} "
             f"{'CXCL':>5} {'CX3CL':>5} {'uncl':>4}  XCL entries")
    by_sp = defaultdict(list)
    for d in kept:
        by_sp[d["species"]].append(d)
    no_xcl = []
    for s in species:
        ds = by_sp.get(s["name"], [])
        c = Counter(d["family"] for d in ds)
        xcl = [f"{d['accession']}:{d['gene'] or 'NA'}({d['family_source']})"
               for d in ds if d["family"] == "XCL"]
        if not xcl:
            no_xcl.append(s["name"])
        L.append(f"{s['name']:<30} {s['group']:<18} {len(ds):>4} {c.get('XCL', 0):>4} "
                 f"{c.get('CCL', 0):>4} {c.get('CXCL', 0):>5} {c.get('CX3CL', 0):>5} "
                 f"{c.get('unclassified', 0):>4}  {' '.join(xcl)}")
    L.append("")
    L.append(f"species WITHOUT any XCL entry ({len(no_xcl)}): {', '.join(no_xcl) or '-'}")
    return "\n".join(L) + "\n"


def cmd_fetch(args) -> Dict[str, str]:
    species = read_species(args.species)
    os.makedirs(args.outdir, exist_ok=True)
    raw_by_species: Dict[str, List[Dict[str, str]]] = {}
    print(f"[fetch] {len(species)} species from {args.species}")
    for s in species:
        rows = fetch_uniprot_tsv(uniprot_query_for_taxon(s["taxon_id"]))
        raw_by_species[s["taxon_id"]] = rows
        print(f"[fetch]   {s['name']:<30} taxon={s['taxon_id']:<7} UniProt hits={len(rows)}")
    with open(os.path.join(args.outdir, "uniprot_raw.json"), "w") as fh:
        json.dump(raw_by_species, fh)

    kept, drops = curate(raw_by_species, species, args.min_len, args.max_len,
                         args.max_within_species_id)
    fasta = [(record_header(d), d["mature"]) for d in kept]
    fasta_path = os.path.join(args.outdir, "chemokines_curated.fasta")
    write_fasta(fasta, fasta_path)
    man_path = os.path.join(args.outdir, "manifest.tsv")
    with open(man_path, "w", newline="\n") as fh:
        w = csv.DictWriter(fh, fieldnames=MANIFEST_COLS, delimiter="\t",
                           extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for d, (h, _) in zip(kept, fasta):
            w.writerow({**d, "seq_id": h})

    refs = {}
    for pdb, tag in ((args.query_pdb, "query"), (args.partner_pdb, "partner")):
        h, s = fetch_rcsb_entity_seq(pdb)
        refs[tag] = (h, s)
        write_fasta([(f"{pdb}_rcsb|{_safe(h)}", s)],
                    os.path.join(args.outdir, f"rcsb_{pdb}.fasta"))
    hx = [d for d in kept if d["accession"] == HUMAN_XCL1_ACC]
    if hx:
        write_fasta([(record_header(hx[0]), hx[0]["mature"])],
                    os.path.join(args.outdir, "human_XCL1_mature.fasta"))
        print(f"[fetch] human XCL1 {HUMAN_XCL1_ACC}: signal 1..{hx[0]['signal_end']} "
              f"trimmed, mature len {hx[0]['mature_len']}")
    else:
        print(f"[fetch] WARN human XCL1 {HUMAN_XCL1_ACC} not in the curated set")
    print(f"[fetch] query {args.query_pdb} entity len={len(refs['query'][1])}; "
          f"partner {args.partner_pdb} entity len={len(refs['partner'][1])}")
    report = summarize(kept, species, drops)
    with open(os.path.join(args.outdir, "fetch_summary.txt"), "w", newline="\n") as fh:
        fh.write(report)
    print(report)
    print(f"[fetch] wrote {fasta_path} ({len(fasta)} seqs) and {man_path}")
    return {"fasta": fasta_path,
            "query": os.path.join(args.outdir, f"rcsb_{args.query_pdb}.fasta"),
            "partner": os.path.join(args.outdir, f"rcsb_{args.partner_pdb}.fasta")}


# ---------------------------------------------------------------------------
# align -> DeepMsa.a3m
# ---------------------------------------------------------------------------
def project_to_query(q_aln: str, t_aln: str, trim_ends: bool = True) -> str:
    """Project an aligned (query, target) pair into one a3m row: query columns
    uppercase or '-', target insertions lowercase. With trim_ends, target
    residues before the first / after the last query-aligned residue are
    dropped (as a local MMseqs2 hit would be)."""
    assert len(q_aln) == len(t_aln)
    cols = list(zip(q_aln, t_aln))
    matched = [i for i, (q, t) in enumerate(cols) if q != "-" and t != "-"]
    lo, hi = (matched[0], matched[-1]) if matched else (0, -1)
    out = []
    for i, (q, t) in enumerate(cols):
        if q != "-":
            out.append(t.upper() if t != "-" else "-")
        elif t != "-":
            if trim_ends and (i < lo or i > hi):
                continue
            out.append(t.lower())
    return "".join(out)


def align_pairwise(query: str, seqs: List[str]) -> List[Tuple[str, str]]:
    al = _make_aligner("global")
    return [tuple(al.align(query, s)[0]) for s in seqs]


def align_mafft(query: str, seqs: List[str], mafft: str, mafft_args: str,
                threads: int) -> List[Tuple[str, str]]:
    with tempfile.TemporaryDirectory() as td:
        inp = os.path.join(td, "in.fa")
        write_fasta([("q", query)] + [(f"s{i}", s) for i, s in enumerate(seqs)], inp)
        cmd = [mafft] + mafft_args.split() + ["--thread", str(threads), inp]
        print("[align] " + " ".join(cmd))
        res = subprocess.run(cmd, check=True, capture_output=True, text=True)
    recs, h, cur = [], None, []
    for line in res.stdout.splitlines():
        if line.startswith(">"):
            if h is not None:
                recs.append((h, "".join(cur)))
            h, cur = line[1:].strip(), []
        else:
            cur.append(line.strip())
    if h is not None:
        recs.append((h, "".join(cur)))
    aln = {h: s.upper() for h, s in recs}
    q = aln["q"]
    return [(q, aln[f"s{i}"]) for i in range(len(seqs))]


def cmd_align(args) -> None:
    q_recs = read_fasta(args.query_fasta)
    p_recs = read_fasta(args.partner_fasta) if args.partner_fasta else []
    if not q_recs:
        sys.exit(f"[align] empty query FASTA {args.query_fasta}")
    query = re.sub(r"[^A-Za-z]", "", q_recs[0][1]).upper()
    partner = re.sub(r"[^A-Za-z]", "", p_recs[0][1]).upper() if p_recs else ""
    recs = read_fasta(args.fasta)
    print(f"[align] query {args.query_tag} len={len(query)}; partner "
          f"{args.partner_tag} len={len(partner)}; {len(recs)} curated sequences")

    aligner = args.aligner
    mafft = shutil.which("mafft")
    if aligner == "auto":
        aligner = "mafft" if mafft else "pairwise"
    if aligner == "mafft" and not mafft:
        sys.exit("[align] --aligner mafft but mafft is not on PATH")
    seqs = [s for _, s in recs]
    if aligner == "mafft":
        pairs = align_mafft(query, seqs, mafft, args.mafft_args, args.threads)
    else:
        print("[align] mafft not on PATH -> Biopython pairwise-to-query projection "
              "(BLOSUM62, gaps -10/-0.5, free end gaps)")
        pairs = align_pairwise(query, seqs)
    rows = [project_to_query(qa, ta, trim_ends=True) for qa, ta in pairs]

    out = [(args.query_tag, query)]
    if partner:
        # The partner row must keep ALL of its residues: the pipeline matches
        # the insertion-stripped row against the chain FASTA. Align pairwise
        # and do not trim; warn if residues land outside the query columns.
        qa, ta = align_pairwise(query, [partner])[0]
        prow = project_to_query(qa, ta, trim_ends=False)
        n_ins = sum(ch.islower() for ch in prow)
        if n_ins:
            print(f"[align] WARN partner {args.partner_tag}: {n_ins} residue(s) are "
                  f"insertions relative to the query; the pipeline's DeepMsa lookup "
                  f"for {args.partner_tag} strips insertions and may not find it.")
        out.append((args.partner_tag, prow))
    n_bad = 0
    for (h, _), r in zip(recs, rows):
        if sum(ch.isupper() for ch in r) < args.min_query_cov * len(query):
            n_bad += 1
            continue
        out.append((h, r))
    if n_bad:
        print(f"[align] dropped {n_bad} sequences covering "
              f"< {args.min_query_cov:.0%} of the query columns")
    write_fasta(out, args.out)
    cov = sorted(sum(ch.isupper() for ch in r) / len(query) for _, r in out[2:])
    med = f"{cov[len(cov) // 2]:.2f}" if cov else "n/a"
    print(f"[align] wrote {args.out}: {len(out)} rows "
          f"(query + partner + {len(out) - 2} homologs); aligner={aligner}; "
          f"median query coverage={med}")


def cmd_build(args) -> None:
    paths = cmd_fetch(args)
    args.fasta = paths["fasta"]
    args.query_fasta = args.query_fasta or paths["query"]
    args.partner_fasta = args.partner_fasta or paths["partner"]
    args.out = args.out or os.path.join(args.outdir, "DeepMsa.a3m")
    cmd_align(args)


# ---------------------------------------------------------------------------
# annotate: cluster composition + clade-shift interpretation
# ---------------------------------------------------------------------------
def _strip_ins(s: str) -> str:
    return "".join(ch for ch in s if ch == "-" or ch.isupper())


def _pid(a: str, b: str) -> float:
    n = m = 0
    for x, y in zip(a, b):
        if x != "-" and y != "-":
            n += 1
            m += (x == y)
    return m / n if n else float("nan")


def _consensus(rows: List[str]) -> str:
    if not rows:
        return ""
    L = len(rows[0])
    out = []
    for j in range(L):
        c = Counter(r[j] for r in rows if len(r) == L)
        out.append(c.most_common(1)[0][0] if c else "-")
    return "".join(out)


def cmd_annotate(args) -> None:
    with open(args.manifest) as fh:
        man = {r["accession"]: r for r in csv.DictReader(fh, delimiter="\t")}
    pair_dir = args.pair_dir
    q_tag, p_tag = os.path.basename(os.path.normpath(pair_dir)).split("_", 1)
    deep = read_fasta(os.path.join(pair_dir, "output_get_msa", "DeepMsa.a3m"))
    ref = None
    for h, s in deep:
        if h.split()[0].split("|")[0] == HUMAN_XCL1_ACC:
            ref = _strip_ins(s)
    if ref is None:
        sys.exit(f"[annotate] human XCL1 ({HUMAN_XCL1_ACC}) row not found in DeepMsa.a3m")

    def meta_of(h: str) -> Dict[str, str]:
        tok = h.split()[0]
        if tok in (q_tag, p_tag):
            return {"family": "XCL", "group": "placental", "species": "Homo_sapiens_PDB"}
        return man.get(tok.split("|")[0],
                       {"family": "unknown", "group": "unknown", "species": "?"})

    n_xcl_total = sum(1 for h, _ in deep if meta_of(h)["family"] == "XCL")
    clusters: Dict[str, List[Tuple[str, str]]] = {}
    pat = os.path.join(pair_dir, args.cluster_subdir, "ShallowMsa_*.a3m")
    for f in sorted(glob.glob(pat)):
        clusters[os.path.splitext(os.path.basename(f))[0]] = [
            (h, _strip_ins(s)) for h, s in read_fasta(f)]
    if not clusters:
        sys.exit(f"[annotate] no cluster MSAs matching {pat}")

    rows = []
    for tag, members in clusters.items():
        fams = Counter(meta_of(h)["family"] for h, _ in members)
        groups = Counter(meta_of(h)["group"] for h, _ in members)
        n = len(members)
        pids = [p for p in (_pid(s, ref) for _, s in members) if p == p]
        cons = _consensus([s for _, s in members])
        n_mam = sum(groups.get(g, 0) for g in MAMMAL_GROUPS)
        rows.append({
            "cluster": tag, "n": n,
            **{f"n_{f}": fams.get(f, 0) for f in FAMILIES},
            "frac_XCL": round(fams.get("XCL", 0) / n, 3) if n else 0.0,
            "n_mammal": n_mam, "n_nonmammal": n - n_mam,
            "dominant_family": fams.most_common(1)[0][0] if fams else "NA",
            "groups": ";".join(f"{g}:{c}" for g, c in groups.most_common()),
            "mean_pid_to_human_XCL1": (round(sum(pids) / len(pids), 3) if pids
                                       else float("nan")),
            "consensus_pid_to_human_XCL1": round(_pid(cons, ref), 3),
            "contains_pdb_query": any(h.split()[0] in (q_tag, p_tag) for h, _ in members),
        })
    out_csv = args.output or os.path.join(pair_dir, "Analysis",
                                          "cluster_family_composition.csv")
    os.makedirs(os.path.dirname(os.path.abspath(out_csv)), exist_ok=True)
    with open(out_csv, "w", newline="\n") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"[annotate] {len(rows)} clusters; {n_xcl_total} XCL rows in the curated "
          f"DeepMsa -> {out_csv}")
    hdr = (f"{'cluster':<16}{'n':>5}{'XCL':>5}{'CCL':>5}{'CXCL':>6}{'CX3CL':>6}"
           f"{'%XCL':>7}{'mam':>5}{'nonmam':>7}{'meanPID':>9}")
    print(hdr)
    for r in sorted(rows, key=lambda r: -r["frac_XCL"]):
        print(f"{r['cluster']:<16}{r['n']:>5}{r['n_XCL']:>5}{r['n_CCL']:>5}"
              f"{r['n_CXCL']:>6}{r['n_CX3CL']:>6}{r['frac_XCL'] * 100:>6.0f}%"
              f"{r['n_mammal']:>5}{r['n_nonmammal']:>7}"
              f"{r['mean_pid_to_human_XCL1']:>9.3f}")
    comp = {r["cluster"]: r for r in rows}
    if args.phylo_csv:
        _report_clade_shift(args.phylo_csv, comp, args.min_frac_xcl)
    else:
        print("[annotate] no --phylo_csv given: skipping the clade-shift interpretation")


def _report_clade_shift(phylo_csv: str, comp: Dict[str, Dict], min_frac: float) -> None:
    """For each method's best clade-shift clade, say whether it is the XCL
    clade (or a mammalian / non-mammalian XCL sub-clade) and report the
    identity to human XCL1 (compare with Dishman's Anc.3 60% / Anc.4 68%)."""
    if not os.path.isfile(phylo_csv):
        print(f"[annotate] phylo CSV not found: {phylo_csv}")
        return
    with open(phylo_csv) as fh:
        rows = list(csv.DictReader(fh))
    if not rows or "clade_best_members" not in rows[0]:
        print(f"[annotate] {phylo_csv} has no clade_best_members column "
              f"(run scripts/phylo_placement.py --tests clade_shift,...)")
        return
    print()
    print(f"{'method':<12}{'clade_p':>9}{'size':>6}{'%XCL':>7}{'mam':>6}{'nonmam':>8}"
          f"{'meanPID':>9}  verdict")
    for r in rows:
        members = [m for m in (r.get("clade_best_members") or "").split(";") if m]
        if not members:
            continue
        sel = [comp[m] for m in members if m in comp]
        if not sel:
            continue
        n = sum(c["n"] for c in sel)
        n_xcl = sum(c["n_XCL"] for c in sel)
        n_mam = sum(c["n_mammal"] for c in sel)
        n_non = sum(c["n_nonmammal"] for c in sel)
        pid = sum(c["mean_pid_to_human_XCL1"] * c["n"] for c in sel) / n if n else float("nan")
        frac = n_xcl / n if n else 0.0
        xcl_total = sum(c["n_XCL"] for c in comp.values())
        if frac >= min_frac:
            recall = n_xcl / xcl_total if xcl_total else 0.0
            if recall >= 0.8:
                verdict = "XCL clade"
            elif n_mam == 0:
                verdict = "XCL sub-clade (non-mammalian)"
            elif n_non == 0:
                verdict = "XCL sub-clade (mammalian)"
            else:
                verdict = f"XCL sub-clade ({recall:.0%} of XCL seqs)"
        else:
            dom = Counter()
            for c in sel:
                for f in FAMILIES:
                    dom[f] += c[f"n_{f}"]
            verdict = f"NOT XCL (mostly {dom.most_common(1)[0][0]})"
        try:
            p = float(r.get("clade_p", "nan"))
        except ValueError:
            p = float("nan")
        print(f"{r.get('method', '?'):<12}{p:>9.4f}{len(members):>6}{frac * 100:>6.0f}%"
              f"{n_mam:>6}{n_non:>8}{pid:>9.3f}  {verdict}")
    print("  (Dishman et al. 2021: Anc.3 is 60% and Anc.4 68% identical to human XCL1)")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    def add_fetch_args(q):
        q.add_argument("--species", default="data/pilots/xcl1_dishman_species.tsv",
                       help="TSV: taxon_id, name, group (swap in Dishman's list).")
        q.add_argument("--outdir", default="data/pilots/xcl1_dishman")
        q.add_argument("--query_pdb", default="2N54", help="PDB entry of the F1 chain.")
        q.add_argument("--partner_pdb", default="2HDM", help="PDB entry of the F2 chain.")
        q.add_argument("--min_len", type=int, default=60,
                       help="Min precursor length (default 60).")
        q.add_argument("--max_len", type=int, default=450,
                       help="Max precursor length (default 450; CX3CL1/CXCL16 are long).")
        q.add_argument("--max_within_species_id", type=float, default=0.98,
                       help="Collapse within-species sequences at >= this identity.")

    def add_align_args(q, require_inputs: bool):
        q.add_argument("--fasta", required=require_inputs,
                       help="Curated FASTA from the fetch step.")
        q.add_argument("--query_fasta", default=None,
                       help="FASTA of the F1 chain (default: the fetched RCSB entity).")
        q.add_argument("--partner_fasta", default=None,
                       help="FASTA of the F2 chain (row 2 of the a3m).")
        q.add_argument("--query_tag", default="2n54B")
        q.add_argument("--partner_tag", default="2hdmA")
        q.add_argument("--out", default=None, help="Output DeepMsa.a3m path.")
        q.add_argument("--aligner", choices=["auto", "mafft", "pairwise"], default="auto")
        q.add_argument("--mafft_args", default="--maxiterate 1000 --localpair --quiet",
                       help="MAFFT arguments (default: L-INS-i).")
        q.add_argument("--threads", type=int, default=4)
        q.add_argument("--min_query_cov", type=float, default=0.3,
                       help="Drop rows covering less than this fraction of query columns.")

    f = sub.add_parser("fetch", help="UniProt -> curated FASTA + manifest")
    add_fetch_args(f)
    f.set_defaults(func=cmd_fetch)

    a = sub.add_parser("align", help="curated FASTA -> DeepMsa.a3m")
    add_align_args(a, require_inputs=True)
    a.set_defaults(func=cmd_align)

    b = sub.add_parser("build", help="fetch + align")
    add_fetch_args(b)
    add_align_args(b, require_inputs=False)
    b.set_defaults(func=cmd_build)

    n = sub.add_parser("annotate", help="cluster family composition + clade-shift verdict")
    n.add_argument("--pair_dir", required=True,
                   help="Pilot pair dir, e.g. Pipeline/Pilots/xcl1_dishman/2n54B_2hdmA")
    n.add_argument("--manifest", default="data/pilots/xcl1_dishman/manifest.tsv")
    n.add_argument("--cluster_subdir", default="output_msa_cluster")
    n.add_argument("--phylo_csv", default=None,
                   help="phylo_placement CSV with clade_best_members / clade_p.")
    n.add_argument("--min_frac_xcl", type=float, default=0.5,
                   help="A clade counts as XCL when this fraction of its sequences are XCL.")
    n.add_argument("--output", default=None,
                   help="Output CSV (default: <pair_dir>/Analysis/cluster_family_composition.csv)")
    n.set_defaults(func=cmd_annotate)
    return p


def main() -> None:
    args = build_argparser().parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
