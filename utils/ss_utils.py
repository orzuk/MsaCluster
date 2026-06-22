"""Single source of truth for secondary-structure (SS) assignment.

EVERY secondary-structure consumer in this project MUST go through ``compute_ss``
so all artifacts agree on one assignment: the Figure-5 per-residue SS track
(``utils/plot_per_residue_ddg``), the 3-D viewers' cartoon (``utils/structure_viewer``,
which authors this SS into the structure so Mol* renders it instead of computing
its own), and any paper main/SI figure. Do NOT re-derive SS anywhere else (neither
by hand - e.g. phi/psi Ramachandran boxes, which over-call helices ~4x - nor
implicitly by letting a viewer/library compute its own).

Engine: pydssp (gold-standard H-bond DSSP, Kabsch-Sander) run on Biopython-parsed
backbone coordinates. Coordinates (not a file path) are handed to pydssp, which
sidesteps mdtraj's failure to parse the project's CIF-as-".pdb" files. Special
cases (e.g. fibril/oligomer inter-chain sheets) belong HERE, in one place.
"""
from __future__ import annotations

import os
from typing import Optional, Tuple

import numpy as np

_AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M",
}


def _load_structure(pdb_path: str):
    """Parse a structure file robustly. The project's ``*_cif.pdb`` files are
    actually mmCIF content despite the extension, so detect by content."""
    from Bio.PDB import PDBParser, MMCIFParser
    head = ""
    try:
        with open(pdb_path, "r", errors="replace") as fh:
            head = fh.read(3000)
    except Exception:
        pass
    is_cif = ("_atom_site." in head) or head.lstrip().startswith("data_")
    parser = MMCIFParser(QUIET=True) if is_cif else PDBParser(QUIET=True)
    return parser.get_structure("x", pdb_path)


def _chain_backbone(ch):
    """(coords (L,4,3) float32, seq, residue_list) for standard-aa, full-backbone
    residues of a Bio.PDB chain, in order."""
    from Bio.PDB.Polypeptide import is_aa
    coords, seq, residues = [], [], []
    for res in ch:
        if not is_aa(res, standard=True):
            continue
        try:
            n = res["N"].coord; ca = res["CA"].coord; c = res["C"].coord
            o = res["O"].coord if "O" in res else res["OXT"].coord
        except KeyError:
            continue
        coords.append([n, ca, c, o])
        seq.append(_AA3TO1.get(res.resname.strip().upper(), "X"))
        residues.append(res)
    return np.asarray(coords, dtype=np.float32), "".join(seq), residues


def _assembly_ss_for_chain(model, target_chain_id) -> str:
    """SS for one chain computed IN THE CONTEXT OF THE WHOLE ASSEMBLY.

    Fibril / oligomer protomers form their β-sheets by hydrogen-bonding to
    *neighbouring* chains (cross-β). DSSP on a single isolated chain finds no
    partner strands, so it calls those residues coil. Running DSSP on the
    concatenated backbone of ALL chains recovers the inter-chain sheet; we then
    return only ``target_chain_id``'s residues. (A few residues at chain joins may
    be mislabelled by the artificial adjacency, but the cross-β core is recovered.)
    Returns '' on failure.
    """
    try:
        import pydssp
        all_coords, mask = [], []
        for ch in model:
            coords, _seq, _res = _chain_backbone(ch)
            same = (ch.id == target_chain_id)
            for row in coords:
                all_coords.append(row)
                mask.append(same)
        if len(all_coords) < 3 or not any(mask):
            return ""
        arr = np.asarray(all_coords, dtype=np.float32)
        raw = pydssp.assign(arr, out_type="c3")
        ss_all = ["H" if str(x) == "H" else "E" if str(x) == "E" else "C" for x in raw]
        return "".join(s for s, m in zip(ss_all, mask) if m)
    except Exception as e:
        print(f"[ss_utils] _assembly_ss_for_chain failed: {e}")
        return ""


def _n_protein_chains(model) -> int:
    from Bio.PDB.Polypeptide import is_aa
    return sum(1 for ch in model
              if any(is_aa(r, standard=True) for r in ch))


def compute_ss_chain(ch) -> str:
    """SS string (H/E/C per standard-aa residue) for an in-memory Bio.PDB chain -
    same engine as compute_ss(). '' on failure. Use this in the 3-D viewer so the
    structure's SS comes from THIS function, not Mol*'s internal DSSP.

    If the chain's parent model holds OTHER protein chains (an oligomer/fibril
    assembly), SS is computed in the full-assembly context so inter-chain β-sheets
    (cross-β) are recovered instead of being mis-called coil - keeping Figure 5 and
    the 3-D cartoon identical for those pairs.
    """
    try:
        import pydssp
        parent = ch.get_parent()
        if parent is not None and _n_protein_chains(parent) > 1:
            ss = _assembly_ss_for_chain(parent, ch.id)
            if ss:
                return ss
        coords, _seq, _res = _chain_backbone(ch)
        if len(coords) < 3:
            return ""
        raw = pydssp.assign(coords, out_type="c3")
        return "".join({"H": "H", "E": "E"}.get(str(x), "C") for x in raw)
    except Exception as e:
        print(f"[ss_utils] compute_ss_chain failed: {e}")
        return ""


def _ss_runs(ss: str, code: str):
    """1-based (beg,end) inclusive runs of `code` in `ss`."""
    out, i, n = [], 0, len(ss)
    while i < n:
        j = i
        while j < n and ss[j] == ss[i]:
            j += 1
        if ss[i] == code:
            out.append((i + 1, j))
        i = j
    return out


def structure_to_mmcif(struct, fold_chain: Optional[str] = None,
                       ss_override: Optional[str] = None) -> str:
    """Serialize a (single-fold-switching-chain, optional ligand) Bio.PDB structure
    to mmCIF text whose secondary structure is OUR compute_ss_chain() result,
    written as _struct_conf (helix) + _struct_sheet_range (strand). Mol* reads
    these natively (verified), so the 3-D cartoon uses our DSSP, not Mol*'s own.
    Protein residues get sequential label_seq_id (1..N) so the SS ranges map
    cleanly; ligand/hetero atoms are emitted as HETATM in their own label_asym_id
    so Mol* still recognises the ligand. Waters are assumed already removed.
    """
    from Bio.PDB.Polypeptide import is_aa
    model = next(iter(struct))
    lines = ["data_x", "#", "loop_",
             "_atom_site.group_PDB", "_atom_site.id", "_atom_site.type_symbol",
             "_atom_site.label_atom_id", "_atom_site.label_comp_id",
             "_atom_site.label_asym_id", "_atom_site.label_seq_id",
             "_atom_site.Cartn_x", "_atom_site.Cartn_y", "_atom_site.Cartn_z",
             "_atom_site.occupancy", "_atom_site.B_iso_or_equiv"]
    aid = 0
    helices, strands = [], []
    het_asym = "L"   # ligand chain label
    for ch in model:
        asym = (str(ch.id) or "A")[:1] or "A"
        # protein residues (sequential seq_id) for this chain
        prot = [r for r in ch if is_aa(r, standard=True)]
        seqid = 0
        for res in prot:
            seqid += 1
            comp = res.resname.strip()[:3]
            for atom in res:
                aid += 1
                x, y, z = atom.coord
                el = (atom.element or atom.get_id()[0]).strip() or "C"
                lines.append(f"ATOM {aid} {el} {atom.get_id()} {comp} {asym} "
                             f"{seqid} {x:.3f} {y:.3f} {z:.3f} 1.00 {atom.get_bfactor():.2f}")
        # SS for this chain's protein residues. ss_override (assembly-context SS
        # for the single displayed chain, precomputed before the viewer reduced
        # the oligomer to one chain) wins when its length matches, so the 3-D
        # cartoon shows the same cross-β as Figure 5.
        ss = compute_ss_chain(ch)
        if ss_override and len(ss_override) == len(prot):
            ss = ss_override
        for b, e in _ss_runs(ss, "H"):
            helices.append((asym, b, e))
        for b, e in _ss_runs(ss, "E"):
            strands.append((asym, b, e))
        # ligand / hetero atoms (no seq_id, HETATM, own asym)
        for res in ch:
            if is_aa(res, standard=True):
                continue
            comp = res.resname.strip()[:3]
            for atom in res:
                aid += 1
                x, y, z = atom.coord
                el = (atom.element or atom.get_id()[0]).strip() or "C"
                lines.append(f"HETATM {aid} {el} {atom.get_id()} {comp} {het_asym} "
                             f". {x:.3f} {y:.3f} {z:.3f} 1.00 {atom.get_bfactor():.2f}")
    if helices:
        lines += ["#", "loop_", "_struct_conf.conf_type_id", "_struct_conf.id",
                  "_struct_conf.beg_label_asym_id", "_struct_conf.beg_label_seq_id",
                  "_struct_conf.end_label_asym_id", "_struct_conf.end_label_seq_id"]
        for k, (a, b, e) in enumerate(helices, 1):
            lines.append(f"HELX_P HELX{k} {a} {b} {a} {e}")
    if strands:
        lines += ["#", "loop_", "_struct_sheet_range.sheet_id",
                  "_struct_sheet_range.id", "_struct_sheet_range.beg_label_asym_id",
                  "_struct_sheet_range.beg_label_seq_id",
                  "_struct_sheet_range.end_label_asym_id",
                  "_struct_sheet_range.end_label_seq_id"]
        for k, (a, b, e) in enumerate(strands, 1):
            lines.append(f"A {k} {a} {b} {a} {e}")
    lines.append("#")
    return "\n".join(lines)


def compute_ss(pdb_path: str, chain: Optional[str] = None) -> Tuple[str, str]:
    """Return ``(ss, seq)`` for one chain - the canonical SS for the whole project.

    ss : str of one letter per residue in {'H','E','C'} (helix/strand/coil),
         assigned by pydssp DSSP. seq : the matching one-letter amino-acid string.
    ``chain`` selects a chain (None -> the longest). Returns ``('', '')`` on any
    failure - callers fall back to drawing no SS. Never raises.
    """
    try:
        import pydssp
        from Bio.PDB.Polypeptide import is_aa
        if not pdb_path or not os.path.isfile(pdb_path):
            return "", ""
        struct = _load_structure(pdb_path)
        model = next(iter(struct))
        if chain and chain in model:
            ch = model[chain]
        else:
            ch = max(model, key=lambda c: sum(1 for _ in c), default=None)
        if ch is None:
            return "", ""
        # Oligomer/fibril: compute SS with all chains present so inter-chain
        # cross-β is recovered (same rule as compute_ss_chain; keeps Fig5 == 3D).
        if _n_protein_chains(model) > 1:
            ss_asm = _assembly_ss_for_chain(model, ch.id)
            if ss_asm:
                _, seq_only, _ = _chain_backbone(ch)
                if len(ss_asm) == len(seq_only):
                    return ss_asm, seq_only
        coords, seq = [], []
        for res in ch:
            if not is_aa(res, standard=True):
                continue
            try:
                n = res["N"].coord
                ca = res["CA"].coord
                c = res["C"].coord
                o = res["O"].coord if "O" in res else res["OXT"].coord
            except KeyError:
                continue  # incomplete backbone -> skip residue
            coords.append([n, ca, c, o])
            seq.append(_AA3TO1.get(res.resname.strip().upper(), "X"))
        if len(coords) < 3:
            return "", ""
        arr = np.asarray(coords, dtype=np.float32)          # (L, 4, 3)
        raw = pydssp.assign(arr, out_type="c3")             # array of '-','H','E'
        ss = "".join({"H": "H", "E": "E"}.get(str(x), "C") for x in raw)
        return ss, "".join(seq)
    except Exception as e:
        print(f"[ss_utils] compute_ss failed for {pdb_path} (chain={chain}): {e}")
        return "", ""
