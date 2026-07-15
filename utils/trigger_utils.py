"""Programmatic access to the fold-switching trigger classification.

The trigger taxonomy is auto-generated from PDB metadata by
``scripts/classify_triggers_from_pdb.py`` and lives at
``docs/triggers_from_pdb.csv``. This module loads that CSV once
(per-process cache) and exposes a small API for any pipeline script
that needs the classification or its underlying evidence (ligand PDB
codes, chain counts, etc.).

Typical usage
-------------
    from utils.trigger_utils import get_trigger_info, is_ligand_triggered

    info = get_trigger_info("1zk9A_3jv6A")
    # -> {"trigger_class": "protein_binding",
    #     "ligands_pdb1": "", "ligands_pdb2": "", ...}

    if is_ligand_triggered("4qhhA_4qhfA"):
        ligand_codes = ligand_codes_for_holo("4qhhA_4qhfA")
        # Use ligand_codes (PDB chemical-component IDs) to compose a
        # holo Boltz-2 YAML.

Canonical column layout of ``docs/triggers_from_pdb.csv`` (column 8 is
``trigger_class``, not column 6 -- earlier diagnostic mistakes have
counted the chain-letter column instead):

    pair_id, pdbid1, chain1_requested, pdbid2, chain2_requested,
    chain1_best_match, chain2_best_match,
    trigger_class,
    evidence, ligands_pdb1, ligands_pdb2,
    n_unique_chains_pdb1, n_unique_chains_pdb2,
    seq_identity, length_diff, notes

Trigger classes (priority order, single-label per pair):
    mutation               pdb1.SEQRES != pdb2.SEQRES (non-trivial diff)
    protein_binding        hetero-complex in one entry, monomer in the other
    ligand                 different non-buffer HETATM sets between entries
    oligomerization        same sequence + ligands, different chain counts
    equilibrium_or_unknown all else
"""

from __future__ import annotations

import csv
import os
from functools import lru_cache
from typing import Dict, List, Optional, Set

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(_THIS_DIR)
DEFAULT_TRIGGERS_CSV = os.path.join(_REPO_ROOT, "docs", "triggers_from_pdb.csv")

# Classes that admit a meaningful holo Boltz-2 experiment.
# (Implementation status: ligand is wired; protein_binding / oligomerization
#  are scaffolded but deferred to follow-up commits.)
TRIGGERED_CLASSES: Set[str] = {
    "ligand", "oligomerization", "protein_binding", "mutation",
}
EQUILIBRIUM_CLASSES: Set[str] = {"equilibrium_or_unknown"}


@lru_cache(maxsize=1)
def _load_triggers_table(path: str = DEFAULT_TRIGGERS_CSV) -> Dict[str, Dict[str, str]]:
    """Load the trigger CSV once and return {pair_id -> row-as-dict}.

    Cached at the module level so repeated calls within one process
    don't re-read the CSV. Pass a non-default path to bypass the cache.
    """
    if not os.path.isfile(path):
        return {}
    out: Dict[str, Dict[str, str]] = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for r in reader:
            pid = r.get("pair_id", "").strip()
            if pid:
                out[pid] = {k: (v or "").strip() for k, v in r.items()}
    return out


def get_trigger_info(pair_id: str,
                     csv_path: Optional[str] = None) -> Optional[Dict[str, str]]:
    """Return the full trigger row for ``pair_id`` or None if absent."""
    table = _load_triggers_table(csv_path or DEFAULT_TRIGGERS_CSV)
    return table.get(pair_id)


def get_trigger_class(pair_id: str,
                      csv_path: Optional[str] = None) -> Optional[str]:
    """Return the trigger_class string for ``pair_id`` or None if absent."""
    info = get_trigger_info(pair_id, csv_path)
    return info["trigger_class"] if info and "trigger_class" in info else None


def pairs_in_class(trigger_class: str,
                   csv_path: Optional[str] = None) -> List[str]:
    """Return all pair_ids assigned to ``trigger_class`` (e.g. "ligand").

    For 'triggered' (the Chakravarty pool) pass 'triggered' as a virtual
    class -- expands to mutation+protein_binding+ligand+oligomerization.
    """
    table = _load_triggers_table(csv_path or DEFAULT_TRIGGERS_CSV)
    if trigger_class == "triggered":
        wanted = TRIGGERED_CLASSES
    elif trigger_class == "equilibrium":
        wanted = EQUILIBRIUM_CLASSES
    else:
        wanted = {trigger_class}
    return [pid for pid, r in table.items() if r.get("trigger_class") in wanted]


def is_ligand_triggered(pair_id: str,
                        csv_path: Optional[str] = None) -> bool:
    return get_trigger_class(pair_id, csv_path) == "ligand"


def is_protein_binding_triggered(pair_id: str,
                                 csv_path: Optional[str] = None) -> bool:
    return get_trigger_class(pair_id, csv_path) == "protein_binding"


def is_oligomerization_triggered(pair_id: str,
                                 csv_path: Optional[str] = None) -> bool:
    return get_trigger_class(pair_id, csv_path) == "oligomerization"


def is_mutation_triggered(pair_id: str,
                          csv_path: Optional[str] = None) -> bool:
    return get_trigger_class(pair_id, csv_path) == "mutation"


def ligand_codes_for_holo(pair_id: str,
                          csv_path: Optional[str] = None) -> List[str]:
    """Return the union of PDB chemical-component IDs across both PDB
    entries for ``pair_id``, suitable for passing as Boltz YAML ``ccd:``
    fields in a holo prediction.

    The trigger CSV stores ligand codes as space- or pipe-separated
    strings under ``ligands_pdb1`` / ``ligands_pdb2``. We strip empties
    and de-duplicate. Returns an empty list if neither column has codes
    (which is the typical situation for non-ligand-triggered pairs).
    """
    info = get_trigger_info(pair_id, csv_path)
    if not info:
        return []
    out: List[str] = []
    for col in ("ligands_pdb1", "ligands_pdb2"):
        raw = info.get(col, "")
        if not raw:
            continue
        # Codes can be space- or pipe-separated; strip + split robustly
        for token in raw.replace("|", " ").split():
            t = token.strip().strip(",")
            if t and t not in out:
                out.append(t)
    return out


def _to_int(x: str, default: int = 1) -> int:
    try:
        return int(str(x).strip())
    except (TypeError, ValueError):
        return default


def oligomer_counts(pair_id: str,
                    csv_path: Optional[str] = None) -> tuple:
    """Return (n_unique_chains_pdb1, n_unique_chains_pdb2) as ints.

    NOTE: these are the chain counts of the entries as recorded by
    classify_triggers_from_pdb.py, which reflect the *extracted* chains,
    not necessarily the full biological assembly. For the 'oligomerization'
    class they are frequently (1, 1) because a single chain was extracted
    from each entry, so they do NOT give the true homo-oligomer stoichiometry
    (see bound_fold / condition_spec below).
    """
    info = get_trigger_info(pair_id, csv_path)
    if not info:
        return (1, 1)
    return (_to_int(info.get("n_unique_chains_pdb1", "1")),
            _to_int(info.get("n_unique_chains_pdb2", "1")))


def bound_fold(pair_id: str,
               csv_path: Optional[str] = None) -> Optional[int]:
    """Return which fold (1 or 2) is the *bound / holo / complexed* state,
    or None when it cannot be determined from the trigger CSV.

    This is the condition->fold mapping needed to score fold preference
    under two conditions (the Delta_switch statistic; see
    scripts/fold_diversity_survey.py). The bound fold is the deposited
    structure that carries the trigger:

      ligand           -> the entry whose ligands_pdb* column is populated
                          (if both, the one with more codes; tie -> None)
      protein_binding  -> the entry with more unique chains (the hetero-complex)
      oligomerization  -> the more oligomeric entry, IF the chain counts differ;
                          the CSV often stores (1, 1) here, in which case the
                          assembly stoichiometry is unknown and we return None
      mutation /        -> None (no apo/holo axis; the two states differ by
      equilibrium          sequence or are both spontaneously populated)
    """
    info = get_trigger_info(pair_id, csv_path)
    if not info:
        return None
    cls = info.get("trigger_class", "")
    if cls == "ligand":
        n1 = len(ligand_codes_for_holo(pair_id, csv_path)) if info.get("ligands_pdb1") else 0
        n2 = 0
        # Count per-entry so we can pick the more-liganded fold.
        c1 = [t for t in info.get("ligands_pdb1", "").replace("|", " ").split() if t]
        c2 = [t for t in info.get("ligands_pdb2", "").replace("|", " ").split() if t]
        if c1 and not c2:
            return 1
        if c2 and not c1:
            return 2
        if c1 and c2:
            if len(c1) > len(c2):
                return 1
            if len(c2) > len(c1):
                return 2
        return None
    if cls in ("protein_binding", "oligomerization"):
        n1, n2 = oligomer_counts(pair_id, csv_path)
        if n1 > n2:
            return 1
        if n2 > n1:
            return 2
        return None
    return None


def condition_spec(pair_id: str,
                   csv_path: Optional[str] = None) -> Dict:
    """Single source of truth for how to run a *holo / condition-aware*
    prediction for ``pair_id`` and how to score it.

    Returns a dict:
        {
          "trigger_class": str,
          "representable": bool,   # can any predictor encode this trigger?
          "bound_fold":    1|2|None,
          "ligand_ccds":   [str],  # PDB CCD codes for the ligand block (ligand class)
          "n_copies":      int|None,  # homo-oligomer copies for the bound fold, if known
          "needs_partner": bool,   # protein_binding: a partner SEQUENCE is required
                                   #   (not stored in the CSV -> must be supplied)
          "note":          str,
        }

    Representability by class (what a structure predictor can actually encode):
      ligand           -> ligand_ccds (fully driven by the CSV)
      protein_binding  -> needs a partner sequence (NOT in the CSV; supply
                          a partner block explicitly, e.g. Boltz --partner_yaml)
      oligomerization  -> homo-oligomer copies; the CSV lacks the true
                          stoichiometry, so n_copies is None unless the chain
                          counts differ, in which case max(n1,n2) is used
      mutation /        -> not representable as an apo/holo axis
      equilibrium
    """
    info = get_trigger_info(pair_id, csv_path)
    cls = (info or {}).get("trigger_class", "") if info else ""
    spec = {
        "trigger_class": cls,
        "representable": False,
        "bound_fold": bound_fold(pair_id, csv_path),
        "ligand_ccds": [],
        "n_copies": None,
        "needs_partner": False,
        "note": "",
    }
    if cls == "ligand":
        spec["ligand_ccds"] = ligand_codes_for_holo(pair_id, csv_path)
        spec["representable"] = bool(spec["ligand_ccds"])
        if not spec["ligand_ccds"]:
            spec["note"] = "ligand class but no CCD codes in CSV"
    elif cls == "protein_binding":
        spec["representable"] = True
        spec["needs_partner"] = True
        spec["note"] = "partner sequence not in CSV; supply a partner block"
    elif cls == "oligomerization":
        spec["representable"] = True
        n1, n2 = oligomer_counts(pair_id, csv_path)
        spec["n_copies"] = max(n1, n2) if max(n1, n2) > 1 else None
        if spec["n_copies"] is None:
            spec["note"] = ("oligomer stoichiometry unknown from CSV "
                            "(chain counts equal); supply --copies")
    elif cls == "mutation":
        spec["note"] = "mutation: the two states differ by sequence, no apo/holo axis"
    else:
        spec["note"] = "equilibrium/unknown: no external trigger"
    return spec


def summarize() -> Dict[str, int]:
    """Return {trigger_class -> count} for the loaded table.

    Useful as a one-line sanity check: ``from utils.trigger_utils
    import summarize; print(summarize())``.
    """
    table = _load_triggers_table()
    counts: Dict[str, int] = {}
    for r in table.values():
        c = r.get("trigger_class", "")
        counts[c] = counts.get(c, 0) + 1
    return counts


if __name__ == "__main__":
    # Smoke test: print the class distribution.
    for cls, n in sorted(summarize().items()):
        print(f"  {cls:<25s} n={n}")
