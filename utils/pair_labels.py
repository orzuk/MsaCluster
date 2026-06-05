"""Single source of truth for human-readable pair labels.

Everything that needs a protein name / organism for a fold-switch pair -- figures,
HTML tables, the paper -- should read it from here, which loads the ONE canonical
table docs/pair_protein_names.csv (built by scripts/fetch_pair_protein_names.py).

    from utils.pair_labels import pair_display, pair_organism, pair_label
    pair_display("2n54B_2hdmA")  -> "XCL1"
    pair_organism("2n54B_2hdmA") -> "Homo sapiens"
    pair_label("2n54B_2hdmA")    -> "XCL1 (2n54B_2hdmA)"
"""
import os
import functools
import pandas as pd

from config import TABLES_RES

NAMES_CSV = os.path.join(TABLES_RES, "pair_protein_names.csv")


@functools.lru_cache(maxsize=1)
def _table():
    """pair_id -> {display_name, organism, ...} (empty dict if the CSV is absent)."""
    if not (os.path.isfile(NAMES_CSV) and os.path.getsize(NAMES_CSV) > 0):
        return {}
    df = pd.read_csv(NAMES_CSV).fillna("")
    return {str(r["pair_id"]): r.to_dict() for _, r in df.iterrows()}


def pair_display(pair_id: str) -> str:
    """Short human name for a pair (alias/cleaned protein name); '' if unknown."""
    return str(_table().get(str(pair_id), {}).get("display_name", "") or "")


def pair_organism(pair_id: str) -> str:
    """Source organism(s) for a pair; '' if unknown."""
    return str(_table().get(str(pair_id), {}).get("organism", "") or "")


def pair_label(pair_id: str, with_id: bool = True) -> str:
    """Display name plus the pair code, e.g. 'XCL1 (2n54B_2hdmA)'.

    Falls back to just the id when no name is known."""
    name = pair_display(pair_id)
    if not name:
        return str(pair_id)
    return f"{name} ({pair_id})" if with_id else name
