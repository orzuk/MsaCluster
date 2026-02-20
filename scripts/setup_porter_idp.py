#!/usr/bin/env python3
"""
Extract Porter 2022 IDP control proteins from the supplementary Excel file
and create per-protein directories under Pipeline/Controls/PorterIDP/.

Source: Porter & Looger 2022, Protein Science (pro4353-sup-0003-tables2.xlsx)
Contains 99 intrinsically disordered proteins from DisProt.

Usage:
    python scripts/setup_porter_idp.py
"""
import os
import sys
import pandas as pd

# Allow running from repo root
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, ROOT)
from config import PORTER_IDP_DIR

EXCEL_PATH = "/mnt/g/My Drive/Research/ProteinStructure/FoldSwitchEvolution/external_data/pro4353-sup-0003-tables2.xlsx"
CSV_OUTPUT = os.path.join(ROOT, "data", "porter_idp_proteins.csv")


def main():
    # Try Excel first (local machine with Google Drive), fall back to CSV
    if os.path.isfile(EXCEL_PATH):
        df = pd.read_excel(EXCEL_PATH)
        print(f"Read {len(df)} rows from {os.path.basename(EXCEL_PATH)}")
        # Save CSV for reference
        os.makedirs(os.path.dirname(CSV_OUTPUT), exist_ok=True)
        df.to_csv(CSV_OUTPUT, index=False)
        print(f"Saved {CSV_OUTPUT}")
    elif os.path.isfile(CSV_OUTPUT):
        df = pd.read_csv(CSV_OUTPUT)
        print(f"Read {len(df)} rows from {os.path.basename(CSV_OUTPUT)} (Excel not available)")
    else:
        sys.exit(f"ERROR: Neither Excel ({EXCEL_PATH}) nor CSV ({CSV_OUTPUT}) found")

    print(f"Columns: {list(df.columns)}")

    # Expected columns: DisprotID, UniprotID, DisorderedRegion, UniprotSequence
    required = {"UniprotID", "UniprotSequence"}
    if not required.issubset(df.columns):
        sys.exit(f"ERROR: Missing columns. Expected {required}, got {set(df.columns)}")

    # Drop rows with missing sequence
    df = df.dropna(subset=["UniprotID", "UniprotSequence"])
    print(f"After dropping NaNs: {len(df)} proteins")

    # Create per-protein directories and FASTA files
    os.makedirs(PORTER_IDP_DIR, exist_ok=True)
    created = 0
    for _, row in df.iterrows():
        uniprot_id = str(row["UniprotID"]).strip()
        sequence = str(row["UniprotSequence"]).strip()
        disprot_id = str(row.get("DisprotID", "")).strip()

        pdir = os.path.join(PORTER_IDP_DIR, uniprot_id)
        os.makedirs(pdir, exist_ok=True)

        # Write FASTA
        fasta_path = os.path.join(pdir, f"{uniprot_id}.fasta")
        header = f">{uniprot_id}"
        if disprot_id:
            header += f" DisProt:{disprot_id}"
        with open(fasta_path, "w") as fh:
            fh.write(f"{header}\n{sequence}\n")
        created += 1

    print(f"Created {created} protein directories under {PORTER_IDP_DIR}")


if __name__ == "__main__":
    main()
