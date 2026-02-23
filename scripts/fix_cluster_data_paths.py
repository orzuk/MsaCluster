#!/usr/bin/env python3
"""
Fix stale data paths on the cluster after the Phase 7 directory restructuring.

Two fixes:
  1) samples_index.tsv: rewrite pdb_path entries that contain the old
     "Pipeline/<pair>/" layout to use "Pipeline/FoldPairs/<pair>/".
  2) output_cmap_esm/: if any pair still has the legacy cmap directory,
     move its contents to output_cmaps/msa_transformer/ (what the code expects).

Run from the repo root on the cluster:
    python scripts/fix_cluster_data_paths.py             # dry-run (report only)
    python scripts/fix_cluster_data_paths.py --apply      # actually fix

Safe to run multiple times (idempotent).
"""
import os
import re
import shutil
import sys
import argparse

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import DATA_DIR

PAIR_RE = re.compile(r"^\w{4,6}_\w{4,6}$")


def fix_samples_index_tsv(pair_dir, pair_id, apply=False):
    """Fix stale absolute paths in samples_index.tsv files."""
    fixes = []
    for model in ("esm2", "esm3"):
        tsv = os.path.join(pair_dir, "output_esm_fold", model, "samples_index.tsv")
        if not os.path.isfile(tsv):
            continue
        with open(tsv, "r") as f:
            lines = f.readlines()

        new_lines = []
        changed = False
        for line in lines:
            orig = line
            # Fix pattern: /Pipeline/<pair_id>/ -> /Pipeline/FoldPairs/<pair_id>/
            # This handles any absolute path prefix
            old_pat = f"/Pipeline/{pair_id}/"
            new_pat = f"/Pipeline/FoldPairs/{pair_id}/"
            if old_pat in line and new_pat not in line:
                line = line.replace(old_pat, new_pat)
                changed = True
            new_lines.append(line)

        if changed:
            n_fixed = sum(1 for o, n in zip(lines, new_lines) if o != n)
            fixes.append((tsv, n_fixed))
            if apply:
                with open(tsv, "w") as f:
                    f.writelines(new_lines)

    return fixes


def fix_cmap_legacy_dir(pair_dir, pair_id, apply=False):
    """Move output_cmap_esm/ contents to output_cmaps/msa_transformer/."""
    old_dir = os.path.join(pair_dir, "output_cmap_esm")
    new_dir = os.path.join(pair_dir, "output_cmaps", "msa_transformer")

    if not os.path.isdir(old_dir):
        return None  # nothing to do

    files = [f for f in os.listdir(old_dir) if f.endswith((".npy", ".npz"))]
    if not files:
        return None

    # Check if new_dir already has files (don't overwrite)
    if os.path.isdir(new_dir):
        existing = set(os.listdir(new_dir))
        conflicts = set(files) & existing
        if conflicts:
            return ("CONFLICT", old_dir, new_dir, len(files), len(conflicts))

    if apply:
        os.makedirs(new_dir, exist_ok=True)
        for f in files:
            src = os.path.join(old_dir, f)
            dst = os.path.join(new_dir, f)
            shutil.move(src, dst)
        # Remove old directory if empty
        remaining = os.listdir(old_dir)
        if not remaining:
            os.rmdir(old_dir)

    return ("MOVE", old_dir, new_dir, len(files), 0)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--apply", action="store_true",
                        help="Actually apply fixes (default: dry-run, report only)")
    args = parser.parse_args()

    if not os.path.isdir(DATA_DIR):
        print(f"ERROR: DATA_DIR not found: {DATA_DIR}")
        sys.exit(1)

    pairs = sorted(d for d in os.listdir(DATA_DIR) if PAIR_RE.match(d))
    print(f"Scanning {len(pairs)} pairs in {DATA_DIR}")
    print(f"Mode: {'APPLY' if args.apply else 'DRY-RUN (use --apply to fix)'}\n")

    # --- Fix 1: samples_index.tsv ---
    print("=" * 60)
    print("FIX 1: samples_index.tsv stale paths")
    print("=" * 60)
    total_tsv_fixes = 0
    for pair_id in pairs:
        pair_dir = os.path.join(DATA_DIR, pair_id)
        fixes = fix_samples_index_tsv(pair_dir, pair_id, apply=args.apply)
        for tsv_path, n_lines in fixes:
            action = "FIXED" if args.apply else "WOULD FIX"
            print(f"  {action}: {tsv_path} ({n_lines} lines)")
            total_tsv_fixes += 1

    if total_tsv_fixes == 0:
        print("  All samples_index.tsv files are OK (no stale paths found)")
    else:
        print(f"\n  Total TSV files {'fixed' if args.apply else 'to fix'}: {total_tsv_fixes}")

    # --- Fix 2: output_cmap_esm/ -> output_cmaps/msa_transformer/ ---
    print()
    print("=" * 60)
    print("FIX 2: Legacy output_cmap_esm/ directories")
    print("=" * 60)
    total_cmap_moves = 0
    for pair_id in pairs:
        pair_dir = os.path.join(DATA_DIR, pair_id)
        result = fix_cmap_legacy_dir(pair_dir, pair_id, apply=args.apply)
        if result is None:
            continue
        status, old, new, n_files, n_conflicts = result
        if status == "CONFLICT":
            print(f"  CONFLICT: {pair_id} - {old} has {n_files} files but "
                  f"{new} already has {n_conflicts} conflicting files. SKIPPED.")
        else:
            action = "MOVED" if args.apply else "WOULD MOVE"
            print(f"  {action}: {pair_id} - {n_files} files from output_cmap_esm/ "
                  f"to output_cmaps/msa_transformer/")
            total_cmap_moves += 1

    if total_cmap_moves == 0:
        print("  No legacy output_cmap_esm/ directories found")
    else:
        print(f"\n  Total pairs {'moved' if args.apply else 'to move'}: {total_cmap_moves}")

    # --- Summary ---
    print()
    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"  TSV path fixes:  {total_tsv_fixes}")
    print(f"  Cmap dir moves:  {total_cmap_moves}")
    if not args.apply and (total_tsv_fixes or total_cmap_moves):
        print("\n  Re-run with --apply to execute these changes.")


if __name__ == "__main__":
    main()
