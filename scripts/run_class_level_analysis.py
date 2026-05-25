"""Thin CLI wrapper for Analysis.class_level.run_all -- see that module's
docstring for details. Run from the MsaCluster repo root:

    python3 scripts/run_class_level_analysis.py
    python3 scripts/run_class_level_analysis.py --n-perm 5000
    python3 scripts/run_class_level_analysis.py --include-interaction
"""
from __future__ import annotations
import os, sys

# Make `Analysis.class_level.run_all` importable when invoked as a script.
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_THIS_DIR))

from Analysis.class_level.run_all import main

if __name__ == "__main__":
    main()
