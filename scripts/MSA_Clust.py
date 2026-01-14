#!/usr/bin/env python3
"""
Legacy wrapper for MSA clustering.
Use run_ClusterMSA.py (single source of truth) instead of this script.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def main() -> int:
    repo_root = Path(__file__).resolve().parents[1]
    target = repo_root / "run_ClusterMSA.py"
    cmd = [sys.executable, str(target)] + sys.argv[1:]
    print("[deprecated] scripts/MSA_Clust.py -> run_ClusterMSA.py", flush=True)
    return subprocess.call(cmd)


if __name__ == "__main__":
    raise SystemExit(main())
