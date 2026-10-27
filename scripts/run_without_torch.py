"""Run a repo tree post-processing script without torch installed.

Usage: QT_QPA_PLATFORM=offscreen python3 scripts/run_without_torch.py scripts/phylo_placement.py [args]


utils.plot_contacts imports torch at module level only to reach
seqs_ids_to_cluster_ids; stub that module with the same function so the tree
code runs without torch.
"""
import os, runpy, sys, types
from glob import glob

sys.path.insert(0, os.getcwd())


def _ids(path):
    with open(path) as fh:
        return [ln[1:].split()[0] for ln in fh if ln.startswith(">") and ln[1:].strip()]


def seqs_ids_to_cluster_ids(msa_cluster_dir, seqs_ids=None):
    files = (glob(os.path.join(msa_cluster_dir, "*.a3m")) if msa_cluster_dir.endswith("/")
             else glob(msa_cluster_dir))
    cluster_ids = {}
    for f in files:
        tag = os.path.basename(f).replace(".a3m", "")[7:]
        cluster_ids.update({s: tag for s in _ids(f)})
    if not seqs_ids:
        return cluster_ids
    return {s: cluster_ids.get(s, "p") for s in seqs_ids}


stub = types.ModuleType("utils.plot_contacts")
stub.seqs_ids_to_cluster_ids = seqs_ids_to_cluster_ids
sys.modules["utils.plot_contacts"] = stub
script = sys.argv[1]
sys.argv = sys.argv[1:]
runpy.run_path(script, run_name="__main__")
