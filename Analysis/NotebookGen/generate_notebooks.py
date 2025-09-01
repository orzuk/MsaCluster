
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import argparse
import subprocess
from pathlib import Path
from typing import Iterable, Optional

# ---------- repo imports ----------
REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from config import *
from utils.utils import list_protein_pairs


# ---------- worker: one pair -> one HTML ----------
def generate_html_for_pair(
    pair: str,
    output_html: Path,
    template_notebook: Path,
    kernel_name: Optional[str] = None,
    timeout_seconds: int = 1800,
) -> int:
    """
    Execute a template notebook for a given protein pair and export to a single HTML.

    Parameters
    ----------
    pair : str
        Pair identifier, e.g. "4n9wA_4nc9C".
    output_html : Path
        Full output path ending with .html (directory is created if needed).
    template_notebook : Path
        Path to the template .ipynb to execute.
    kernel_name : Optional[str]
        Jupyter kernel name to use. If None, uses the notebook's recorded kernel.
    timeout_seconds : int
        Execution timeout passed to nbconvert ExecutePreprocessor.

    Returns
    -------
    int
        Return code from the nbconvert process (0 = success).
    """
    output_html = Path(output_html)
    output_html.parent.mkdir(parents=True, exist_ok=True)

    env = os.environ.copy()
    env["PAIR_ID"] = pair

    cmd = [
        sys.executable, "-m", "jupyter", "nbconvert",
        "--to", "html",
        "--execute", str(template_notebook),
        f"--ExecutePreprocessor.timeout={timeout_seconds}",
        "--output", output_html.stem,
        "--output-dir", str(output_html.parent),
    ]
    if kernel_name:
        cmd.append(f"--ExecutePreprocessor.kernel_name={kernel_name}")

    print(f"Executing {template_notebook.name} for pair '{pair}' -> {output_html.name}")
    return subprocess.run(cmd, env=env, text=True).returncode


# ---------- CLI / main ----------
def _resolve_pairs(args_pairs: Iterable[str]) -> list[str]:
    DEFAULT_PAIRS = ["4n9wA_4nc9C", "1nrjB_2gedB"]
    pairs = list(args_pairs)
    if len(pairs) == 1 and pairs[0].upper() == "ALL":
        return list_protein_pairs()
    if not pairs:
        # default: take from DEFAULT_PAIRS if provided; else all pairs
        default_pairs = getattr(DEFAULT_PAIRS)
        return list(default_pairs) if default_pairs else list_protein_pairs()
    return pairs


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate per-pair HTML reports by executing a template notebook."
    )
    parser.add_argument(
        "pairs",
        nargs="*",
        help="Protein pairs (e.g., 4n9wA_4nc9C 1nrjB_2gedB) or 'ALL' to run all.",
    )
    args = parser.parse_args()

    # Paths & settings from config.py
    TEMPLATE_NOTEBOOK = Path(NOTEBOOK_TEMPLATE)   # e.g., MAIN_DIR / 'TemplateNotebook.ipynb'
    HTML_DIR = Path(HTML_OUTPUT_DIR)              # e.g., MAIN_DIR / 'HTML'
    KERNEL_NAME = NOTEBOOK_KERNEL_NAME

    pairs = _resolve_pairs(args.pairs)

    print("Current directory:", Path.cwd())
    print("Start generate Notebook!!!")
    print("Pairs:", pairs)
    print("Template:", TEMPLATE_NOTEBOOK)
    print("HTML out dir:", HTML_DIR)
    if KERNEL_NAME:
        print("Kernel:", KERNEL_NAME)

    failures: list[str] = []
    for pair in pairs:
        print("#" * 74)
        print(pair)
        print("#" * 74)
        out_html = HTML_DIR / f"{pair}.html"
        rc = generate_html_for_pair(
            pair=pair,
            output_html=out_html,
            template_notebook=TEMPLATE_NOTEBOOK,
            kernel_name=KERNEL_NAME,
        )
        if rc == 0:
            print(f"OK: {out_html}")
        else:
            print(f"FAILED (rc={rc}): {pair}")
            failures.append(pair)

    print("Finish to run !")
    if failures:
        # Non-zero exit to surface failures in CI/IDE run panel
        sys.exit(1)


if __name__ == "__main__":
    main()

