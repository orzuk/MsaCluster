#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import argparse
import subprocess
from pathlib import Path
from typing import Optional, Tuple, Union, Iterable

# --- repo imports (top-level only) ---
REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from config import *            # expects MAIN_DIR, etc.
from utils.utils import list_protein_pairs
from utils.protein_utils import normalize_pair



# ---------- worker: one pair -> one HTML ----------
def generate_html_for_pair(
    pair: Union[Tuple[str, str], str],
    output_html: Path,
    template_notebook: Path,
    kernel_name: Optional[str] = None,
    timeout_seconds: int = 1800,
) -> int:
    """
    Execute a template notebook for a given protein pair and export a single HTML.
    """
    pair_a, pair_b, pair_id = normalize_pair(pair)

    output_html = Path(output_html)
    output_html.parent.mkdir(parents=True, exist_ok=True)

    # Env must be strings on Windows
    env = os.environ.copy()
    env["PAIR_ID"] = pair_id
    env["PAIR_A"] = pair_a
    env["PAIR_B"] = pair_b

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

    print(f"Executing {template_notebook.name} for '{pair_id}' -> {output_html.name}")
    return subprocess.run(cmd, env=env, text=True).returncode


# ---------- CLI / main ----------
def _resolve_pairs(args_pairs: Iterable[Union[str, Tuple[str, str]]]) -> list[Union[str, Tuple[str, str]]]:
    pairs = list(args_pairs)
    if len(pairs) == 1 and isinstance(pairs[0], str) and pairs[0].upper() == "ALL":
        return list_protein_pairs()
    if not pairs:
        # DEFAULT: run ONE pair (the first returned by utils.list_protein_pairs())
        all_pairs = list_protein_pairs()
        return [all_pairs[0]] if all_pairs else []
    return pairs


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate per-pair HTML reports by executing a template notebook."
    )
    parser.add_argument(
        "pairs",
        nargs="*",
        help="Pairs like 4n9wA_4nc9C (or 'ALL' to run all pairs). If no args: run ONE default pair.",
    )
    parser.add_argument(
        "--kernel",
        default=os.environ.get("JUPYTER_KERNEL_NAME", "python3"),
        help="Jupyter kernel name to use (default: 'python3').",
    )
    args = parser.parse_args()

    main_dir = Path(MAIN_DIR)
    template_notebook = main_dir / "TemplateNotebook.ipynb"
    html_dir = main_dir / "HTML"
    html_dir.mkdir(parents=True, exist_ok=True)

    pairs = _resolve_pairs(args.pairs)

#    user = os.environ.get("USERNAME") or os.environ.get("USER") or "unknown"
#    env_name = "Windows" if os.name == "nt" else "Unix"
    print("Current directory: ", Path(__file__).parent)
    print("RUNNING FOLD-SWITCH GENERATE NOTEBOOK WITH USER:", user, " ENVIRONMENT:", platform.system())
    print("Current directory:", Path.cwd())
    print("Start generate Notebook!!!")
    print("Pairs:", pairs)
    print("Template:", template_notebook)
    print("HTML out dir:", html_dir)
    print("Kernel:", args.kernel)

    failures: list[str] = []
    for pair in pairs:
        _, _, pair_id = normalize_pair(pair)
        print("#" * 74)
        print(pair_id)
        print("#" * 74)
        out_html = html_dir / f"{pair_id}.html"

        rc = generate_html_for_pair(
            pair=pair,
            output_html=out_html,
            template_notebook=template_notebook,
            kernel_name=args.kernel)
        if rc == 0:
            print(f"OK: {out_html}")
        else:
            print(f"FAILED (rc={rc}): {pair_id}")
            failures.append(pair_id)

    print("Finish to run !")
    if failures:
        sys.exit(1)


if __name__ == "__main__":
    main()
