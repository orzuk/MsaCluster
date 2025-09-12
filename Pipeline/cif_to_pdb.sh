#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <run_dir> [--mode all|rank1]"; exit 2
fi

RUN_DIR="$1"; shift || true
MODE="all"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --mode) MODE="$2"; shift 2 ;;
    --mode=*) MODE="${1#--mode=}"; shift ;;
    *) echo "[fatal] unknown arg: $1"; exit 2 ;;
  esac
done

if [[ ! -d "$RUN_DIR" ]]; then
  echo "[fatal] run_dir not found: $RUN_DIR"; exit 2
fi

convert_one() {
  local cif="$1" pdb="$2"
  if command -v gemmi >/dev/null 2>&1; then
    gemmi convert "$cif" "$pdb"
  else
    python - "$cif" "$pdb" <<'PY'
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBIO import PDBIO
import sys
cif, pdb = sys.argv[1], sys.argv[2]
parser = MMCIFParser(QUIET=True)
s = parser.get_structure("model", cif)
io = PDBIO(); io.set_structure(s); io.save(pdb)
PY
  fi
}

# get job name from ranking CSV if present
CSV="$(ls -1 "$RUN_DIR"/*_ranking_scores.csv 2>/dev/null | head -n1 || true)"
JOB=""
if [[ -n "$CSV" ]]; then
  JOB="$(basename "$CSV" | sed 's/_ranking_scores\.csv$//')"
fi

if [[ "$MODE" == "rank1" ]]; then
  if [[ -z "$CSV" ]]; then
    echo "[warn] ranking CSV not found in $RUN_DIR; skipping rank1 export"; exit 0
  fi
  BEST_MODEL="$(python - "$CSV" <<'PY'
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1])))
keys = ("ranking_score","mean_plddt","plddt","confidence")
k = next((x for x in keys if x in rows[0]), list(rows[0].keys())[-1])
best = max(rows, key=lambda r: float(r[k]))
print(best.get("model","seed-1_sample-0"))
PY
)"
  CIF="$RUN_DIR/$BEST_MODEL/${JOB}_${BEST_MODEL}_model.cif"
  PDB="$RUN_DIR/$BEST_MODEL/${JOB}_${BEST_MODEL}_model.pdb"
  if [[ -f "$CIF" ]]; then
    convert_one "$CIF" "$PDB" && echo "[ok] wrote $PDB"
    cp -f "$PDB" "$RUN_DIR/${JOB}_rank1_model.pdb"
  else
    echo "[warn] best CIF not found: $CIF"
  fi
else
  found=0
  while IFS= read -r -d '' CIF; do
    PDB="${CIF%.cif}.pdb"
    convert_one "$CIF" "$PDB" && echo "[ok] wrote $PDB"
    found=1
  done < <(find "$RUN_DIR" -type f -name "*_model.cif" -print0)
  if [[ $found -eq 0 ]]; then
    echo "[warn] no *_model.cif files under $RUN_DIR"
  fi
fi
