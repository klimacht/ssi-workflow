#!/usr/bin/env bash
# Master analysis script for SSI workflow.
# Runs all analysis steps in order after production jobs complete.
#
# Usage:
#   bash run_analysis.sh 1DPX [phi_star]
#   bash run_analysis.sh 1DPX          # auto-detects phi* from N_vs_phi.tsv
#   bash run_analysis.sh 1DPX phi_06p0 # use given phi*

set -euo pipefail

ROOT="__SSI_ROOT__"
SRC="${ROOT}/src"

PDBID=$(echo "${1:-}" | tr '[:lower:]' '[:upper:]')
if [ -z "$PDBID" ]; then
  echo "Usage: $0 PDBID [phi_star]"
  exit 1
fi

PHI_STAR_ARG="${2:-}"

source "${ROOT}/config/workflow.env"
module load "${MODULE_PYTHON}" 2>/dev/null || true

RUN_DIR="${ROOT}/runs/${PDBID}"
ANA_DIR="${RUN_DIR}/analysis"
mkdir -p "${ANA_DIR}"
T_START_PS=$(awk "BEGIN{printf \"%d\", ${DISCARD_NS}*1000}")

echo "=========================================="
echo "SSI Analysis: ${PDBID}"
echo "$(date)"
echo "=========================================="

# ---- Step 1: N vs phi from COLVAR files (fast) ----
echo
echo "[1] Computing N vs phi from COLVAR files..."
python3 "${SRC}/compute_N_vs_phi.py" \
  --pdbid "${PDBID}" \
  --root  "${ROOT}" \
  --t-start "${T_START_PS}" \
  --phi-list "${PHI_LIST}"

# ---- Determine phi* ----
if [ -n "$PHI_STAR_ARG" ]; then
  PHI_STAR_DIR="$PHI_STAR_ARG"
  echo
  echo "Using user-supplied phi* = ${PHI_STAR_DIR}"
else
  # Pick phi with max N_var from TSV
  PHI_STAR_DIR=$(python3 - <<PY
import csv
tsv = "${ANA_DIR}/N_vs_phi.tsv"
result = None
with open(tsv) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        if row.get("phi_star", "0") == "1":
            result = row["dirname"]
            break
if result is None:
    raise SystemExit("phi_star not found in N_vs_phi.tsv")
print(result)
PY
)
  echo
  echo "Auto-detected phi* = ${PHI_STAR_DIR}"
fi

# ---- Step 2: Compute n0 from unbiased trajectory ----
echo
echo "[2] Computing n0 (unbiased per-atom hydration)..."
if [ -f "${ANA_DIR}/n0_per_atom.npy" ]; then
  echo "  n0_per_atom.npy already exists — skipping."
else
  python3 "${SRC}/compute_n0.py" \
    --pdbid   "${PDBID}" \
    --root    "${ROOT}" \
    --t-start "${T_START_PS}" \
    --radius  "${SSI_RADIUS_NM}"
fi

# ---- Step 3: Compute nS at phi* ----
echo
echo "[3] Computing nS at phi* = ${PHI_STAR_DIR}..."
if [ -f "${ANA_DIR}/nS_${PHI_STAR_DIR}.npy" ]; then
  echo "  nS_${PHI_STAR_DIR}.npy already exists — skipping."
else
  python3 "${SRC}/compute_nS.py" \
    --pdbid   "${PDBID}" \
    --phi     "${PHI_STAR_DIR}" \
    --root    "${ROOT}" \
    --t-start "${T_START_PS}" \
    --radius  "${SSI_RADIUS_NM}"
fi

# ---- Step 4: Combine into dewetting map ----
echo
echo "[4] Building dewetting map..."
python3 "${SRC}/combine_dewetting_map.py" \
  --pdbid              "${PDBID}" \
  --phi                "${PHI_STAR_DIR}" \
  --root               "${ROOT}" \
  --surface-threshold  4.0

echo
echo "=========================================="
echo "Analysis complete for ${PDBID}"
echo "phi* used: ${PHI_STAR_DIR}"
echo
echo "Outputs:"
echo "  ${ANA_DIR}/N_vs_phi.tsv"
echo "  ${ANA_DIR}/n0_per_atom.npy"
echo "  ${ANA_DIR}/nS_${PHI_STAR_DIR}.npy"
echo "  ${RUN_DIR}/dewetting/dewetting_${PHI_STAR_DIR}.dat"
echo "  ${RUN_DIR}/dewetting/dewetting_${PHI_STAR_DIR}.pdb"
echo "=========================================="
