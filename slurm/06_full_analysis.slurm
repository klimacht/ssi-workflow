#!/bin/bash -l
#SBATCH -o __SSI_ROOT__/logs/06_full_analysis.%j.out
#SBATCH -e __SSI_ROOT__/logs/06_full_analysis.%j.err
#SBATCH -J ssi_fullana
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=400G
#SBATCH --time=04:00:00
#SBATCH --mail-type=none

# Full per-atom dewetting analysis using ALL phi windows.
#
#   1. compute_n0           (if missing)
#   2. compute_nS for every phi window (if missing)
#   3. compute_N_vs_phi     (susceptibility, phi*)
#   4. combine_dewetting_map.py            (single-window eta map at phi*)
#   5. combine_dewetting_phistar_map.py    (per-atom transition map, all windows)
#
# Usage:
#   sbatch slurm/06_full_analysis.slurm PHLP5

set -euo pipefail

if [ $# -lt 1 ]; then
  echo "Usage: sbatch 06_full_analysis.slurm PDBID"
  exit 1
fi

PDBID=$(echo "$1" | tr '[:lower:]' '[:upper:]')

source __SSI_ROOT__/config/workflow.env
source "${SSI_ROOT}/bin/ssi_lib.sh"
ssi_load_modules python
ssi_setup_env

SRC="${SSI_ROOT}/src"
ANA_DIR="${SSI_ROOT}/runs/${PDBID}/analysis"
PROD_DIR="${SSI_ROOT}/runs/${PDBID}/ssi_production"
mkdir -p "${ANA_DIR}" "${SSI_ROOT}/logs"

T_START_PS=$(awk "BEGIN{printf \"%d\", ${DISCARD_NS}*1000}")

echo "============================================================"
echo "SSI full per-atom analysis"
echo "PDBID        = ${PDBID}"
echo "PHI_LIST     = ${PHI_LIST}"
echo "t_start      = ${T_START_PS} ps"
echo "radius       = ${SSI_RADIUS_NM} nm"
echo "============================================================"

# ---- 1. n0 ----
echo
echo "[1] n0 (unbiased per-atom hydration)"
if [ -f "${ANA_DIR}/n0_per_atom.npy" ]; then
  echo "    already present — skipping"
else
  python3 "${SRC}/compute_n0.py" \
    --pdbid "${PDBID}" --root "${SSI_ROOT}" \
    --t-start "${T_START_PS}" --radius "${SSI_RADIUS_NM}"
fi

# ---- 2. nS for every phi window ----
echo
echo "[2] nS for every phi window (parallel, one process per window)"
# Each compute_nS is single-threaded; run all windows concurrently on the
# reserved node. Limit BLAS/OMP threads per process to avoid oversubscription.
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# Build the list of windows that still need computing.
todo=()
for PHI in ${PHI_LIST}; do
  PHIDIR=$(ssi_phidir "${PHI}")
  if [ -f "${ANA_DIR}/nS_${PHIDIR}.npy" ]; then
    echo "    [skip] ${PHIDIR}"
    continue
  fi
  if [ ! -f "${PROD_DIR}/${PHIDIR}/ssi_phi.xtc" ]; then
    echo "    [warn] ${PHIDIR}: no trajectory, skipping"
    continue
  fi
  todo+=( "${PHIDIR}" )
done

if [ ${#todo[@]} -gt 0 ]; then
  echo "    launching ${#todo[@]} windows in parallel: ${todo[*]}"
  # Use up to as many concurrent jobs as windows (<=13), each pinned to 1 thread.
  printf '%s\n' "${todo[@]}" | xargs -P "${#todo[@]}" -I {} \
    bash -c '
      PHIDIR="$1"
      echo "    [run ] ${PHIDIR} (pid $$)"
      python3 "'"${SRC}"'/compute_nS.py" \
        --pdbid "'"${PDBID}"'" --phi "${PHIDIR}" --root "'"${SSI_ROOT}"'" \
        --t-start "'"${T_START_PS}"'" --radius "'"${SSI_RADIUS_NM}"'" \
        > "'"${ANA_DIR}"'/nS_${PHIDIR}.log" 2>&1 \
      && echo "    [done] ${PHIDIR}" \
      || echo "    [FAIL] ${PHIDIR} (see nS_${PHIDIR}.log)"
    ' _ {}
  echo "    all parallel nS jobs finished"
else
  echo "    nothing to compute — all windows present"
fi

# ---- 3. N vs phi + susceptibility + phi* ----
echo
echo "[3] N vs phi / susceptibility"
python3 "${SRC}/compute_N_vs_phi.py" \
  --pdbid "${PDBID}" --root "${SSI_ROOT}" \
  --t-start "${T_START_PS}" --phi-list "${PHI_LIST}"

# Determine phi* dirname from the TSV (phi_star column == 1)
PHI_STAR_DIR=$(python3 - "${ANA_DIR}/N_vs_phi.tsv" <<'PYEOF'
import csv, sys
tsv = sys.argv[1]
res = None
with open(tsv) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        if row.get("phi_star", "0") == "1":
            res = row["dirname"]; break
if res is None:
    raise SystemExit("phi_star not found in N_vs_phi.tsv")
print(res)
PYEOF
)
echo "    phi* = ${PHI_STAR_DIR}"

# ---- 4. single-window eta map (reference method) ----
echo
echo "[4] single-window eta map at ${PHI_STAR_DIR}"
python3 "${SRC}/combine_dewetting_map.py" \
  --pdbid "${PDBID}" --phi "${PHI_STAR_DIR}" --root "${SSI_ROOT}" \
  --surface-threshold 4.0

# ---- 5. per-atom transition map (all windows) ----
echo
echo "[5] per-atom transition (phi*_i) map using all windows"
python3 "${SRC}/combine_dewetting_phistar_map.py" \
  --pdbid "${PDBID}" --root "${SSI_ROOT}" \
  --surface-threshold 4.0 --phi-list "${PHI_LIST}"

echo
echo "============================================================"
echo "Full analysis complete for ${PDBID}"
echo "Outputs in runs/${PDBID}/dewetting/:"
echo "  dewetting_${PHI_STAR_DIR}.pdb   (single-window eta)"
echo "  dewetting_phistar.pdb           (per-atom transition, all windows)"
echo "============================================================"
