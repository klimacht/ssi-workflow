#!/bin/bash -l
#SBATCH -o __SSI_ROOT__/logs/03_ssi_preflight.%j.out
#SBATCH -e __SSI_ROOT__/logs/03_ssi_preflight.%j.err
#SBATCH -J ssi_preflight
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=8
#SBATCH --mem=400G
#SBATCH --time=02:00:00
#SBATCH --mail-type=none

set -euo pipefail

if [ $# -lt 1 ]; then
  echo "Usage: sbatch 03_ssi_preflight.slurm PDBID [PHI]"
  echo "Example: sbatch 03_ssi_preflight.slurm 1DPX 0"
  exit 1
fi

PDBID=$(echo "$1" | tr '[:lower:]' '[:upper:]')
PHI="${2:-0}"

source __SSI_ROOT__/config/workflow.env

source "${SSI_ROOT}/bin/ssi_lib.sh"
ssi_load_modules python
ssi_setup_env

RUN_DIR="${SSI_ROOT}/runs/${PDBID}"
WORK_DIR="${RUN_DIR}/ssi_preflight/phi_00p0"

mkdir -p "${WORK_DIR}" "${SSI_ROOT}/logs"

echo "============================================================"
echo "SSI preflight job"
echo "PDBID                 = ${PDBID}"
echo "PHI                   = ${PHI}"
echo "RUN_DIR               = ${RUN_DIR}"
echo "WORK_DIR              = ${WORK_DIR}"
echo "GMX                   = ${GMX}"
echo "UNION_MODE            = ${SSI_UNION_MODE}"
echo "BIAS_MODE             = ${SSI_BIAS_MODE}"
echo "RADIUS_NM             = ${SSI_RADIUS_NM}"
echo "CELL_SIZE_NM          = ${SSI_CELL_SIZE_NM}"
echo "SLURM_NTASKS          = ${SLURM_NTASKS}"
echo "SLURM_CPUS_PER_TASK   = ${SLURM_CPUS_PER_TASK}"
echo "OMP_NUM_THREADS       = ${OMP_NUM_THREADS}"
echo "============================================================"

echo
echo "[0] Check unbiased input"
ssi_require "${RUN_DIR}/unbiased/prod_unbiased.gro" \
            "${RUN_DIR}/unbiased/prod_unbiased.cpt" \
            "${RUN_DIR}/unbiased/topol.top"
ssi_require_posre "${RUN_DIR}/unbiased"

echo
echo "[1] Run SSI setup script"
python3 "${SSI_ROOT}/src/setup_ssi_case.py" \
  --pdbid "${PDBID}" \
  --run-dir "${RUN_DIR}" \
  --work-dir "${WORK_DIR}" \
  --gmx "${GMX}" \
  --phi "${PHI}" \
  --radius-nm "${SSI_RADIUS_NM}" \
  --cell-size-nm "${SSI_CELL_SIZE_NM}" \
  --union-mode "${SSI_UNION_MODE}" \
  --bias-mode "${SSI_BIAS_MODE}" \
  --mdp "${SSI_ROOT}/mdp/ssi_biased.mdp"

cd "${WORK_DIR}"

echo
echo "[2] File summary"
ls -lh

echo
echo "[3] SSI setup summary"
cat ssi_setup_summary.txt

echo
echo "[4] INDUS input"
cat indus_union_runtime.input

echo
echo "[5] PLUMED input"
cat plumed_runtime.dat

echo
echo "[6] Centers file sanity"
echo "Number of centers:"
wc -l heavy_atom_centers_fixed.dat
echo "First 5 centers:"
head -n 5 heavy_atom_centers_fixed.dat
echo "Last 5 centers:"
tail -n 5 heavy_atom_centers_fixed.dat

echo
echo "============================================================"
echo "SSI preflight finished successfully for ${PDBID}"
echo "============================================================"
