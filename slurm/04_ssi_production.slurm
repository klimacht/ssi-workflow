#!/bin/bash -l
#SBATCH -o __SSI_ROOT__/logs/04_ssi_production.%A_%a.out
#SBATCH -e __SSI_ROOT__/logs/04_ssi_production.%A_%a.err
#SBATCH -J ssi_prod
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=8
#SBATCH --mem=400G
#SBATCH --time=12:00:00
#SBATCH --mail-type=none

set -euo pipefail

if [ $# -lt 1 ]; then
  echo "Usage: sbatch --array=0-12 04_ssi_production.slurm PDBID"
  exit 1
fi

PDBID=$(echo "$1" | tr '[:lower:]' '[:upper:]')

source __SSI_ROOT__/config/workflow.env

source "${SSI_ROOT}/bin/ssi_lib.sh"
ssi_load_modules python
ssi_setup_env

# ---- Resolve phi from array index ----
PHI_ARRAY=(${PHI_LIST})
PHI="${PHI_ARRAY[${SLURM_ARRAY_TASK_ID}]}"

# phi dir name via shared helper: 0 -> phi_00p0, 5.5 -> phi_05p5
PHIDIR_NAME=$(ssi_phidir "${PHI}")

RUN_DIR="${SSI_ROOT}/runs/${PDBID}"
WORK_DIR="${RUN_DIR}/ssi_production/${PHIDIR_NAME}"

mkdir -p "${WORK_DIR}" "${SSI_ROOT}/logs"

echo "============================================================"
echo "SSI production job"
echo "PDBID                 = ${PDBID}"
echo "ARRAY_TASK_ID         = ${SLURM_ARRAY_TASK_ID}"
echo "PHI                   = ${PHI}"
echo "PHIDIR_NAME           = ${PHIDIR_NAME}"
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

cd "${WORK_DIR}"

# ---- Checkpoint restart ----
if [ -f "ssi_phi.cpt" ] && [ -f "ssi_phi.tpr" ] && [ -f "plumed_runtime.dat" ]; then
  echo
  echo "[1] Checkpoint found — restarting"
  ls -lh ssi_phi.cpt ssi_phi.tpr plumed_runtime.dat

  mpirun -np "${SLURM_NTASKS}" "${GMX}" mdrun \
    -v \
    -s ssi_phi.tpr \
    -cpi ssi_phi.cpt \
    -deffnm ssi_phi \
    -plumed plumed_runtime.dat \
    -append \
    -ntomp "${SLURM_CPUS_PER_TASK}"

else
  echo
  echo "[1] No checkpoint — running setup"

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

  echo
  echo "[2] File summary after setup"
  ls -lh

  echo
  echo "[3] SSI setup summary"
  cat ssi_setup_summary.txt

  echo
  echo "[4] MDrun"
  mpirun -np "${SLURM_NTASKS}" "${GMX}" mdrun \
    -v \
    -s ssi_phi.tpr \
    -deffnm ssi_phi \
    -plumed plumed_runtime.dat \
    -ntomp "${SLURM_CPUS_PER_TASK}"
fi

echo
echo "[5] Output check"
ls -lh ssi_phi.gro ssi_phi.cpt ssi_phi.edr ssi_phi.log ssi_phi.xtc COLVAR

echo
echo "[6] COLVAR tail"
tail -n 5 COLVAR

echo
echo "============================================================"
echo "SSI production finished: ${PDBID} ${PHIDIR_NAME}"
echo "============================================================"
