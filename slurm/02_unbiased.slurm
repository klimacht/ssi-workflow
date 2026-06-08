#!/bin/bash -l
#SBATCH -o __SSI_ROOT__/logs/02_unbiased.%j.out
#SBATCH -e __SSI_ROOT__/logs/02_unbiased.%j.err
#SBATCH -J ssi_unbiased
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=8
#SBATCH --mem=400G
#SBATCH --time=12:00:00
#SBATCH --mail-type=none

set -euo pipefail

if [ $# -lt 1 ]; then
  echo "Usage: sbatch 02_unbiased.slurm PDBID"
  exit 1
fi

PDBID=$(echo "$1" | tr '[:lower:]' '[:upper:]')
source __SSI_ROOT__/config/workflow.env

source "${SSI_ROOT}/bin/ssi_lib.sh"
ssi_load_modules
ssi_setup_env

RUN_DIR="${SSI_ROOT}/runs/${PDBID}"
EQUIL_DIR="${RUN_DIR}/equil"
UNBIASED_DIR="${RUN_DIR}/unbiased"

mkdir -p "${UNBIASED_DIR}" "${SSI_ROOT}/logs"
cd "${UNBIASED_DIR}"

echo "============================================================"
echo "SSI unbiased production job"
echo "PDBID                 = ${PDBID}"
echo "RUN_DIR               = ${RUN_DIR}"
echo "EQUIL_DIR             = ${EQUIL_DIR}"
echo "UNBIASED_DIR          = ${UNBIASED_DIR}"
echo "GMX                   = ${GMX}"
echo "SLURM_NTASKS          = ${SLURM_NTASKS}"
echo "SLURM_CPUS_PER_TASK   = ${SLURM_CPUS_PER_TASK}"
echo "OMP_NUM_THREADS       = ${OMP_NUM_THREADS}"
echo "============================================================"

echo
echo "[0] Check input files"
ssi_require "${EQUIL_DIR}/npt.gro" "${EQUIL_DIR}/npt.cpt" "${EQUIL_DIR}/topol.top"
ssi_require_posre "${EQUIL_DIR}"

rm -f npt.gro npt.cpt
cp -p "${EQUIL_DIR}/npt.gro" npt.gro
cp -p "${EQUIL_DIR}/npt.cpt" npt.cpt
ssi_stage_topology "${EQUIL_DIR}"

echo
echo "Local topology tail:"
tail -n 30 topol.top

echo
echo "[1] Unbiased production grompp"
mpirun -np 1 "${GMX}" grompp \
  -f "${SSI_ROOT}/mdp/prod_unbiased.mdp" \
  -c npt.gro \
  -t npt.cpt \
  -r npt.gro \
  -p topol.top \
  -o prod_unbiased.tpr \
  -maxwarn 1

echo
echo "[2] Unbiased production mdrun"
if [ -f prod_unbiased.cpt ]; then
  echo "Checkpoint exists: restarting with -cpi prod_unbiased.cpt"
  mpirun -np "${SLURM_NTASKS}" "${GMX}" mdrun \
    -v \
    -s prod_unbiased.tpr \
    -cpi prod_unbiased.cpt \
    -deffnm prod_unbiased \
    -append \
    -ntomp "${SLURM_CPUS_PER_TASK}"
else
  echo "No checkpoint found: starting new run"
  mpirun -np "${SLURM_NTASKS}" "${GMX}" mdrun \
    -v \
    -s prod_unbiased.tpr \
    -deffnm prod_unbiased \
    -append \
    -ntomp "${SLURM_CPUS_PER_TASK}"
fi

echo
echo "[3] Summary"
ls -lh topol.top posre*.itp prod_unbiased.tpr prod_unbiased.gro prod_unbiased.cpt prod_unbiased.edr prod_unbiased.log prod_unbiased.xtc

echo
echo "[4] Last lines of log"
tail -n 40 prod_unbiased.log

echo
echo "============================================================"
echo "Unbiased production finished successfully for ${PDBID}"
echo "============================================================"
