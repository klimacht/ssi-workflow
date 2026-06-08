#!/bin/bash -l
#SBATCH -o __SSI_ROOT__/logs/01_equilibrate.%j.out
#SBATCH -e __SSI_ROOT__/logs/01_equilibrate.%j.err
#SBATCH -J ssi_equil
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=8
#SBATCH --mem=400G
#SBATCH --time=12:00:00
#SBATCH --mail-type=none

set -euo pipefail

if [ $# -lt 1 ]; then
  echo "Usage: sbatch 01_equilibrate.slurm PDBID"
  exit 1
fi

PDBID=$(echo "$1" | tr '[:lower:]' '[:upper:]')
source __SSI_ROOT__/config/workflow.env

source "${SSI_ROOT}/bin/ssi_lib.sh"
ssi_load_modules
ssi_setup_env

RUN_DIR="${SSI_ROOT}/runs/${PDBID}"
PREP_DIR="${RUN_DIR}/prep"
EQUIL_DIR="${RUN_DIR}/equil"

mkdir -p "${EQUIL_DIR}" "${SSI_ROOT}/logs"
cd "${EQUIL_DIR}"

echo "============================================================"
echo "SSI equilibration job"
echo "PDBID                 = ${PDBID}"
echo "RUN_DIR               = ${RUN_DIR}"
echo "PREP_DIR              = ${PREP_DIR}"
echo "EQUIL_DIR             = ${EQUIL_DIR}"
echo "GMX                   = ${GMX}"
echo "SLURM_NTASKS          = ${SLURM_NTASKS}"
echo "SLURM_CPUS_PER_TASK   = ${SLURM_CPUS_PER_TASK}"
echo "OMP_NUM_THREADS       = ${OMP_NUM_THREADS}"
echo "============================================================"

echo
echo "[0] Check input files"
ssi_require "${PREP_DIR}/em.gro" "${PREP_DIR}/topol.top"
ssi_require_posre "${PREP_DIR}"

rm -f em.gro
cp -p "${PREP_DIR}/em.gro" em.gro
ssi_stage_topology "${PREP_DIR}"

echo
echo "Local topology tail:"
tail -n 30 topol.top

echo
echo "[1] NVT grompp"
mpirun -np 1 "${GMX}" grompp \
  -f "${SSI_ROOT}/mdp/nvt.mdp" \
  -c em.gro \
  -r em.gro \
  -p topol.top \
  -o nvt.tpr \
  -maxwarn 1

echo
echo "[2] NVT mdrun"
mpirun -np "${SLURM_NTASKS}" "${GMX}" mdrun \
  -v \
  -deffnm nvt \
  -ntomp "${SLURM_CPUS_PER_TASK}"

echo
echo "[3] NPT grompp"
mpirun -np 1 "${GMX}" grompp \
  -f "${SSI_ROOT}/mdp/npt.mdp" \
  -c nvt.gro \
  -t nvt.cpt \
  -r nvt.gro \
  -p topol.top \
  -o npt.tpr \
  -maxwarn 1

echo
echo "[4] NPT mdrun"
mpirun -np "${SLURM_NTASKS}" "${GMX}" mdrun \
  -v \
  -deffnm npt \
  -ntomp "${SLURM_CPUS_PER_TASK}"

echo
echo "[5] Summary"
ls -lh topol.top posre*.itp nvt.tpr nvt.gro nvt.cpt nvt.edr nvt.log nvt.xtc npt.tpr npt.gro npt.cpt npt.edr npt.log npt.xtc

echo
echo "============================================================"
echo "Equilibration finished successfully for ${PDBID}"
echo "============================================================"
