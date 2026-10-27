#!/bin/bash -l
#SBATCH -o __SSI_ROOT__/logs/dopa_qm.%A_%a.out
#SBATCH -e __SSI_ROOT__/logs/dopa_qm.%A_%a.err
#SBATCH -J dopa_qm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=60G
#SBATCH --time=08:00:00
#SBATCH --array=0-5
#SBATCH --mail-type=none

# HF/6-31G* optimization (phi/psi constrained) + ESP on MK points for one
# capped dipeptide per array task (PySCF + geomeTRIC from DOPA_QM_VENV).
# Fit afterwards with 11_resp_fit.slurm.

set -euo pipefail

source __SSI_ROOT__/config/workflow.env
source "${SSI_ROOT}/dopa/dopa.env"

CONFS=(dah_alpha_a dah_alpha_b dah_beta_a dah_beta_b tyr_alpha tyr_beta)
NAME="${CONFS[${SLURM_ARRAY_TASK_ID}]}"
WORK="${DOPA_RUNS}/qm/${NAME}"

module purge
module load "${MODULE_PYTHON}"
source "${DOPA_QM_VENV}/bin/activate"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

mkdir -p "${WORK}/tmp" "${SSI_ROOT}/logs"
export TMPDIR="${WORK}/tmp"
export PYSCF_TMPDIR="${WORK}/tmp"
cd "${WORK}"

echo "============================================================"
echo "DOPA RESP QM job"
echo "CONFORMER        = ${NAME}"
echo "WORK             = ${WORK}"
echo "PYTHON           = $(command -v python3)"
echo "OMP_NUM_THREADS  = ${OMP_NUM_THREADS}"
echo "============================================================"

python3 "${SSI_ROOT}/dopa/params/resp_dopa.py" qm \
  --name "${NAME}" --outdir "${WORK}" --memory-mb 50000

rm -rf "${WORK}/tmp"
echo
echo "QM done for ${NAME}"
