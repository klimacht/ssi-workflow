#!/bin/bash -l
#SBATCH -o __SSI_ROOT__/logs/%x.%j.out
#SBATCH -e __SSI_ROOT__/logs/%x.%j.err
#SBATCH -J dopa_qm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=400G
#SBATCH --time=08:00:00
#SBATCH --mail-type=none

# HF/6-31G* optimization (phi/psi constrained) + ESP on MK points for the
# capped dipeptides, all conformers in parallel on one node (PySCF + geomeTRIC
# from DOPA_QM_VENV). Fit afterwards with 11_resp_fit.slurm.
#
#   sbatch 10_resp_qm.slurm            all conformers (TYR, DAH, TOP)
#   sbatch 10_resp_qm.slurm TOP        only the conformers of one residue
# Already optimized conformers are reused (only SCF + ESP are recomputed).

set -euo pipefail

RES_FILTER="${1:-}"

source __SSI_ROOT__/dopa/bin/dopa_env.sh __SSI_ROOT__
RESP="${SSI_ROOT}/dopa/params/resp_dopa.py"

module purge
module load "${MODULE_PYTHON}"
source "${DOPA_QM_VENV}/bin/activate"

mapfile -t CONFS < <(python3 "${RESP}" list ${RES_FILTER:+--res "${RES_FILTER}"})
N=${#CONFS[@]}
THREADS=$(( SLURM_CPUS_PER_TASK / N ))
[ "${THREADS}" -ge 1 ] || THREADS=1
MEM_MB=$(( 380000 / N ))
mkdir -p "${DOPA_RUNS}/qm" "${SSI_ROOT}/logs"

echo "============================================================"
echo "DOPA RESP QM job"
echo "RESIDUES        = ${RES_FILTER:-all}"
echo "CONFORMERS      = ${CONFS[*]}"
echo "THREADS / conf  = ${THREADS}   memory / conf = ${MEM_MB} MB"
echo "PYTHON          = $(command -v python3)"
echo "============================================================"

pids=()
for NAME in "${CONFS[@]}"; do
  WORK="${DOPA_RUNS}/qm/${NAME}"
  mkdir -p "${WORK}/tmp"
  (
    cd "${WORK}"
    export OMP_NUM_THREADS="${THREADS}" TMPDIR="${WORK}/tmp" PYSCF_TMPDIR="${WORK}/tmp"
    python3 "${RESP}" qm --name "${NAME}" --outdir "${WORK}" --memory-mb "${MEM_MB}" > "${WORK}/qm.log" 2>&1
    rm -rf "${WORK}/tmp"
  ) &
  pids+=($!)
done

fail=0
for i in "${!pids[@]}"; do
  if wait "${pids[$i]}"; then
    tail -n 1 "${DOPA_RUNS}/qm/${CONFS[$i]}/qm.log"
  else
    echo "ERROR: ${CONFS[$i]} failed, see ${DOPA_RUNS}/qm/${CONFS[$i]}/qm.log"
    tail -n 20 "${DOPA_RUNS}/qm/${CONFS[$i]}/qm.log"
    fail=1
  fi
done
exit "${fail}"
