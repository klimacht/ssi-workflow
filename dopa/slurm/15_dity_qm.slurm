#!/bin/bash -l
#SBATCH -o __SSI_ROOT__/logs/%x.%j.out
#SBATCH -e __SSI_ROOT__/logs/%x.%j.err
#SBATCH -J dopa_dity_qm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=400G
#SBATCH --time=08:00:00
#SBATCH --mail-type=none

# 3,3'-dityrosine parameters: torsion scan (MP2//HF/6-31G*) and RESP conformers of
# the model compounds, all in parallel on one node, then the fit
# (dopa/params/dity_params.py) -> dopa/ff/dty.*, dity_*.itp/.atp, params/results/dity_*.

set -euo pipefail

source __SSI_ROOT__/dopa/bin/dopa_env.sh __SSI_ROOT__
DITY="${SSI_ROOT}/dopa/params/dity_params.py"
WORK="${DOPA_RUNS}/qm_dity"
FFBASE="$(dirname "$(readlink -f "${GMX}")")/../share/gromacs/top/${DOPA_BASE_FF}.ff"

module purge
module load "${MODULE_PYTHON}"
source "${DOPA_QM_VENV}/bin/activate"

mapfile -t TASKS < <(python3 "${DITY}" list)
N=${#TASKS[@]}
THREADS=$(( SLURM_CPUS_PER_TASK / N ))
[ "${THREADS}" -ge 1 ] || THREADS=1
MEM_MB=$(( 380000 / N ))
mkdir -p "${WORK}" "${SSI_ROOT}/logs"

echo "============================================================"
echo "Dityrosine QM job: ${N} tasks, ${THREADS} threads and ${MEM_MB} MB each"
echo "WORK   = ${WORK}"
echo "FFBASE = ${FFBASE}"
echo "============================================================"

pids=()
for T in "${TASKS[@]}"; do
  D="${WORK}/${T}"
  mkdir -p "${D}/tmp"
  (
    cd "${D}"
    export OMP_NUM_THREADS="${THREADS}" TMPDIR="${D}/tmp" PYSCF_TMPDIR="${D}/tmp"
    python3 "${DITY}" qm --task "${T}" --outdir "${D}" --memory-mb "${MEM_MB}" > "${D}/qm.log" 2>&1
    rm -rf "${D}/tmp"
  ) &
  pids+=($!)
done

fail=0
for i in "${!pids[@]}"; do
  if wait "${pids[$i]}"; then
    tail -n 1 "${WORK}/${TASKS[$i]}/qm.log"
  else
    echo "ERROR: ${TASKS[$i]} failed, see ${WORK}/${TASKS[$i]}/qm.log"
    tail -n 20 "${WORK}/${TASKS[$i]}/qm.log"
    fail=1
  fi
done
[ "${fail}" -eq 0 ] || exit 1

echo
echo "=== Fit"
python3 "${DITY}" fit --workdir "${WORK}" --ffbase "${FFBASE}" \
  --outdir "${SSI_ROOT}/dopa/params/results" --ffdir "${SSI_ROOT}/dopa/ff"
