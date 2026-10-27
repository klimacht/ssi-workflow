#!/bin/bash -l
#SBATCH -o __SSI_ROOT__/logs/%x.%j.out
#SBATCH -e __SSI_ROOT__/logs/%x.%j.err
#SBATCH -J dopa_md
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=2
#SBATCH --mem=120G
#SBATCH --time=12:00:00
#SBATCH --mail-type=none

# Unrestrained production for all replicas of one system (-multidir).
# Continues from checkpoints and resubmits itself until every replica has
# written md.gro (at most DOPA_MAX_RESUBMIT times).
#
#   sbatch 13_production.slurm SYSTEM          production
#   sbatch 13_production.slurm SYSTEM bench    20k-step benchmark (bench.log)

set -euo pipefail

if [ $# -lt 1 ]; then
  echo "Usage: sbatch 13_production.slurm SYSTEM [bench]"
  exit 1
fi
SYSTEM="$1"
MODE="${2:-prod}"

source __SSI_ROOT__/config/workflow.env
source "${SSI_ROOT}/dopa/dopa.env"

source "${SSI_ROOT}/bin/ssi_lib.sh"
ssi_load_modules
ssi_setup_env

SYS_DIR="${DOPA_RUNS}/${SYSTEM}"
MAX_RESUBMIT="${DOPA_MAX_RESUBMIT:-20}"
REPS=()
for r in $(seq 1 "${DOPA_NREP}"); do REPS+=("rep${r}"); done

cd "${SYS_DIR}"
for rep in "${REPS[@]}"; do ssi_require "${rep}/md.tpr"; done

echo "============================================================"
echo "DOPA production job (${MODE})"
echo "SYSTEM                = ${SYSTEM}"
echo "SYS_DIR               = ${SYS_DIR}"
echo "REPLICAS              = ${REPS[*]}"
echo "SLURM_NTASKS          = ${SLURM_NTASKS}"
echo "SLURM_CPUS_PER_TASK   = ${SLURM_CPUS_PER_TASK}"
echo "============================================================"

if [ "${MODE}" = "bench" ]; then
  mpirun -np "${SLURM_NTASKS}" "${GMX}" mdrun -multidir "${REPS[@]}" \
    -s md.tpr -deffnm bench -nsteps 20000 -resethway -noconfout \
    -ntomp "${SLURM_CPUS_PER_TASK}"
  grep -H "Performance:" rep*/bench.log
  exit 0
fi

mpirun -np "${SLURM_NTASKS}" "${GMX}" mdrun -multidir "${REPS[@]}" \
  -s md.tpr -deffnm md -cpi md.cpt -maxh 11.5 \
  -ntomp "${SLURM_CPUS_PER_TASK}"

n_done=0
for rep in "${REPS[@]}"; do
  [ -s "${rep}/md.gro" ] && n_done=$((n_done + 1))
  grep -H "Performance:" "${rep}/md.log" | tail -n 1 || true
done

if [ "${n_done}" -eq "${#REPS[@]}" ]; then
  echo "All ${n_done} replicas finished for ${SYSTEM}."
  exit 0
fi

n_sub=$(cat resubmit.count 2>/dev/null || echo 0)
if [ "${n_sub}" -ge "${MAX_RESUBMIT}" ]; then
  echo "ERROR: ${n_sub} resubmissions reached, not resubmitting (${n_done}/${#REPS[@]} done)"
  exit 1
fi
echo $((n_sub + 1)) > resubmit.count
echo "${n_done}/${#REPS[@]} replicas finished; resubmitting (#$((n_sub + 1)))"
sbatch -J "dopa_md_${SYSTEM}" "${SSI_ROOT}/dopa/slurm/13_production.slurm" "${SYSTEM}"
