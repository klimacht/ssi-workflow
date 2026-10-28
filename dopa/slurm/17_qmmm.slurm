#!/bin/bash -l
#SBATCH -o __SSI_ROOT__/logs/%x.%j.out
#SBATCH -e __SSI_ROOT__/logs/%x.%j.err
#SBATCH -J dopa_qmmm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=60G
#SBATCH --time=03:00:00
#SBATCH --mail-type=none

# QM/MM single points on MD frames: ionisation and deprotonation energy of each Tyr side chain in
# the protein field, and where the hole sits in the Y65-Y127-Y112 triad (dopa/src/qmmm_tyr.py).
# Needs the output of 14_analysis.slurm (md_center.xtc, md_protein.tpr, protein_ref.gro).
#
#   sbatch 17_qmmm.slurm prep                        frames + charges (MDAnalysis), DOPA_QMMM_FRAMES per replica
#   sbatch --array=0-39 17_qmmm.slurm run sites      five Tyr of one frame per array task (PySCF from DOPA_QM_VENV)
#   sbatch --array=0-39 --time=10:00:00 17_qmmm.slurm run triad
#                                                    the triad of one frame per array task (3x the basis, ~3x longer)
#   sbatch 17_qmmm.slurm summary                     mean +- SEM per site, triad spin per ring
# dopa_submit.sh qmmm chains them (time limits DOPA_QMMM_TIME_SITES / _TRIAD). A task that hits
# its limit loses at most the state it was on: resubmit the same index and it continues.

set -euo pipefail

MODE="${1:?Usage: 17_qmmm.slurm prep|run [sites|triad|all]|summary}"
PART="${2:-all}"
source __SSI_ROOT__/dopa/bin/dopa_env.sh __SSI_ROOT__
SRC="${SSI_ROOT}/dopa/src"
OUT="${DOPA_RUNS}/qmmm/${DOPA_QMMM_SYSTEM}"
ARGS=(--qmmm "${OUT}" --sites ${DOPA_QMMM_SITES} --triad ${DOPA_QMMM_TRIAD})

module purge
module load "${MODULE_PYTHON}"

case "${MODE}" in
  prep)
    python3 "${SRC}/qmmm_prep.py" --runs "${DOPA_RUNS}" --system "${DOPA_QMMM_SYSTEM}" --nrep "${DOPA_NREP}" \
      --frames "${DOPA_QMMM_FRAMES}" --skip-ns "${DOPA_SKIP_NS}" --sites ${DOPA_QMMM_SITES} --out "${OUT}"
    python3 "${SRC}/qmmm_tyr.py" check "${ARGS[@]}"
    ;;
  run)
    source "${DOPA_QM_VENV}/bin/activate"
    WORK="${OUT}/tmp/${PART}_${SLURM_ARRAY_TASK_ID:-0}"     # sites and triad of a frame may run together
    mkdir -p "${WORK}"
    export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}" TMPDIR="${WORK}" PYSCF_TMPDIR="${WORK}"
    python3 "${SRC}/qmmm_tyr.py" run "${ARGS[@]}" --index "${SLURM_ARRAY_TASK_ID:-0}" --part "${PART}" \
      --xc "${DOPA_QMMM_XC}" --basis "${DOPA_QMMM_BASIS}" \
      --threads "${SLURM_CPUS_PER_TASK}" --memory-mb $(( 55000 ))
    rm -rf "${WORK}"
    ;;
  summary)
    python3 "${SRC}/qmmm_tyr.py" summary "${ARGS[@]}"
    ;;
  *)
    echo "unknown mode ${MODE}" >&2
    exit 1
    ;;
esac
