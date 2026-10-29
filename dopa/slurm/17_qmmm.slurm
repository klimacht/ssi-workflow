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
#   sbatch --array=0-39 --time=04:00:00 17_qmmm.slurm run triad
#                                                    the triad of one frame per array task (~1 h on a workstation)
#   sbatch 17_qmmm.slurm summary                     mean +- SEM per site, triad spin per ring
#   sbatch --nodes=N ... 17_qmmm.slurm pack          all of it in one job: prep, then every frame x part as
#                                                    its own job step (DOPA_QMMM_STEPS_PER_NODE per node), summary
# dopa_submit.sh qmmm chains them (time limits DOPA_QMMM_TIME_SITES / _TRIAD), or submits one pack job
# when DOPA_QMMM_NODES > 1. A task that hits its limit loses at most the state it was on: resubmit the
# same index (or the pack job) and it continues.

set -euo pipefail

MODE="${1:?Usage: 17_qmmm.slurm prep|run [sites|triad|all]|summary|pack}"
PART="${2:-all}"
source __SSI_ROOT__/dopa/bin/dopa_env.sh __SSI_ROOT__
SRC="${SSI_ROOT}/dopa/src"
OUT="${DOPA_RUNS}/qmmm/${DOPA_QMMM_OUT:-${DOPA_QMMM_SYSTEM}}"     # e.g. wt_water for full-system frames
ARGS=(--qmmm "${OUT}" --sites ${DOPA_QMMM_SITES} --triad ${DOPA_QMMM_TRIAD})
INDEX="${QMMM_INDEX:-${SLURM_ARRAY_TASK_ID:-0}}"                 # QMMM_INDEX: set by the pack mode

module purge
module load "${MODULE_PYTHON}"

case "${MODE}" in
  prep)
    python3 "${SRC}/qmmm_prep.py" --runs "${DOPA_RUNS}" --system "${DOPA_QMMM_SYSTEM}" --nrep "${DOPA_NREP}" \
      --frames "${DOPA_QMMM_FRAMES}" --skip-ns "${DOPA_SKIP_NS}" --sites ${DOPA_QMMM_SITES} --out "${OUT}" \
      --top "${DOPA_QMMM_TOP:-protein_ref.gro}" --tpr "${DOPA_QMMM_TPR:-md_protein.tpr}" \
      --xtc "${DOPA_QMMM_XTC:-md_center.xtc}"
    python3 "${SRC}/qmmm_tyr.py" check "${ARGS[@]}"
    ;;
  run)
    source "${DOPA_QM_VENV}/bin/activate"
    WORK="${OUT}/tmp/${PART}_${INDEX}"     # sites and triad of a frame may run together
    mkdir -p "${WORK}"
    export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}" TMPDIR="${WORK}" PYSCF_TMPDIR="${WORK}"
    python3 "${SRC}/qmmm_tyr.py" run "${ARGS[@]}" --index "${INDEX}" --part "${PART}" \
      --xc "${DOPA_QMMM_XC}" --basis "${DOPA_QMMM_BASIS}" \
      --threads "${SLURM_CPUS_PER_TASK}" --memory-mb $(( 55000 ))
    rm -rf "${WORK}"
    ;;
  summary)
    python3 "${SRC}/qmmm_tyr.py" summary "${ARGS[@]}"
    ;;
  pack)
    SELF="${SSI_ROOT}/dopa/slurm/17_qmmm.slurm"
    [ -s "${OUT}/frames.tsv" ] || bash "${SELF}" prep
    N=$(( $(wc -l < "${OUT}/frames.tsv") - 1 ))
    PER_NODE="${DOPA_QMMM_STEPS_PER_NODE:-8}"
    STEP_CPUS=$(( SLURM_CPUS_PER_TASK / PER_NODE ))         # one task per node holds the whole node
    SLOTS=$(( SLURM_NNODES * PER_NODE ))
    PARTS="sites"
    [ -n "${DOPA_QMMM_TRIAD// /}" ] && PARTS="triad sites"          # the long ones first
    mkdir -p "${OUT}/logs"
    : > "${OUT}/logs/failed.txt"
    echo "QM/MM pack: ${N} frames x (${PARTS}) as job steps, ${SLOTS} at a time, ${STEP_CPUS} CPUs each"
    for p in ${PARTS}; do
      for i in $(seq 0 $(( N - 1 ))); do
        while [ "$(jobs -rp | wc -l)" -ge "${SLOTS}" ]; do wait -n || true; done
        (
          QMMM_INDEX="${i}" srun --exact --nodes=1 --ntasks=1 --cpus-per-task="${STEP_CPUS}" --mem=55G \
            bash "${SELF}" run "${p}" > "${OUT}/logs/${p}_${i}.log" 2>&1 \
            || echo "${p} ${i}" >> "${OUT}/logs/failed.txt"
        ) &
      done
    done
    wait
    bash "${SELF}" summary || echo "summary failed" >&2
    if [ -s "${OUT}/logs/failed.txt" ]; then
      echo "ERROR: failed steps (resubmit the pack job to continue them):" >&2
      cat "${OUT}/logs/failed.txt" >&2
      exit 1
    fi
    ;;
  *)
    echo "unknown mode ${MODE}" >&2
    exit 1
    ;;
esac
