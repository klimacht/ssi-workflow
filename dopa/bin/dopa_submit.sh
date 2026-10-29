#!/usr/bin/env bash
# dopa_submit.sh — submit the stages of the DOPA mutant study.
#
#   dopa_submit.sh resp   [RES]             QM (all conformers, or one residue) + RESP fit
#   dopa_submit.sh setup  [SYSTEM...]       build + equilibrate (default: all systems)
#   dopa_submit.sh bench  SYSTEM            20k-step benchmark after setup
#   dopa_submit.sh prod   [SYSTEM...]       production (waits for a running setup job)
#   dopa_submit.sh all    [SYSTEM...]       setup -> production chain
#
# Production uses DOPA_PROD_NODES nodes per system (dopa.env); set
# DOPA_TIME_MIN (e.g. 01:00:00) to let Slurm shorten jobs to fit before a
# maintenance reservation.
#   dopa_submit.sh analysis                 post-processing + comparison, then the figures
#   dopa_submit.sh render                   figures and movies only (after an analysis)
#   dopa_submit.sh qmmm                     QM/MM single points per Tyr on MD frames (after an analysis)
#   dopa_submit.sh sample [SYSTEM...]       full-system frames with velocities after the production
#                                           (default DOPA_SAMPLE_SYSTEMS)
#
# Analysis and QM/MM run on DOPA_ANALYSIS_NODES / DOPA_QMMM_NODES nodes; with more than one, the
# replicas or frames run as job steps side by side (on Viper all 1-node jobs share one CPU cap).
# DOPA_AFTER=JOBID makes the first job of the stage wait for that job (afterok).
# Add --dry-run as the last argument to print the sbatch calls only.

set -euo pipefail

ROOT="__SSI_ROOT__"
source "${ROOT}/dopa/bin/dopa_env.sh" "${ROOT}"

STAGE="${1:-}"
shift || true
DRY_RUN=""
ARGS=()
for a in "$@"; do
  if [ "$a" = "--dry-run" ]; then DRY_RUN=1; else ARGS+=("$a"); fi
done
SYSTEMS=("${ARGS[@]:-}")
[ -z "${SYSTEMS[0]:-}" ] && read -r -a SYSTEMS <<< "${DOPA_SYSTEMS}"

SL="${ROOT}/dopa/slurm"
PROD=(--nodes="${DOPA_PROD_NODES:-1}")
[ -n "${DOPA_TIME_MIN:-}" ] && PROD+=(--time-min="${DOPA_TIME_MIN}")
AFTER=()
[ -n "${DOPA_AFTER:-}" ] && AFTER=(--dependency=afterok:${DOPA_AFTER})
multinode() {   # multinode N -> sbatch args for N whole nodes, one task per node holding all its CPUs
  echo "--nodes=$1 --ntasks=$1 --ntasks-per-node=1 --cpus-per-task=${DOPA_NODE_CPUS:-128} --mem=0"
}
[ -z "${DRY_RUN}" ] && mkdir -p "${DOPA_RUNS}" "${ROOT}/logs"

submit() {   # submit DESC [sbatch args...] -> prints job id
  local desc="$1"; shift
  if [ -n "${DRY_RUN}" ]; then
    echo "[DRY-RUN] ${desc}: sbatch $*" >&2
    echo 99999999
  else
    sbatch --parsable "$@"
  fi
}

case "${STAGE}" in
  resp)
    JID=$(submit qm "${SL}/10_resp_qm.slurm" "${ARGS[@]:-}")
    echo "QM:         ${JID} ${ARGS[*]:-}"
    JID2=$(submit fit --dependency=afterok:${JID} "${SL}/11_resp_fit.slurm")
    echo "RESP fit:   ${JID2} (after ${JID})"
    ;;
  setup|all)
    for s in "${SYSTEMS[@]}"; do
      [ -z "${DRY_RUN}" ] && mkdir -p "${DOPA_RUNS}/${s}"
      JID=$(submit "setup ${s}" -J "dopa_setup_${s}" "${SL}/12_setup.slurm" "${s}")
      [ -z "${DRY_RUN}" ] && echo "${JID}" > "${DOPA_RUNS}/${s}/setup.jobid"
      echo "setup ${s}: ${JID}"
      if [ "${STAGE}" = "all" ]; then
        JID2=$(submit "prod ${s}" -J "dopa_md_${s}" "${PROD[@]}" --dependency=afterok:${JID} "${SL}/13_production.slurm" "${s}")
        echo "prod  ${s}: ${JID2} (after ${JID})"
      fi
    done
    ;;
  bench)
    s="${SYSTEMS[0]}"
    JID=$(submit "bench ${s}" -J "dopa_bench_${s}" "${SL}/13_production.slurm" "${s}" bench)
    echo "bench ${s}: ${JID}"
    ;;
  prod)
    for s in "${SYSTEMS[@]}"; do
      DEP=()
      if [ -s "${DOPA_RUNS}/${s}/setup.jobid" ]; then
        SJ=$(cat "${DOPA_RUNS}/${s}/setup.jobid")
        if squeue -h -j "${SJ}" 2>/dev/null | grep -q .; then
          DEP=(--dependency=afterok:${SJ})
        fi
      fi
      JID=$(submit "prod ${s}" -J "dopa_md_${s}" "${PROD[@]}" "${DEP[@]}" "${SL}/13_production.slurm" "${s}")
      echo "prod ${s}: ${JID} ${DEP[*]}"
    done
    ;;
  analysis)
    AN=()
    if [ "${DOPA_ANALYSIS_NODES:-1}" -gt 1 ]; then
      read -r -a AN <<< "$(multinode "${DOPA_ANALYSIS_NODES}") --time=${DOPA_ANALYSIS_TIME:-02:00:00}"
    fi
    JID=$(submit analysis "${AN[@]}" "${AFTER[@]}" "${SL}/14_analysis.slurm")
    echo "analysis: ${JID} (${DOPA_ANALYSIS_NODES:-1} node(s))"
    JID2=$(submit render --dependency=afterok:${JID} "${SL}/16_render.slurm")
    echo "render:   ${JID2} (after ${JID})"
    ;;
  render)
    JID=$(submit render "${SL}/16_render.slurm")
    echo "render: ${JID}"
    ;;
  qmmm)
    N=$(( DOPA_NREP * DOPA_QMMM_FRAMES ))
    if [ "${DOPA_QMMM_NODES:-1}" -gt 1 ]; then
      read -r -a QN <<< "$(multinode "${DOPA_QMMM_NODES}") --time=${DOPA_QMMM_TIME_PACK:-06:00:00}"
      JID=$(submit "qmmm pack" -J dopa_qmmm "${QN[@]}" "${AFTER[@]}" "${SL}/17_qmmm.slurm" pack)
      echo "qmmm:         ${JID} (prep, ${N} frames x sites/triad, summary; ${DOPA_QMMM_NODES} nodes)"
      exit 0
    fi
    ARRAY="0-$(( N - 1 ))%${DOPA_QMMM_PARALLEL:-8}"
    JID=$(submit "qmmm prep" -J dopa_qmmm_prep --cpus-per-task=4 --mem=16G --time=01:00:00 "${AFTER[@]}" "${SL}/17_qmmm.slurm" prep)
    echo "qmmm prep:    ${JID}"
    JID2=$(submit "qmmm sites" -J dopa_qmmm_sites --dependency=afterok:${JID} --array=${ARRAY} \
      --time="${DOPA_QMMM_TIME_SITES:-03:00:00}" "${SL}/17_qmmm.slurm" run sites)
    echo "qmmm sites:   ${JID2} (${N} frames, after ${JID})"
    DEP="afterany:${JID2}"
    if [ -n "${DOPA_QMMM_TRIAD// /}" ]; then
      JIDT=$(submit "qmmm triad" -J dopa_qmmm_triad --dependency=afterok:${JID} --array=${ARRAY} \
        --time="${DOPA_QMMM_TIME_TRIAD:-04:00:00}" "${SL}/17_qmmm.slurm" run triad)
      echo "qmmm triad:   ${JIDT} (${N} frames, after ${JID})"
      DEP="${DEP}:${JIDT}"
    fi
    JID3=$(submit "qmmm summary" -J dopa_qmmm_sum --dependency=${DEP} --cpus-per-task=1 --mem=4G --time=00:20:00 "${SL}/17_qmmm.slurm" summary)
    echo "qmmm summary: ${JID3} (after ${DEP#afterany:})"
    ;;
  sample)
    [ -z "${ARGS[0]:-}" ] && read -r -a SYSTEMS <<< "${DOPA_SAMPLE_SYSTEMS:-wt}"
    for s in "${SYSTEMS[@]}"; do
      JID=$(submit "sample ${s}" -J "dopa_sample_${s}" "${PROD[@]}" "${SL}/18_sample.slurm" "${s}")
      echo "sample ${s}: ${JID}"
    done
    ;;
  *)
    sed -n '2,/^$/p' "$0"
    exit 1
    ;;
esac
