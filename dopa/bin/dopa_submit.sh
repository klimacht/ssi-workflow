#!/usr/bin/env bash
# dopa_submit.sh — submit the stages of the DOPA mutant study.
#
#   dopa_submit.sh resp                     QM array (6 conformers) + RESP fit
#   dopa_submit.sh setup  [SYSTEM...]       build + equilibrate (default: all systems)
#   dopa_submit.sh bench  SYSTEM            20k-step benchmark after setup
#   dopa_submit.sh prod   [SYSTEM...]       production (waits for a running setup job)
#   dopa_submit.sh all    [SYSTEM...]       setup -> production chain
#   dopa_submit.sh analysis                 post-processing + comparison
#
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
    JID=$(submit qm "${SL}/10_resp_qm.slurm")
    echo "QM array:   ${JID}"
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
        JID2=$(submit "prod ${s}" -J "dopa_md_${s}" --dependency=afterok:${JID} "${SL}/13_production.slurm" "${s}")
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
      JID=$(submit "prod ${s}" -J "dopa_md_${s}" "${DEP[@]}" "${SL}/13_production.slurm" "${s}")
      echo "prod ${s}: ${JID} ${DEP[*]}"
    done
    ;;
  analysis)
    JID=$(submit analysis "${SL}/14_analysis.slurm")
    echo "analysis: ${JID}"
    ;;
  *)
    sed -n '2,12p' "$0"
    exit 1
    ;;
esac
