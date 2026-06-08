#!/usr/bin/env bash
set -euo pipefail

ROOT="__SSI_ROOT__"
source "${ROOT}/config/workflow.env"

PDBID=$(echo "${1:-}" | tr '[:lower:]' '[:upper:]')
DRY_RUN="${2:-}"

if [ -z "$PDBID" ]; then
  echo "Usage: $0 PDBID [--dry-run]"
  exit 1
fi

RUN_DIR="${ROOT}/runs/${PDBID}"
mkdir -p "${RUN_DIR}" "${ROOT}/logs"

PHI_COUNT=$(echo ${PHI_LIST} | wc -w)
ARRAY_MAX=$((PHI_COUNT - 1))

do_sbatch() {
  # do_sbatch JOBID_VAR DESC [sbatch args...]
  local VARNAME="$1"; shift
  local DESC="$1"; shift
  if [ "$DRY_RUN" = "--dry-run" ]; then
    printf "    [DRY-RUN] %-14s %s\n" "$DESC" "$*"
    eval "${VARNAME}=99999999"
  else
    local JID
    JID=$(sbatch --parsable "$@")
    eval "${VARNAME}=${JID}"
  fi
}

echo "=========================================="
echo "SSI Full Workflow: ${PDBID}"
echo "$(date)"
echo "=========================================="

echo
echo "[1] Prepare"
do_sbatch JID_PREP prepare \
  "${ROOT}/slurm/00_prepare.slurm" "${PDBID}"
echo "${JID_PREP}" > "${RUN_DIR}/prepare.jobid"
echo "    Job ID: ${JID_PREP}"

echo
echo "[2] Equilibrate (after ${JID_PREP})"
do_sbatch JID_EQUIL equilibrate \
  --dependency=afterok:${JID_PREP} \
  "${ROOT}/slurm/01_equilibrate.slurm" "${PDBID}"
echo "${JID_EQUIL}" > "${RUN_DIR}/equilibrate.jobid"
echo "    Job ID: ${JID_EQUIL}"

echo
echo "[3] Unbiased production (after ${JID_EQUIL})"
do_sbatch JID_UNBIAS unbiased \
  --dependency=afterok:${JID_EQUIL} \
  "${ROOT}/slurm/02_unbiased.slurm" "${PDBID}"
echo "${JID_UNBIAS}" > "${RUN_DIR}/unbiased.jobid"
echo "    Job ID: ${JID_UNBIAS}"

echo
echo "[4] Preflight (after ${JID_UNBIAS})"
do_sbatch JID_PRE preflight \
  --dependency=afterok:${JID_UNBIAS} \
  "${ROOT}/slurm/03_ssi_preflight.slurm" "${PDBID}"
echo "${JID_PRE}" > "${RUN_DIR}/ssi_preflight.jobid"
echo "    Job ID: ${JID_PRE}"

echo
echo "[5] SSI production array 0-${ARRAY_MAX} (after ${JID_PRE})"
do_sbatch JID_PROD ssi_production \
  --dependency=afterok:${JID_PRE} \
  --array=0-${ARRAY_MAX} \
  "${ROOT}/slurm/04_ssi_production.slurm" "${PDBID}"
echo "${JID_PROD}" > "${RUN_DIR}/ssi_production.jobid"
echo "    Job ID: ${JID_PROD}"

echo
echo "[6] Analysis (after ${JID_PROD})"
do_sbatch JID_ANA analysis \
  --dependency=afterok:${JID_PROD} \
  --job-name=ssi_analysis \
  --nodes=1 --ntasks-per-node=1 --cpus-per-task=8 \
  --mem=32G --time=02:00:00 \
  --output="${ROOT}/logs/05_analysis.%j.out" \
  --error="${ROOT}/logs/05_analysis.%j.err" \
  --wrap="source ${ROOT}/config/workflow.env && module load ${MODULE_PYTHON} && bash ${ROOT}/bin/run_analysis.sh ${PDBID}"
echo "${JID_ANA}" > "${RUN_DIR}/analysis.jobid"
echo "    Job ID: ${JID_ANA}"

echo
echo "=========================================="
echo "Job chain for ${PDBID}:"
printf "  %-14s %s\n" "prepare"     "${JID_PREP}"
printf "  %-14s %s\n" "equilibrate" "${JID_EQUIL}"
printf "  %-14s %s\n" "unbiased"    "${JID_UNBIAS}"
printf "  %-14s %s\n" "preflight"   "${JID_PRE}"
printf "  %-14s %s  (array 0-${ARRAY_MAX})\n" "production"  "${JID_PROD}"
printf "  %-14s %s\n" "analysis"    "${JID_ANA}"
echo
echo "Monitor pipeline:"
echo "  watch -n 30 \"squeue -u \$USER -o '%.10i %.12j %.8T %.10M %R'\""
echo
echo "Monitor SSI production (once running):"
echo "  watch -n 30 bash ${ROOT}/bin/monitor_ssi_production.sh ${PDBID}"
echo "=========================================="
