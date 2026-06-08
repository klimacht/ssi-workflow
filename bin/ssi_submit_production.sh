#!/usr/bin/env bash
set -euo pipefail

if [ $# -lt 1 ]; then
  echo "Usage: $0 PDBID [PHI_INDICES]"
  echo "  PHI_INDICES: optional sbatch array spec, default 0-12"
  echo "Example: $0 1DPX"
  echo "Example: $0 1DPX 0-6"
  exit 1
fi

PDBID=$(echo "$1" | tr '[:lower:]' '[:upper:]')
ARRAY_SPEC="${2:-0-12}"
ROOT="__SSI_ROOT__"

source "${ROOT}/config/workflow.env"

if [ ! -f "${SSI_ROOT}/runs/${PDBID}/unbiased/prod_unbiased.gro" ]; then
  echo "ERROR: missing prod_unbiased.gro — run unbiased production first."
  exit 1
fi

mkdir -p "${ROOT}/logs" "${ROOT}/runs/${PDBID}"

jid=$(sbatch --parsable \
  --array="${ARRAY_SPEC}" \
  "${ROOT}/slurm/04_ssi_production.slurm" \
  "${PDBID}")

echo "${jid}" > "${ROOT}/runs/${PDBID}/ssi_production.jobid"

echo "Submitted SSI production array job for ${PDBID}"
echo "Job ID:     ${jid}"
echo "Array spec: ${ARRAY_SPEC}"
echo "PHI_LIST:   ${PHI_LIST}"
echo
echo "Monitor:    squeue -j ${jid}"
echo "Logs:       ${ROOT}/logs/04_ssi_production.${jid}_*.out"
