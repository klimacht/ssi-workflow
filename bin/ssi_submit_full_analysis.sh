#!/usr/bin/env bash
# Submit the full per-atom dewetting analysis as a Slurm job.
#   bash bin/ssi_submit_full_analysis.sh PHLP5
set -euo pipefail
ROOT="__SSI_ROOT__"
source "${ROOT}/config/workflow.env"
PDBID=$(echo "${1:-}" | tr '[:lower:]' '[:upper:]')
if [ -z "$PDBID" ]; then echo "Usage: $0 PDBID"; exit 1; fi
jid=$(sbatch --parsable "${ROOT}/slurm/06_full_analysis.slurm" "${PDBID}")
echo "${jid}" > "${ROOT}/runs/${PDBID}/full_analysis.jobid"
echo "Submitted full analysis for ${PDBID}: job ${jid}"
echo "Log: ${ROOT}/logs/06_full_analysis.${jid}.out"
