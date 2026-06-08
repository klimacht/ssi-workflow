#!/usr/bin/env bash
set -euo pipefail

if [ $# -lt 1 ]; then
  echo "Usage: $0 PDBID"
  echo "Example: $0 1DPX"
  exit 1
fi

PDBID=$(echo "$1" | tr '[:lower:]' '[:upper:]')
ROOT="__SSI_ROOT__"

mkdir -p "${ROOT}/logs" "${ROOT}/runs/${PDBID}"

jid=$(sbatch --parsable "${ROOT}/slurm/00_prepare.slurm" "${PDBID}")
echo "${jid}" > "${ROOT}/runs/${PDBID}/prepare.jobid"

echo "Submitted preparation job for ${PDBID}"
echo "Job ID: ${jid}"
echo
echo "Monitor with:"
echo "  squeue -j ${jid}"
echo
echo "When finished, inspect:"
echo "  less ${ROOT}/logs/00_prepare.${jid}.out"
echo "  less ${ROOT}/logs/00_prepare.${jid}.err"
