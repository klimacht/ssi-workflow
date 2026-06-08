#!/usr/bin/env bash
set -euo pipefail

if [ $# -lt 1 ]; then
  echo "Usage: $0 PDBID [PHI]"
  echo "Example: $0 1DPX 0"
  exit 1
fi

PDBID=$(echo "$1" | tr '[:lower:]' '[:upper:]')
PHI="${2:-0}"
ROOT="__SSI_ROOT__"

if [ ! -s "${ROOT}/runs/${PDBID}/unbiased/prod_unbiased.gro" ]; then
  echo "ERROR: missing unbiased production output."
  echo "Run first:"
  echo "  bash ${ROOT}/bin/ssi_submit_unbiased.sh ${PDBID}"
  exit 1
fi

mkdir -p "${ROOT}/logs" "${ROOT}/runs/${PDBID}"

jid=$(sbatch --parsable "${ROOT}/slurm/03_ssi_preflight.slurm" "${PDBID}" "${PHI}")
echo "${jid}" > "${ROOT}/runs/${PDBID}/ssi_preflight.jobid"

echo "Submitted SSI preflight job for ${PDBID}, phi=${PHI}"
echo "Job ID: ${jid}"
echo
echo "Monitor with:"
echo "  squeue -j ${jid}"
echo
echo "When finished, inspect:"
echo "  less ${ROOT}/logs/03_ssi_preflight.${jid}.out"
echo "  less ${ROOT}/logs/03_ssi_preflight.${jid}.err"
