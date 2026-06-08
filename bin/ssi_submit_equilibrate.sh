#!/usr/bin/env bash
set -euo pipefail

if [ $# -lt 1 ]; then
  echo "Usage: $0 PDBID"
  echo "Example: $0 1DPX"
  exit 1
fi

PDBID=$(echo "$1" | tr '[:lower:]' '[:upper:]')
ROOT="__SSI_ROOT__"

PREP_DIR="${ROOT}/runs/${PDBID}/prep"

if [ ! -f "${PREP_DIR}/em.gro" ]; then
  echo "ERROR: missing ${PREP_DIR}/em.gro"
  echo "Run preparation first:"
  echo "  bash ${ROOT}/bin/ssi_submit_prepare.sh ${PDBID}"
  exit 1
fi

mkdir -p "${ROOT}/logs" "${ROOT}/runs/${PDBID}"

jid=$(sbatch --parsable "${ROOT}/slurm/01_equilibrate.slurm" "${PDBID}")
echo "${jid}" > "${ROOT}/runs/${PDBID}/equilibrate.jobid"

echo "Submitted equilibration job for ${PDBID}"
echo "Job ID: ${jid}"
echo
echo "Monitor with:"
echo "  squeue -j ${jid}"
echo
echo "When finished, inspect:"
echo "  less ${ROOT}/logs/01_equilibrate.${jid}.out"
echo "  less ${ROOT}/logs/01_equilibrate.${jid}.err"
