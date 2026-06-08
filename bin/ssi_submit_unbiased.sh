#!/usr/bin/env bash
set -euo pipefail

if [ $# -lt 1 ]; then
  echo "Usage: $0 PDBID"
  echo "Example: $0 1DPX"
  exit 1
fi

PDBID=$(echo "$1" | tr '[:lower:]' '[:upper:]')
ROOT="__SSI_ROOT__"

EQUIL_DIR="${ROOT}/runs/${PDBID}/equil"

if [ ! -f "${EQUIL_DIR}/npt.gro" ]; then
  echo "ERROR: missing ${EQUIL_DIR}/npt.gro"
  echo "Run equilibration first:"
  echo "  bash ${ROOT}/bin/ssi_submit_equilibrate.sh ${PDBID}"
  exit 1
fi

mkdir -p "${ROOT}/logs" "${ROOT}/runs/${PDBID}"

jid=$(sbatch --parsable "${ROOT}/slurm/02_unbiased.slurm" "${PDBID}")
echo "${jid}" > "${ROOT}/runs/${PDBID}/unbiased.jobid"

echo "Submitted unbiased production job for ${PDBID}"
echo "Job ID: ${jid}"
echo
echo "Monitor with:"
echo "  squeue -j ${jid}"
echo
echo "When finished, inspect:"
echo "  less ${ROOT}/logs/02_unbiased.${jid}.out"
echo "  less ${ROOT}/logs/02_unbiased.${jid}.err"
