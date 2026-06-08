#!/usr/bin/env bash
# ssi_rerun_chain.sh — rerun the full pipeline for a multimeric crystal
# structure, keeping only one chain (the paper's monomer).
#
# Usage:
#   bash bin/ssi_rerun_chain.sh 3LDJ A
#   bash bin/ssi_rerun_chain.sh 1RBB A
#
# Wipes previous stage outputs for that PDBID, then resubmits the whole
# chain with Slurm dependencies. The PDB is re-downloaded (or reused if a
# local copy is present) and cleaned to the single requested chain.

set -euo pipefail

ROOT="__SSI_ROOT__"
source "${ROOT}/config/workflow.env"

PDBID=$(echo "${1:-}" | tr '[:lower:]' '[:upper:]')
CHAIN="${2:-}"

if [ -z "$PDBID" ] || [ -z "$CHAIN" ]; then
  echo "Usage: $0 PDBID CHAIN"
  echo "Example: $0 3LDJ A"
  exit 1
fi

RUN_DIR="${ROOT}/runs/${PDBID}"

echo "Wiping previous stage outputs for ${PDBID} (keeping pdb/)..."
rm -rf "${RUN_DIR}/prep" "${RUN_DIR}/equil" "${RUN_DIR}/unbiased" \
       "${RUN_DIR}/ssi_preflight" "${RUN_DIR}/ssi_production" \
       "${RUN_DIR}/analysis" "${RUN_DIR}/dewetting"
mkdir -p "${RUN_DIR}"

PHI_COUNT=$(echo ${PHI_LIST} | wc -w)
ARRAY_MAX=$((PHI_COUNT - 1))

echo "Submitting full pipeline for ${PDBID}, chain ${CHAIN}..."

JID_PREP=$(sbatch --parsable "${ROOT}/slurm/00_prepare.slurm" "${PDBID}" "${CHAIN}")
echo "${JID_PREP}" > "${RUN_DIR}/prepare.jobid"

JID_EQUIL=$(sbatch --parsable --dependency=afterok:${JID_PREP} \
  "${ROOT}/slurm/01_equilibrate.slurm" "${PDBID}")
echo "${JID_EQUIL}" > "${RUN_DIR}/equilibrate.jobid"

JID_UNBIAS=$(sbatch --parsable --dependency=afterok:${JID_EQUIL} \
  "${ROOT}/slurm/02_unbiased.slurm" "${PDBID}")
echo "${JID_UNBIAS}" > "${RUN_DIR}/unbiased.jobid"

JID_PRE=$(sbatch --parsable --dependency=afterok:${JID_UNBIAS} \
  "${ROOT}/slurm/03_ssi_preflight.slurm" "${PDBID}")
echo "${JID_PRE}" > "${RUN_DIR}/ssi_preflight.jobid"

JID_PROD=$(sbatch --parsable --dependency=afterok:${JID_PRE} \
  --array=0-${ARRAY_MAX} \
  "${ROOT}/slurm/04_ssi_production.slurm" "${PDBID}")
echo "${JID_PROD}" > "${RUN_DIR}/ssi_production.jobid"

JID_ANA=$(sbatch --parsable --dependency=afterok:${JID_PROD} \
  --job-name=ssi_analysis \
  --nodes=1 --ntasks-per-node=1 --cpus-per-task=8 \
  --mem=32G --time=02:00:00 \
  --output="${ROOT}/logs/05_analysis.%j.out" \
  --error="${ROOT}/logs/05_analysis.%j.err" \
  --wrap="source ${ROOT}/config/workflow.env && module load ${MODULE_PYTHON} && bash ${ROOT}/bin/run_analysis.sh ${PDBID}")
echo "${JID_ANA}" > "${RUN_DIR}/analysis.jobid"

echo "Submitted ${PDBID} (chain ${CHAIN}):"
echo "  prepare=${JID_PREP} equil=${JID_EQUIL} unbias=${JID_UNBIAS}"
echo "  preflight=${JID_PRE} prod=${JID_PROD} analysis=${JID_ANA}"
