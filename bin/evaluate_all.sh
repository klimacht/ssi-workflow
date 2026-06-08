#!/usr/bin/env bash
# Run full per-atom analysis for all proteins and collect a comparison table.
set -euo pipefail
ROOT="__SSI_ROOT__"
source "${ROOT}/config/workflow.env"

PROTEINS=("$@")
if [ ${#PROTEINS[@]} -eq 0 ]; then
  PROTEINS=(1DPX 1HRC 3LDJ 1FS3 1RBB 4IBA PHLP5)
fi

echo "Submitting full analysis for: ${PROTEINS[*]}"
for PDBID in "${PROTEINS[@]}"; do
  [ -d "${ROOT}/runs/${PDBID}/ssi_production" ] || { echo "  skip ${PDBID} (no production)"; continue; }
  jid=$(sbatch --parsable "${ROOT}/slurm/06_full_analysis.slurm" "${PDBID}")
  echo "${jid}" > "${ROOT}/runs/${PDBID}/full_analysis.jobid"
  echo "  ${PDBID}: job ${jid}"
done

echo
echo "When all jobs finish, build the comparison table with:"
echo "  bash ${ROOT}/bin/collect_comparison.sh ${PROTEINS[*]}"
