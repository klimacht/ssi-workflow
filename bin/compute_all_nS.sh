#!/usr/bin/env bash
# compute_all_nS.sh — compute <n_i>_phi for every phi window of a protein,
# so the per-atom transition map (combine_dewetting_phistar_map.py) can use
# the full data set instead of a single window.
#
# Usage:
#   bash bin/compute_all_nS.sh PHLP5
#   bash bin/compute_all_nS.sh 1DPX
#
# Skips windows whose nS_phi_*.npy already exists.

set -euo pipefail

ROOT="__SSI_ROOT__"
source "${ROOT}/config/workflow.env"
module load "${MODULE_PYTHON}" 2>/dev/null || true

PDBID=$(echo "${1:-}" | tr '[:lower:]' '[:upper:]')
if [ -z "$PDBID" ]; then
  echo "Usage: $0 PDBID"
  exit 1
fi

ANA_DIR="${ROOT}/runs/${PDBID}/analysis"
T_START_PS=$(awk "BEGIN{printf \"%d\", ${DISCARD_NS}*1000}")

# Ensure n0 exists
if [ ! -f "${ANA_DIR}/n0_per_atom.npy" ]; then
  echo "[*] n0 missing — computing it first"
  python3 "${ROOT}/src/compute_n0.py" \
    --pdbid "${PDBID}" --root "${ROOT}" \
    --t-start "${T_START_PS}" --radius "${SSI_RADIUS_NM}"
fi

for PHI in ${PHI_LIST}; do
  # phi dir name via the same convention as the rest of the pipeline
  W="${PHI%.*}"; FRAC="${PHI#*.}"; [ "$FRAC" = "$PHI" ] && FRAC=0
  PHIDIR=$(printf 'phi_%02dp%d' "$W" "${FRAC:0:1}")

  if [ -f "${ANA_DIR}/nS_${PHIDIR}.npy" ]; then
    echo "[skip] ${PHIDIR} (already computed)"
    continue
  fi

  if [ ! -f "${ROOT}/runs/${PDBID}/ssi_production/${PHIDIR}/ssi_phi.xtc" ]; then
    echo "[warn] ${PHIDIR}: no trajectory, skipping"
    continue
  fi

  echo "[run ] ${PHIDIR}"
  python3 "${ROOT}/src/compute_nS.py" \
    --pdbid "${PDBID}" --phi "${PHIDIR}" --root "${ROOT}" \
    --t-start "${T_START_PS}" --radius "${SSI_RADIUS_NM}"
done

echo
echo "All available nS windows computed for ${PDBID}."
echo "Now build the per-atom transition map:"
echo "  python3 ${ROOT}/src/combine_dewetting_phistar_map.py --pdbid ${PDBID}"
