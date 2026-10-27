#!/usr/bin/env bash
# make_qm_venv.sh — Python venv with PySCF + geomeTRIC for the RESP QM step.
#
# Built on top of the cluster Python module (numpy/scipy from there,
# --system-site-packages). Pinned versions; no pip cache and no /tmp use.
#
# Usage: make_qm_venv.sh VENV_DIR      (e.g. /ptmp/$USER/venvs/qm)

set -euo pipefail

VENV="${1:?Usage: $0 VENV_DIR}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
source "${ROOT}/config/workflow.env"

type module >/dev/null 2>&1 || source /etc/profile.d/modules.sh
module purge
module load "${MODULE_PYTHON}"

if [ -e "${VENV}" ]; then
  echo "ERROR: ${VENV} already exists" >&2
  exit 1
fi
mkdir -p "${VENV}"
python3 -m venv --system-site-packages "${VENV}"
mkdir -p "${VENV}/pip-tmp"
TMPDIR="${VENV}/pip-tmp" "${VENV}/bin/python" -m pip install --no-cache-dir \
  pyscf==2.14.0 geometric==1.1.1
rm -rf "${VENV}/pip-tmp"
"${VENV}/bin/python" -c "import pyscf, geometric; print('pyscf', pyscf.__version__, 'geometric', geometric.__version__)"
du -sh "${VENV}"
