#!/usr/bin/env bash
set -euo pipefail

if [ -z "${SSI_ROOT:-}" ]; then
  # shellcheck disable=SC1091
  source __SSI_ROOT__/config/workflow.env
fi

module purge
module load "${MODULE_GCC}"
module load "${MODULE_IMPI}"

export LD_LIBRARY_PATH="${PLUMED_LIB}:${LD_LIBRARY_PATH:-}"
export PLUMED_KERNEL="${PLUMED_KERNEL}"
export GMX_MAXBACKUP=-1
