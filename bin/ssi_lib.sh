#!/usr/bin/env bash
# ssi_lib.sh — shared functions for the SSI workflow.
# Source this after sourcing config/workflow.env.
#
#   source "${SSI_ROOT}/bin/ssi_lib.sh"
#
# Provides:
#   ssi_load_modules            — purge + load gcc/impi (+python if $1=python)
#   ssi_setup_env               — plumed paths, OMP, GMX_MAXBACKUP
#   ssi_require FILE...         — exit if any file missing/empty
#   ssi_require_posre DIR       — exit unless DIR has at least one posre*.itp
#   ssi_stage_topology SRC DST  — copy topol.top + all posre*.itp + topol_*.itp
#
# All functions assume `set -euo pipefail` is active in the caller.

ssi_load_modules() {
  module purge
  module load "${MODULE_GCC}"
  module load "${MODULE_IMPI}"
  if [ "${1:-}" = "python" ]; then
    module load "${MODULE_PYTHON}"
  fi
}

ssi_setup_env() {
  export LD_LIBRARY_PATH="${PLUMED_LIB}:${LD_LIBRARY_PATH:-}"
  export PLUMED_KERNEL="${PLUMED_KERNEL}"
  export GMX_MAXBACKUP=-1
  export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
}

# ssi_require FILE [FILE...] — every file must exist and be non-empty.
ssi_require() {
  local missing=0
  local f
  for f in "$@"; do
    if [ ! -s "$f" ]; then
      echo "ERROR: missing or empty required file: $f" >&2
      missing=1
    else
      ls -lh "$f"
    fi
  done
  [ "$missing" -eq 0 ] || exit 1
}

# ssi_require_posre DIR — at least one posre*.itp must exist in DIR.
# Prints the list; exits on failure.
ssi_require_posre() {
  local dir="$1"
  local files=( $(ls "${dir}"/posre*.itp 2>/dev/null) )
  if [ ${#files[@]} -eq 0 ]; then
    echo "ERROR: no posre*.itp files found in ${dir}" >&2
    exit 1
  fi
  echo "Found ${#files[@]} posre file(s) in ${dir}:"
  ls -lh "${files[@]}"
}

# ssi_stage_topology SRC_DIR DST_DIR
# Copies topol.top and all chain include files (posre*.itp, topol_*.itp)
# from SRC_DIR into the current directory (DST_DIR is informational).
# Removes any stale copies first.
ssi_stage_topology() {
  local src="$1"
  rm -f topol.top posre*.itp topol_*.itp
  cp -p "${src}/topol.top" topol.top
  local f
  for f in "${src}"/posre*.itp "${src}"/topol_*.itp; do
    [ -f "$f" ] && cp -p "$f" "$(basename "$f")"
  done
  if [ ! -s topol.top ]; then
    echo "ERROR: local topol.top is empty after copy from ${src}" >&2
    exit 1
  fi
}

# ssi_phidir PHI — echo the directory name for a phi value (e.g. 6 -> phi_06p0)
ssi_phidir() {
  local phi="$1"
  local w="${phi%.*}"
  local frac="${phi#*.}"
  [ "$frac" = "$phi" ] && frac=0
  printf 'phi_%02dp%d' "$w" "${frac:0:1}"
}
