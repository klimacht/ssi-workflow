#!/usr/bin/env bash
set -euo pipefail

source __SSI_ROOT__/config/workflow.env

echo "============================================================"
echo "SSI workflow doctor"
echo "============================================================"
echo "SSI_ROOT              = ${SSI_ROOT}"
echo "GMX                   = ${GMX}"
echo "PLUMED_LIB            = ${PLUMED_LIB}"
echo "PLUMED_KERNEL         = ${PLUMED_KERNEL}"
echo "MODULE_GCC            = ${MODULE_GCC}"
echo "MODULE_IMPI           = ${MODULE_IMPI}"
echo "MODULE_PYTHON         = ${MODULE_PYTHON}"
echo "SSI_RADIUS_NM         = ${SSI_RADIUS_NM}"
echo "SSI_UNION_MODE        = ${SSI_UNION_MODE}"
echo "SSI_BIAS_MODE         = ${SSI_BIAS_MODE}"
echo "PHI_LIST              = ${PHI_LIST}"
echo "SLURM layout          = ${SLURM_NODES} node, ${SLURM_NTASKS_PER_NODE} MPI ranks/node, ${SLURM_CPUS_PER_TASK} threads/rank"
echo "============================================================"

echo
echo "[1] Checking module command..."
if command -v module >/dev/null 2>&1; then
  echo "OK: module command found"
else
  echo "ERROR: module command not found"
  exit 1
fi

echo
echo "[2] Loading compiler/MPI modules..."
module purge
module load "${MODULE_GCC}"
module load "${MODULE_IMPI}"
echo "OK: loaded ${MODULE_GCC} and ${MODULE_IMPI}"

echo
echo "[3] Setting runtime exports exactly as in production jobs..."
export LD_LIBRARY_PATH=__SOFTWARE_ROOT__/plumed/lib:${LD_LIBRARY_PATH:-}
export PLUMED_KERNEL=__SOFTWARE_ROOT__/plumed/lib/libplumedKernel.so
export GMX_MAXBACKUP=-1
export OMP_NUM_THREADS=1

echo "LD_LIBRARY_PATH=${LD_LIBRARY_PATH}"
echo "PLUMED_KERNEL=${PLUMED_KERNEL}"
echo "GMX_MAXBACKUP=${GMX_MAXBACKUP}"
echo "OMP_NUM_THREADS=${OMP_NUM_THREADS}"

echo
echo "[4] Checking GROMACS executable without launching MPI..."
if [ -x "${GMX}" ]; then
  echo "OK: ${GMX} exists and is executable"
else
  echo "ERROR: ${GMX} missing or not executable"
  exit 1
fi

echo
echo "ldd check for libplumedKernel:"
if ldd "${GMX}" | grep -i plumed; then
  echo "OK: GROMACS binary links to PLUMED and linker can see it."
else
  echo "WARNING: ldd did not show PLUMED. This may be fine if GROMACS loads PLUMED dynamically."
fi

echo
echo "[5] Checking PLUMED library..."
if [ -f "${PLUMED_KERNEL}" ]; then
  echo "OK: ${PLUMED_KERNEL} exists"
else
  echo "ERROR: ${PLUMED_KERNEL} missing"
  exit 1
fi

echo
echo "[6] Checking Python module..."
module load "${MODULE_PYTHON}"

python3 - <<'PY'
import sys
print("Python:", sys.version)
mods = ["numpy", "MDAnalysis"]
for m in mods:
    mod = __import__(m)
    print(f"OK: {m} {getattr(mod, '__version__', '')}")
PY

echo
echo "[7] Checking network access to RCSB PDB..."
tmp="$(mktemp)"
if wget -q -O "$tmp" https://files.rcsb.org/download/1DPX.pdb; then
  echo "OK: downloaded test PDB 1DPX"
  head -n 3 "$tmp"
else
  echo "WARNING: could not download from RCSB using wget"
fi
rm -f "$tmp"

echo
echo "============================================================"
echo "Doctor finished."
echo "Note: MPI/GROMACS execution is intentionally tested only inside Slurm."
echo "============================================================"
