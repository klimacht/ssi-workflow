#!/usr/bin/env bash
# setup.sh — configure the SSI workflow for your system after cloning.
#
# Replaces the __SSI_ROOT__ and __SOFTWARE_ROOT__ placeholders in all
# scripts with absolute paths, and creates config/workflow.env from the
# example template.
#
# Usage:
#   ./setup.sh /abs/path/to/this/checkout  /abs/path/to/software
#
# Example:
#   ./setup.sh /ptmp/myuser/ssi_workflow  /home/myuser/software

set -euo pipefail

SSI_ROOT="${1:-}"
SOFTWARE_ROOT="${2:-}"

if [ -z "$SSI_ROOT" ]; then
  echo "Usage: $0 SSI_ROOT [SOFTWARE_ROOT]"
  echo
  echo "  SSI_ROOT        absolute path where this repo lives and runs/ logs/ will be created"
  echo "  SOFTWARE_ROOT   absolute path to GROMACS/PLUMED install (optional; edit config later if omitted)"
  echo
  echo "Example: $0 /ptmp/myuser/ssi_workflow /home/myuser/software"
  exit 1
fi

SSI_ROOT="${SSI_ROOT%/}"
SOFTWARE_ROOT="${SOFTWARE_ROOT%/}"

echo "Configuring SSI workflow:"
echo "  SSI_ROOT      = ${SSI_ROOT}"
echo "  SOFTWARE_ROOT = ${SOFTWARE_ROOT:-<edit config manually>}"

# 1. Replace placeholders in all tracked scripts
echo "Replacing __SSI_ROOT__ in scripts..."
grep -rl "__SSI_ROOT__" . 2>/dev/null | while read -r f; do
  sed -i "s|__SSI_ROOT__|${SSI_ROOT}|g" "$f"
done

if [ -n "$SOFTWARE_ROOT" ]; then
  echo "Replacing __SOFTWARE_ROOT__ in scripts..."
  grep -rl "__SOFTWARE_ROOT__" . 2>/dev/null | while read -r f; do
    sed -i "s|__SOFTWARE_ROOT__|${SOFTWARE_ROOT}|g" "$f"
  done
fi

# 2. Create config/workflow.env from example if it doesn't exist
if [ ! -f config/workflow.env ]; then
  echo "Creating config/workflow.env from example..."
  sed -e "s|/path/to/ssi_workflow|${SSI_ROOT}|g" \
      -e "s|/path/to/software|${SOFTWARE_ROOT:-/path/to/software}|g" \
      config/workflow.env.example > config/workflow.env
  echo "  -> review config/workflow.env and set module names for your cluster"
else
  echo "config/workflow.env already exists — leaving it untouched."
fi

# 3. Make scripts executable
chmod +x bin/*.sh 2>/dev/null || true

# 4. Prepare PLUMED with INDUS patch (if SOFTWARE_ROOT is provided)
if [ -n "$SOFTWARE_ROOT" ]; then
  echo
  echo "=== Preparing PLUMED 2.9.0 with INDUS patch ==="
  mkdir -p "${SOFTWARE_ROOT}/src"
  cd "${SOFTWARE_ROOT}/src"
  
  if [ ! -d "plumed-2.9.0" ]; then
    echo "Downloading PLUMED 2.9.0..."
    wget -q --show-progress https://github.com/plumed/plumed2/releases/download/v2.9.0/plumed-src-2.9.0.tgz || true
    tar -xzf plumed-src-2.9.0.tgz
    
    echo "Applying INDUS union_spheres performance patch..."
    cd plumed-2.9.0
    patch -p1 < "${SSI_ROOT}/indus-patches/0001-nearest-union-mode-and-derivative-sharing.patch"
    
    echo "Configuring PLUMED prefix..."
    ./configure --prefix="${SOFTWARE_ROOT}/plumed" > configure_output.log 2>&1
    echo "  -> Configuration complete (see ${SOFTWARE_ROOT}/src/plumed-2.9.0/configure_output.log)"
  else
    echo "plumed-2.9.0 source directory already exists in ${SOFTWARE_ROOT}/src. Skipping download/patch."
  fi
  cd "${SSI_ROOT}"
fi

echo
echo "Done. Next steps:"
echo "  1. Review config/workflow.env (module names, force field, phi list)"
if [ -n "$SOFTWARE_ROOT" ]; then
  echo "  2. Load your compiler/MPI modules, then build PLUMED:"
  echo "       cd ${SOFTWARE_ROOT}/src/plumed-2.9.0"
  echo "       make -j 4 && make install"
else
  echo "  2. Apply INDUS patches and build PLUMED (see indus-patches/README.md)"
fi
echo "  3. Run a protein:  bash bin/ssi_run_full.sh 1DPX"
