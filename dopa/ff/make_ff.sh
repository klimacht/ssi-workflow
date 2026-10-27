#!/usr/bin/env bash
# make_ff.sh — create a copy of a GROMACS force field with the hydroxylated
# tyrosines (every <res>.rtp/<res>.hdb in this directory: DAH, TOP, ...).
#
# Writes into DEST_DIR:
#   <BASE_FF>-dopa.ff/   copy of BASE_FF with the residues appended
#   residuetypes.dat     GROMACS default + "<RES> Protein" for each residue
#
# GROMACS searches the working directory first, so pdb2gmx/grompp/analysis
# tools run in DEST_DIR pick up both files.
#
# Usage: make_ff.sh GMX_TOPDIR BASE_FF DEST_DIR
#   e.g. make_ff.sh /path/gromacs/share/gromacs/top amber99sb-ildn runs/dopa/wt

set -euo pipefail

if [ $# -ne 3 ]; then
  echo "Usage: $0 GMX_TOPDIR BASE_FF DEST_DIR" >&2
  exit 1
fi

TOPDIR="$1"
BASE_FF="$2"
DEST="$3"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FF="${DEST}/${BASE_FF}-dopa.ff"

RTPS=( "${HERE}"/*.rtp )
for f in "${TOPDIR}/${BASE_FF}.ff/aminoacids.rtp" "${TOPDIR}/residuetypes.dat" "${RTPS[@]}"; do
  if [ ! -s "$f" ]; then
    echo "ERROR: missing $f" >&2
    exit 1
  fi
done

mkdir -p "${DEST}"
rm -rf "${FF}"
cp -r "${TOPDIR}/${BASE_FF}.ff" "${FF}"
chmod -R u+w "${FF}"
cp "${TOPDIR}/residuetypes.dat" "${DEST}/residuetypes.dat"
chmod u+w "${DEST}/residuetypes.dat"

names=()
for rtp in "${RTPS[@]}"; do
  hdb="${rtp%.rtp}.hdb"
  [ -s "${hdb}" ] || { echo "ERROR: missing ${hdb}" >&2; exit 1; }
  res=$(awk '/^\[ *[A-Z0-9]+ *\]/ { gsub(/[][ ]/, ""); print; exit }' "${rtp}")
  cat "${rtp}" >> "${FF}/aminoacids.rtp"
  cat "${hdb}" >> "${FF}/aminoacids.hdb"
  grep -q "^${res}[[:space:]]" "${DEST}/residuetypes.dat" || printf '%s\tProtein\n' "${res}" >> "${DEST}/residuetypes.dat"
  names+=("${res}")
done
{
  echo "${BASE_FF} + hydroxylated Tyr (${names[*]}; RESP charges, see dopa/README.md)"
  tail -n +2 "${TOPDIR}/${BASE_FF}.ff/forcefield.doc"
} > "${FF}/forcefield.doc"

echo "Created ${FF} (${names[*]}) and ${DEST}/residuetypes.dat"
