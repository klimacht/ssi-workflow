#!/usr/bin/env bash
# make_ff.sh — create a DOPA-enabled copy of a GROMACS force field.
#
# Writes into DEST_DIR:
#   <BASE_FF>-dopa.ff/   copy of BASE_FF with residue DAH (dah.rtp, dah.hdb) appended
#   residuetypes.dat     GROMACS default + "DAH Protein"
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

for f in "${TOPDIR}/${BASE_FF}.ff/aminoacids.rtp" "${TOPDIR}/residuetypes.dat" \
         "${HERE}/dah.rtp" "${HERE}/dah.hdb"; do
  if [ ! -s "$f" ]; then
    echo "ERROR: missing $f" >&2
    exit 1
  fi
done

mkdir -p "${DEST}"
rm -rf "${FF}"
cp -r "${TOPDIR}/${BASE_FF}.ff" "${FF}"
chmod -R u+w "${FF}"

cat "${HERE}/dah.rtp" >> "${FF}/aminoacids.rtp"
cat "${HERE}/dah.hdb" >> "${FF}/aminoacids.hdb"
{
  echo "${BASE_FF} + L-DOPA (DAH, RESP charges, see dopa/README.md)"
  tail -n +2 "${TOPDIR}/${BASE_FF}.ff/forcefield.doc"
} > "${FF}/forcefield.doc"

cp "${TOPDIR}/residuetypes.dat" "${DEST}/residuetypes.dat"
chmod u+w "${DEST}/residuetypes.dat"
grep -q '^DAH[[:space:]]' "${DEST}/residuetypes.dat" || printf 'DAH\tProtein\n' >> "${DEST}/residuetypes.dat"

echo "Created ${FF} and ${DEST}/residuetypes.dat"
