#!/bin/bash -l
#SBATCH -o __SSI_ROOT__/logs/%x.%j.out
#SBATCH -e __SSI_ROOT__/logs/%x.%j.err
#SBATCH -J dopa_analysis
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=60G
#SBATCH --time=08:00:00
#SBATCH --mail-type=none

# Post-process every replica (whole protein, centered; reference frame;
# side-chain SASA of the mutation site) and compare all systems.

set -euo pipefail

source __SSI_ROOT__/config/workflow.env
source "${SSI_ROOT}/dopa/dopa.env"

source "${SSI_ROOT}/bin/ssi_lib.sh"
ssi_load_modules python
ssi_setup_env

gmx1() { mpirun -np 1 "${GMX}" "$@"; }

for SYSTEM in ${DOPA_SYSTEMS}; do
  for r in $(seq 1 "${DOPA_NREP}"); do
    REP_DIR="${DOPA_RUNS}/${SYSTEM}/rep${r}"
    if [ ! -s "${REP_DIR}/md.xtc" ]; then
      echo "skip ${SYSTEM}/rep${r}: no md.xtc"
      continue
    fi
    echo "=== ${SYSTEM}/rep${r}"
    cd "${DOPA_RUNS}/${SYSTEM}"      # residuetypes.dat with DAH lives here
    echo Protein | gmx1 convert-tpr -s "rep${r}/md.tpr" -n index.ndx -o "rep${r}/md_protein.tpr"
    echo Protein | gmx1 trjconv -s "rep${r}/md.tpr" -f "rep${r}/npt.gro" -n index.ndx \
      -o "rep${r}/protein_ref.gro" -pbc mol
    printf 'Protein\nProtein\n' | gmx1 trjconv -s "rep${r}/md_protein.tpr" -f "rep${r}/md.xtc" \
      -o "rep${r}/md_center.xtc" -pbc mol -center
    gmx1 sasa -s "rep${r}/md_protein.tpr" -f "rep${r}/md_center.xtc" -surface all \
      -output "resnr ${DOPA_MUT} and not name N H CA HA C O" \
      -o "rep${r}/sasa_res${DOPA_MUT}.xvg"
  done
done

echo
echo "=== Comparison"
python3 "${SSI_ROOT}/dopa/src/analyze_dopa.py" \
  --runs "${DOPA_RUNS}" --systems ${DOPA_SYSTEMS} --nrep "${DOPA_NREP}" \
  --first "${DOPA_FIRST}" --last "${DOPA_LAST}" --mut "${DOPA_MUT}" \
  --skip-ns "${DOPA_SKIP_NS}" --out "${DOPA_RUNS}/analysis"
