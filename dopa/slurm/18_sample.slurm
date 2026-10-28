#!/bin/bash -l
#SBATCH -o __SSI_ROOT__/logs/%x.%j.out
#SBATCH -e __SSI_ROOT__/logs/%x.%j.err
#SBATCH -J dopa_sample
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=2
#SBATCH --mem=120G
#SBATCH --time=04:00:00
#SBATCH --mail-type=none

# Full-system frames after the production: every replica of one system continues from its final
# state (md.gro, md.cpt) for DOPA_SAMPLE_NS ns and writes all atoms with velocities every
# DOPA_SAMPLE_DT_PS ps (rep*/sample/sample.trr). Starting points for D-NEMD (modification
# inserted into a frame) and for QM/MM with waters and ions (17_qmmm.slurm with
# DOPA_QMMM_TOP/_TPR/_XTC pointing at sample/). The production itself is not touched.
#
#   sbatch [--nodes=4] 18_sample.slurm SYSTEM
#
# Per replica afterwards: sample/sample.trr (raw, x + v), sample/sample_center.xtc (all atoms,
# whole, centred on the protein, compact unit cell) and sample/sample_ref.gro (first frame).
# Restartable: a finished tpr is kept, mdrun continues from sample.cpt.

set -euo pipefail

SYSTEM="${1:?Usage: sbatch 18_sample.slurm SYSTEM}"
source __SSI_ROOT__/dopa/bin/dopa_env.sh __SSI_ROOT__

source "${SSI_ROOT}/bin/ssi_lib.sh"
ssi_load_modules
ssi_setup_env

gmx1() { mpirun -np 1 "${GMX}" "$@"; }

NS="${DOPA_SAMPLE_NS:-10}"
DT="${DOPA_SAMPLE_DT_PS:-100}"
cd "${DOPA_RUNS}/${SYSTEM}"         # topol.top, index.ndx, residuetypes.dat
DIRS=()
for r in $(seq 1 "${DOPA_NREP}"); do
  rep="rep${r}"
  ssi_require "${rep}/md.gro" "${rep}/md.cpt"       # md.gro: production finished
  mkdir -p "${rep}/sample"
  if [ ! -s "${rep}/sample/sample.tpr" ]; then
    T0=$(head -n 1 "${rep}/md.gro" | sed -n 's/.*t= *\([0-9.]*\).*/\1/p')
    sed -e "s/^nsteps .*/nsteps                  = $(( NS * 500000 ))/" \
        -e "s/^nstxout .*/nstxout                 = $(( DT * 500 ))/" \
        -e "s/^nstvout .*/nstvout                 = $(( DT * 500 ))/" \
        "${SSI_ROOT}/dopa/mdp/md.mdp" > "${rep}/sample/sample.mdp"
    echo "tinit                   = ${T0:-0}" >> "${rep}/sample/sample.mdp"
    gmx1 grompp -f "${rep}/sample/sample.mdp" -c "${rep}/md.gro" -t "${rep}/md.cpt" \
      -p topol.top -n index.ndx -o "${rep}/sample/sample.tpr" -po "${rep}/sample/mdout.mdp"
  fi
  DIRS+=("${rep}/sample")
done

echo "============================================================"
echo "DOPA full-system sampling: ${SYSTEM}, ${NS} ns per replica, frames every ${DT} ps"
echo "DIRS                  = ${DIRS[*]}"
echo "SLURM_NTASKS          = ${SLURM_NTASKS}"
echo "============================================================"

mpirun -np "${SLURM_NTASKS}" "${GMX}" mdrun -multidir "${DIRS[@]}" \
  -s sample.tpr -deffnm sample -cpi sample.cpt -ntomp "${SLURM_CPUS_PER_TASK}"

for d in "${DIRS[@]}"; do
  ssi_require "${d}/sample.gro"
  grep -H "Performance:" "${d}/sample.log" | tail -n 1 || true
  printf 'Protein\nSystem\n' | gmx1 trjconv -s "${d}/sample.tpr" -f "${d}/sample.trr" -n index.ndx \
    -o "${d}/sample_center.xtc" -pbc mol -center -ur compact
  printf 'Protein\nSystem\n' | gmx1 trjconv -s "${d}/sample.tpr" -f "${d}/sample.trr" -n index.ndx \
    -o "${d}/sample_ref.gro" -pbc mol -center -ur compact -dump 0
done
echo "Full-system frames of ${SYSTEM}: ${DIRS[*]}"
