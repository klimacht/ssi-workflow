#!/bin/bash -l
#SBATCH -o __SSI_ROOT__/logs/%x.%j.out
#SBATCH -e __SSI_ROOT__/logs/%x.%j.err
#SBATCH -J dopa_setup
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=2
#SBATCH --mem=120G
#SBATCH --time=04:00:00
#SBATCH --mail-type=none

# Build one system (domain + caps, optional Y->DOPA), solvate, add 0.15 M NaCl,
# minimize, then NVT + NPT with restraints for all replicas (-multidir) and
# write rep*/md.tpr for 13_production.slurm.

set -euo pipefail

if [ $# -lt 1 ]; then
  echo "Usage: sbatch 12_setup.slurm SYSTEM   (wt | dah_ce1 | y112_dah_ce2 | y65_top_ce1 | ...)"
  exit 1
fi
SYSTEM="$1"

source __SSI_ROOT__/dopa/bin/dopa_env.sh __SSI_ROOT__

source "${SSI_ROOT}/bin/ssi_lib.sh"
ssi_load_modules python
ssi_setup_env

MDP="${SSI_ROOT}/dopa/mdp"
SYS_DIR="${DOPA_RUNS}/${SYSTEM}"
GMX_TOP="$(dirname "$(readlink -f "${GMX}")")/../share/gromacs/top"
REPS=()
for r in $(seq 1 "${DOPA_NREP}"); do REPS+=("rep${r}"); done

gmx1() { mpirun -np 1 "${GMX}" "$@"; }

# mutation site and residue of this system (see make_system_pdb.py)
read -r SITE NEWRES ORIENT < <(python3 "${SSI_ROOT}/dopa/src/make_system_pdb.py"   --parse "${SYSTEM}" --default-mut "${DOPA_MUT}")

if [ $((SLURM_NTASKS % DOPA_NREP)) -ne 0 ]; then
  echo "ERROR: SLURM_NTASKS=${SLURM_NTASKS} not divisible by DOPA_NREP=${DOPA_NREP}"
  exit 1
fi

mkdir -p "${SYS_DIR}" "${SSI_ROOT}/logs"
cd "${SYS_DIR}"

echo "============================================================"
echo "DOPA setup job"
echo "SYSTEM                = ${SYSTEM}"
echo "SYS_DIR               = ${SYS_DIR}"
echo "DOMAIN                = ${DOPA_FIRST}-${DOPA_LAST}"
echo "MODIFICATION          = ${NEWRES} at residue ${SITE} (ring side ${ORIENT})"
echo "FORCE FIELD / WATER   = ${DOPA_BASE_FF}-dopa / ${DOPA_WATER}"
echo "REPLICAS              = ${REPS[*]}"
echo "GMX                   = ${GMX}"
echo "SLURM_NTASKS          = ${SLURM_NTASKS}"
echo "SLURM_CPUS_PER_TASK   = ${SLURM_CPUS_PER_TASK}"
echo "============================================================"

echo
echo "[1] Force field with the hydroxylated Tyr residues"
if [ "${NEWRES}" != "-" ]; then
  RES_LC=$(echo "${NEWRES}" | tr '[:upper:]' '[:lower:]')
  ssi_require "${SSI_ROOT}/dopa/ff/${RES_LC}.rtp" "${SSI_ROOT}/dopa/ff/${RES_LC}.hdb"
fi
bash "${SSI_ROOT}/dopa/ff/make_ff.sh" "${GMX_TOP}" "${DOPA_BASE_FF}" "${SYS_DIR}"

echo
echo "[2] Domain structure"
python3 "${SSI_ROOT}/dopa/src/make_system_pdb.py" \
  --input "${DOPA_SRC_PDB}" --out start.pdb \
  --first "${DOPA_FIRST}" --last "${DOPA_LAST}" \
  --default-mut "${DOPA_MUT}" --system "${SYSTEM}"

echo
echo "[3] pdb2gmx (ACE/NME caps, no charged termini)"
# pdb2gmx recognizes ACE/NME as capping residues and adds no NH3+/COO- termini
gmx1 pdb2gmx -f start.pdb -o processed.gro -p topol.top -i posre.itp \
  -ff "${DOPA_BASE_FF}-dopa" -water "${DOPA_WATER}" -ignh > pdb2gmx.log 2>&1
grep -E "Identified residue|Total charge" pdb2gmx.log || true
N_CHARGED_TER=$(awk 'NR>2 { a=substr($0,11,5); gsub(/ /,"",a); if (a ~ /^(H1|H2|H3|OC1|OC2|OXT)$/) n++ } END{print n+0}' processed.gro)
if ! grep -q "Identified residue ACE.* as a starting terminus" pdb2gmx.log || \
   ! grep -q "Identified residue NME.* as a ending terminus" pdb2gmx.log || \
   [ "${N_CHARGED_TER}" -ne 0 ]; then
  echo "ERROR: caps not recognized or charged termini added (${N_CHARGED_TER} atoms), see ${SYS_DIR}/pdb2gmx.log"
  exit 1
fi
if [ "${NEWRES}" != "-" ]; then
  N_RES=$(awk -v r="${NEWRES}" 'NR>2 { n0=substr($0,6,5); gsub(/ /,"",n0); if (n0 == r) n++ } END{print n+0}' processed.gro)
  N_RTP=$(awk '/\[ atoms \]/{f=1; next} /\[ bonds \]/{f=0} f && NF>=4 {n++} END{print n+0}' "${SSI_ROOT}/dopa/ff/${RES_LC}.rtp")
  echo "${NEWRES} atoms in processed.gro: ${N_RES} (template: ${N_RTP})"
  if [ "${N_RES}" -ne "${N_RTP}" ]; then
    echo "ERROR: ${NEWRES} atom count does not match the template"
    exit 1
  fi
fi

echo
echo "[4] Box, water, ions"
gmx1 editconf -f processed.gro -o boxed.gro -bt dodecahedron -d "${DOPA_BOX_NM}" -c
gmx1 solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top
gmx1 grompp -f "${SSI_ROOT}/mdp/ions.mdp" -c solvated.gro -p topol.top -o ions.tpr -maxwarn 1
echo "SOL" | gmx1 genion -s ions.tpr -o solv_ions.gro -p topol.top \
  -pname NA -nname CL -neutral -conc "${DOPA_SALT_M}"
tail -n 8 topol.top

echo
echo "[5] Energy minimization"
gmx1 grompp -f "${MDP}/em.mdp" -c solv_ions.gro -p topol.top -o em.tpr
mpirun -np 16 "${GMX}" mdrun -deffnm em -ntomp "${SLURM_CPUS_PER_TASK}"
grep -E "Potential Energy|Maximum force" em.log | tail -n 2

echo
echo "[6] Index; check that the Protein group covers every protein atom"
echo q | gmx1 make_ndx -f em.gro -o index.ndx > make_ndx.log 2>&1
N_PROT_NDX=$(awk '/^\[ Protein \]/{f=1; next} /^\[/{f=0} f{n+=NF} END{print n+0}' index.ndx)
N_PROT_GRO=$(awk 'NR>2 { rn=substr($0,6,5); gsub(/ /,"",rn); if (rn=="SOL"||rn=="NA"||rn=="CL") exit; n++ } END{print n+0}' em.gro)
echo "Protein group: ${N_PROT_NDX} atoms; protein atoms in em.gro: ${N_PROT_GRO}"
if [ "${N_PROT_NDX}" -ne "${N_PROT_GRO}" ]; then
  echo "ERROR: Protein group incomplete (is DAH in residuetypes.dat?)"
  exit 1
fi

echo
echo "[7] NVT (100 ps, restrained), one seed per replica"
for i in "${!REPS[@]}"; do
  rep="${REPS[$i]}"
  mkdir -p "${rep}"
  sed "s/__SEED__/$((DOPA_SEED_BASE + i + 1))/" "${MDP}/nvt.mdp" > "${rep}/nvt.mdp"
  gmx1 grompp -f "${rep}/nvt.mdp" -c em.gro -r em.gro -p topol.top -n index.ndx \
    -o "${rep}/nvt.tpr" -po "${rep}/mdout_nvt.mdp"
done
mpirun -np "${SLURM_NTASKS}" "${GMX}" mdrun -multidir "${REPS[@]}" -deffnm nvt -ntomp "${SLURM_CPUS_PER_TASK}"

echo
echo "[8] NPT (1 ns, restrained)"
for rep in "${REPS[@]}"; do
  gmx1 grompp -f "${MDP}/npt.mdp" -c "${rep}/nvt.gro" -t "${rep}/nvt.cpt" -r em.gro \
    -p topol.top -n index.ndx -o "${rep}/npt.tpr" -po "${rep}/mdout_npt.mdp"
done
mpirun -np "${SLURM_NTASKS}" "${GMX}" mdrun -multidir "${REPS[@]}" -deffnm npt -ntomp "${SLURM_CPUS_PER_TASK}"

echo
echo "[9] Production input (unrestrained)"
for rep in "${REPS[@]}"; do
  gmx1 grompp -f "${MDP}/md.mdp" -c "${rep}/npt.gro" -t "${rep}/npt.cpt" \
    -p topol.top -n index.ndx -o "${rep}/md.tpr" -po "${rep}/mdout_md.mdp"
done

echo
echo "[10] Summary"
ls -lh start.pdb processed.gro topol.top em.gro index.ndx rep*/md.tpr
echo "Atoms in system: $(sed -n 2p em.gro)"

echo
echo "============================================================"
echo "Setup finished for ${SYSTEM}; next: sbatch 13_production.slurm ${SYSTEM}"
echo "============================================================"
