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
# side-chain SASA of every site in DOPA_SITES) and compare, per site, the
# wild type with all systems modified at that site (runs/dopa/analysis/site<N>).
# Then the time-resolved exposure of the reactive Tyr atoms and the Tyr-Tyr
# contacts of DOPA_EXPOSURE_SYSTEMS vs. the paper (runs/dopa/analysis/tyr_exposure).
#
#   sbatch 14_analysis.slurm                  one node, replicas one after another
#   sbatch --nodes=N ... 14_analysis.slurm    every replica as its own job step, side by side
#                                             (dopa_submit.sh analysis with DOPA_ANALYSIS_NODES=N)
#   14_analysis.slurm prep SYSTEM REP         one replica; what each job step runs

set -euo pipefail

source __SSI_ROOT__/dopa/bin/dopa_env.sh __SSI_ROOT__

source "${SSI_ROOT}/bin/ssi_lib.sh"
ssi_load_modules python
ssi_setup_env

gmx1() { mpirun -np 1 "${GMX}" "$@"; }

prep_rep() {   # prep_rep SYSTEM REP
  local SYSTEM="$1" r="$2" SITE
  echo "=== ${SYSTEM}/rep${r}"
  cd "${DOPA_RUNS}/${SYSTEM}"      # residuetypes.dat with DAH/TOP lives here
  echo Protein | gmx1 convert-tpr -s "rep${r}/md.tpr" -n index.ndx -o "rep${r}/md_protein.tpr"
  echo Protein | gmx1 trjconv -s "rep${r}/md.tpr" -f "rep${r}/npt.gro" -n index.ndx \
    -o "rep${r}/protein_ref.gro" -pbc mol
  printf 'Protein\nProtein\n' | gmx1 trjconv -s "rep${r}/md_protein.tpr" -f "rep${r}/md.xtc" \
    -o "rep${r}/md_center.xtc" -pbc mol -center
  for SITE in ${DOPA_SITES}; do
    gmx1 sasa -s "rep${r}/md_protein.tpr" -f "rep${r}/md_center.xtc" -surface all \
      -output "resnr ${SITE} and not name N H CA HA C O" \
      -o "rep${r}/sasa_res${SITE}.xvg"
  done
}

if [ "${1:-}" = "prep" ]; then
  # inside a job step Slurm's PMI allows one MPI start only, so each gmx call gets its own local
  # mpirun (fork bootstrap instead of srun)
  gmx1() { I_MPI_HYDRA_BOOTSTRAP=fork mpirun -n 1 -hosts "$(hostname)" "${GMX}" "$@"; }
  prep_rep "${2:?}" "${3:?}"
  exit 0
fi

NODES="${SLURM_NNODES:-1}"
REPS=()
for SYSTEM in ${DOPA_SYSTEMS}; do
  for r in $(seq 1 "${DOPA_NREP}"); do
    if [ -s "${DOPA_RUNS}/${SYSTEM}/rep${r}/md.xtc" ]; then
      REPS+=("${SYSTEM}:${r}")
    else
      echo "skip ${SYSTEM}/rep${r}: no md.xtc"
    fi
  done
done

if [ "${NODES}" -gt 1 ]; then
  # one job step per replica, spread evenly over the nodes (gmx sasa/trjconv are serial)
  STEP_CPUS=$(( NODES * SLURM_CPUS_PER_TASK / ${#REPS[@]} / 2 * 2 ))
  [ "${STEP_CPUS}" -ge 2 ] || STEP_CPUS=2
  echo "${#REPS[@]} replicas as job steps on ${NODES} nodes, ${STEP_CPUS} CPUs each"
  PIDS=()
  for sr in "${REPS[@]}"; do
    SYSTEM="${sr%%:*}"; r="${sr##*:}"
    srun --exact --nodes=1 --ntasks=1 --cpus-per-task="${STEP_CPUS}" --mem=8G \
      bash "${SSI_ROOT}/dopa/slurm/14_analysis.slurm" prep "${SYSTEM}" "${r}" \
      > "${DOPA_RUNS}/${SYSTEM}/rep${r}/analysis_prep.log" 2>&1 &
    PIDS+=($!)
  done
  FAIL=0
  for p in "${PIDS[@]}"; do wait "$p" || FAIL=1; done
  if [ "${FAIL}" != 0 ]; then
    echo "ERROR: a replica failed, see <system>/rep<i>/analysis_prep.log" >&2
    exit 1
  fi
else
  for sr in "${REPS[@]}"; do
    prep_rep "${sr%%:*}" "${sr##*:}"
  done
fi

# per site WT vs. all systems modified there, and the Tyr exposure, side by side
export OMP_NUM_THREADS=$(( SLURM_CPUS_PER_TASK / 4 > 16 ? 16 : (SLURM_CPUS_PER_TASK / 4 > 0 ? SLURM_CPUS_PER_TASK / 4 : 1) ))
mkdir -p "${DOPA_RUNS}/analysis"
LOGS=()
PIDS=()
for SITE in ${DOPA_SITES}; do
  GROUP="wt"
  for s in ${DOPA_SYSTEMS}; do
    [ "$s" = "wt" ] && continue
    read -r site _ _ < <(python3 "${SSI_ROOT}/dopa/src/make_system_pdb.py" --parse "$s" --default-mut "${DOPA_MUT}")
    [ "$site" = "$SITE" ] && GROUP+=" $s"
  done
  [ "$GROUP" = "wt" ] && continue
  PARTNERS=""
  for c in ${DOPA_CLUSTER}; do [ "$c" != "$SITE" ] && PARTNERS+=" $c"; done
  HELIX=$(echo "${DOPA_SITE_HELIX:-}" | tr ' ' '\n' | awk -F= -v s="${SITE}" '$1 == s {print $2}')
  mkdir -p "${DOPA_RUNS}/analysis/site${SITE}"
  LOGS+=("${DOPA_RUNS}/analysis/site${SITE}/log.txt")
  echo "=== Comparison for site ${SITE}: ${GROUP}"
  python3 "${SSI_ROOT}/dopa/src/analyze_dopa.py" \
    --runs "${DOPA_RUNS}" --systems ${GROUP} --nrep "${DOPA_NREP}" \
    --first "${DOPA_FIRST}" --last "${DOPA_LAST}" --mut "${SITE}" --partners ${PARTNERS} \
    ${HELIX:+--helix "${HELIX}"} --skip-ns "${DOPA_SKIP_NS}" --out "${DOPA_RUNS}/analysis/site${SITE}" \
    > "${LOGS[-1]}" 2>&1 &
  PIDS+=($!)
done

# exposure of the reactive Tyr atoms and Tyr-Tyr contacts vs. the paper (Table 1, modification degrees)
echo "=== Tyrosine exposure and Tyr-Tyr contacts: ${DOPA_EXPOSURE_SYSTEMS:-wt}"
mkdir -p "${DOPA_RUNS}/analysis/tyr_exposure"
LOGS+=("${DOPA_RUNS}/analysis/tyr_exposure/log.txt")
{
  python3 "${SSI_ROOT}/dopa/src/tyr_exposure.py" static --pdb "${DOPA_SRC_PDB}"
  python3 "${SSI_ROOT}/dopa/src/tyr_exposure.py" md \
    --runs "${DOPA_RUNS}" --systems ${DOPA_EXPOSURE_SYSTEMS:-wt} --nrep "${DOPA_NREP}" \
    --skip-ns "${DOPA_SKIP_NS}" --out "${DOPA_RUNS}/analysis/tyr_exposure"
} > "${LOGS[-1]}" 2>&1 &
PIDS+=($!)

FAIL=0
for p in "${PIDS[@]}"; do wait "$p" || FAIL=1; done
for l in "${LOGS[@]}"; do
  echo
  echo "--- ${l}"
  cat "${l}"
done
exit "${FAIL}"
