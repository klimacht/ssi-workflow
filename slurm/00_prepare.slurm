#!/bin/bash -l
#SBATCH -o __SSI_ROOT__/logs/00_prepare.%j.out
#SBATCH -e __SSI_ROOT__/logs/00_prepare.%j.err
#SBATCH -J ssi_prepare
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=8
#SBATCH --mem=400G
#SBATCH --time=04:00:00
#SBATCH --mail-type=none

set -euo pipefail

if [ $# -lt 1 ]; then
  echo "Usage: sbatch 00_prepare.slurm PDBID [CHAIN]"
  echo "  CHAIN: optional single chain to keep, e.g. A (for multimeric crystals)"
  exit 1
fi

export PDBID=$(echo "$1" | tr '[:lower:]' '[:upper:]')
# Optional 2nd arg: chain to keep (e.g. A). Overrides SSI_CHAIN from config.
if [ $# -ge 2 ]; then
  export SSI_CHAIN="$2"
fi

source __SSI_ROOT__/config/workflow.env
# Allow SSI_CHAIN from config if not given on command line
export SSI_CHAIN="${SSI_CHAIN:-}"

source "${SSI_ROOT}/bin/ssi_lib.sh"
ssi_load_modules
ssi_setup_env

export RUN_DIR="${SSI_ROOT}/runs/${PDBID}"
export PDB_DIR="${RUN_DIR}/pdb"
export PREP_DIR="${RUN_DIR}/prep"

mkdir -p "${PDB_DIR}" "${PREP_DIR}" "${SSI_ROOT}/logs"
cd "${PREP_DIR}"

echo "============================================================"
echo "SSI prepare job"
echo "PDBID                 = ${PDBID}"
echo "SSI_CHAIN             = ${SSI_CHAIN:-all}"
echo "RUN_DIR               = ${RUN_DIR}"
echo "PREP_DIR              = ${PREP_DIR}"
echo "GMX                   = ${GMX}"
echo "SLURM_NTASKS          = ${SLURM_NTASKS}"
echo "SLURM_CPUS_PER_TASK   = ${SLURM_CPUS_PER_TASK}"
echo "OMP_NUM_THREADS       = ${OMP_NUM_THREADS}"
echo "============================================================"

echo
echo "[1] Obtain PDB"
if [ -s "${PDB_DIR}/${PDBID}.pdb" ]; then
  echo "Local PDB already present — skipping download:"
  ls -lh "${PDB_DIR}/${PDBID}.pdb"
else
  echo "Downloading ${PDBID} from RCSB..."
  wget -q -O "${PDB_DIR}/${PDBID}.pdb" "https://files.rcsb.org/download/${PDBID}.pdb" || true
fi

if [ ! -s "${PDB_DIR}/${PDBID}.pdb" ]; then
  echo "ERROR: no PDB available (download failed and no local file)."
  echo "       For a custom structure, place it at ${PDB_DIR}/${PDBID}.pdb before submitting."
  exit 1
fi

echo "Using PDB:"
ls -lh "${PDB_DIR}/${PDBID}.pdb"

echo
echo "[2] Clean PDB: keep ATOM records; resolve alternate locations; trim incomplete termini"
python3 - <<'INNERPY'
from pathlib import Path
from collections import defaultdict
import os
pdbid = os.environ["PDBID"]
inp = Path(os.environ["PDB_DIR"]) / f"{pdbid}.pdb"
out = Path("protein_only.pdb")
# Optional: keep only one chain (for multimeric crystal structures).
# Set SSI_CHAIN in the environment, e.g. SSI_CHAIN=A. Empty = keep all chains.
keep_chain = os.environ.get("SSI_CHAIN", "").strip()
BACKBONE_REQUIRED = {"N", "CA", "C"}
atoms = []
ter_lines = []
n_alt_skipped = 0
n_chain_skipped = 0
with inp.open() as fh:
    for line in fh:
        rec = line[0:6]
        if rec == "ATOM  ":
            chain = line[21] if len(line) > 21 else " "
            if keep_chain and chain != keep_chain:
                n_chain_skipped += 1
                continue
            altloc = line[16] if len(line) > 16 else " "
            if altloc not in (" ", "A"):
                n_alt_skipped += 1
                continue
            if len(line) >= 17:
                line = line[:16] + " " + line[17:]
            resseq = int(line[22:26])
            atoms.append((resseq, line.rstrip("\n") + "\n"))
        elif rec == "TER   ":
            chain = line[21] if len(line) > 21 else " "
            if keep_chain and chain != keep_chain:
                continue
            ter_lines.append(line)
res_atoms = defaultdict(set)
for resseq, line in atoms:
    res_atoms[resseq].add(line[12:16].strip())
all_resseqs = sorted(res_atoms.keys())
n_trimmed = 0
while all_resseqs:
    last = all_resseqs[-1]
    if {"N", "CA", "C"}.issubset(res_atoms[last]):
        break
    print(f"WARNING: Dropping C-terminal residue {last} "
          f"(missing backbone: {BACKBONE_REQUIRED - res_atoms[last]})")
    all_resseqs.pop()
    n_trimmed += 1
if not all_resseqs:
    raise SystemExit("ERROR: No residue with complete backbone found.")
keep = set(all_resseqs)
n_in = len(atoms)
n_skip_res = sum(1 for resseq, _ in atoms if resseq not in keep)
with out.open("w") as oh:
    for resseq, line in atoms:
        if resseq in keep:
            oh.write(line)
    if ter_lines:
        oh.write(ter_lines[0])
    oh.write("END\n")
print(f"Input ATOM records:          {n_in}")
print(f"Output ATOM records:         {n_in - n_skip_res}")
print(f"Skipped altloc records:      {n_alt_skipped}")
print(f"Skipped other-chain records: {n_chain_skipped}  (SSI_CHAIN={keep_chain or 'all'})")
print(f"Trimmed C-terminal residues: {n_trimmed}")
print(f"Last residue kept:           {max(keep)}")
INNERPY

echo "Number of ATOM lines:"
grep -c '^ATOM' protein_only.pdb || true

echo "Last 20 ATOM/TER/END lines of cleaned PDB:"
grep -E '^(ATOM|TER|END)' protein_only.pdb | tail -n 20

echo
echo "[3] pdb2gmx"
echo "Force field: ${FORCEFIELD}"
echo "Water model: ${WATERMODEL}"
mpirun -np 1 "${GMX}" pdb2gmx \
  -f protein_only.pdb \
  -o processed.gro \
  -p topol.top \
  -ff "${FORCEFIELD}" \
  -water "${WATERMODEL}" \
  -ignh -missing
  

echo
echo "[4] Create cubic box"
mpirun -np 1 "${GMX}" editconf \
  -f processed.gro \
  -o boxed.gro \
  -c \
  -d "${BOX_DISTANCE_NM}" \
  -bt cubic

echo
echo "[5] Solvate"
mpirun -np 1 "${GMX}" solvate \
  -cp boxed.gro \
  -cs spc216.gro \
  -o solvated.gro \
  -p topol.top

if [ ! -s topol.top ]; then
  echo "ERROR: topol.top missing or empty after solvate"
  exit 1
fi
cp -p topol.top topol.after_solvate.top
echo "Topology after solvate:"
ls -lh topol.top topol.after_solvate.top

echo
echo "[6] grompp for ions"
mpirun -np 1 "${GMX}" grompp \
  -f "${SSI_ROOT}/mdp/ions.mdp" \
  -c solvated.gro \
  -p topol.top \
  -o ions.tpr \
  -maxwarn 1

echo
echo "[7] Neutralize"
# Usually the solvent group is called SOL.
# If this fails, post the error output and we will adjust the group selection.
echo "SOL" | mpirun -np 1 "${GMX}" genion \
  -s ions.tpr \
  -o solv_ions.gro \
  -p topol.top \
  -pname NA \
  -nname CL \
  -neutral

if [ ! -s topol.top ]; then
  echo "ERROR: topol.top missing or empty after genion"
  exit 1
fi
cp -p topol.top topol.after_genion.top
echo "Topology after genion:"
ls -lh topol.top topol.after_genion.top
tail -n 30 topol.top

echo
echo "[8] Energy minimization grompp"
mpirun -np 1 "${GMX}" grompp \
  -f "${SSI_ROOT}/mdp/minim.mdp" \
  -c solv_ions.gro \
  -p topol.top \
  -o em.tpr \
  -maxwarn 1

echo
echo "[9] Energy minimization mdrun"
mpirun -np "${SLURM_NTASKS}" "${GMX}" mdrun \
  -v \
  -deffnm em \
  -ntomp "${SLURM_CPUS_PER_TASK}"

echo
echo "[10] Summary"
ls -lh protein_only.pdb processed.gro boxed.gro solvated.gro solv_ions.gro em.gro topol.top em.tpr em.edr em.log

echo
echo "============================================================"
echo "Preparation finished successfully for ${PDBID}"
echo "============================================================"
