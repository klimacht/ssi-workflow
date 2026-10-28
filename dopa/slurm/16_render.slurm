#!/bin/bash -l
#SBATCH -o __SSI_ROOT__/logs/%x.%j.out
#SBATCH -e __SSI_ROOT__/logs/%x.%j.err
#SBATCH -J dopa_render
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=30G
#SBATCH --time=01:30:00
#SBATCH --mail-type=none

# PyMOL figures and movies from the analysis (runs/dopa/analysis -> analysis/figures), in three
# steps as for the SSI maps: structures and movie frames (viz_prep.py, MDAnalysis), raw images
# (render_dopa_pymol.py, PyMOL from DOPA_PYMOL_ENV), panels, colorbars and MP4 (render_dopa_compose.py,
# matplotlib + ffmpeg). Needs the output of 14_analysis.slurm.

set -euo pipefail

source __SSI_ROOT__/dopa/bin/dopa_env.sh __SSI_ROOT__
SRC="${SSI_ROOT}/dopa/src"
AN="${DOPA_RUNS}/analysis"
VIZ="${AN}/viz"

echo "=== Structures and movie frames"
module purge
module load "${MODULE_PYTHON}"
python3 "${SRC}/viz_prep.py" --runs "${DOPA_RUNS}" --analysis "${AN}" --systems ${DOPA_SYSTEMS} \
  --nrep "${DOPA_NREP}" --skip-ns "${DOPA_SKIP_NS}" --out "${VIZ}"

echo "=== PyMOL"
module purge
set +u                              # conda's activation scripts use unset variables
source "${DOPA_CONDA_SH}"
conda activate "${DOPA_PYMOL_ENV}"
set -u
python "${SRC}/render_dopa_pymol.py" --viz "${VIZ}" --out "${VIZ}/raw"
set +u
conda deactivate
set -u

echo "=== Figures and movies"
module purge
module load "${MODULE_GCC}" "${DOPA_MODULE_FFMPEG}" "${MODULE_PYTHON}"
python3 "${SRC}/render_dopa_compose.py" --viz "${VIZ}" --raw "${VIZ}/raw" --out "${AN}/figures"
