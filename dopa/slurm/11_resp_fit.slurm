#!/bin/bash -l
#SBATCH -o __SSI_ROOT__/logs/%x.%j.out
#SBATCH -e __SSI_ROOT__/logs/%x.%j.err
#SBATCH -J dopa_fit
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --mail-type=none

# Two-stage RESP fit over all conformers; writes dopa/params/results/ and
# dopa/ff/dah.rtp + dah.hdb (used by dopa/ff/make_ff.sh).

set -euo pipefail

source __SSI_ROOT__/dopa/bin/dopa_env.sh __SSI_ROOT__

module purge
module load "${MODULE_PYTHON}"

python3 "${SSI_ROOT}/dopa/params/resp_dopa.py" fit \
  --workdir "${DOPA_RUNS}/qm" \
  --outdir "${SSI_ROOT}/dopa/params/results" \
  --ffdir "${SSI_ROOT}/dopa/ff"
