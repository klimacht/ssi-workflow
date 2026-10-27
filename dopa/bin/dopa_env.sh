# dopa_env.sh — load settings for the DOPA scripts.
#
#   source dopa/bin/dopa_env.sh SSI_ROOT
#
# Reads config/workflow.env without its SLURM_* lines: those are defaults for
# the SSI jobs (e.g. SLURM_CPUS_PER_TASK=8) and would overwrite the running
# job's own Slurm variables, which makes srun/mpirun abort. Then dopa/dopa.env.

_dopa_root="${1:?Usage: source dopa_env.sh SSI_ROOT}"
eval "$(grep -E '^[A-Za-z_][A-Za-z0-9_]*=' "${_dopa_root}/config/workflow.env" | grep -v '^SLURM_')"
source "${_dopa_root}/dopa/dopa.env"
unset _dopa_root
