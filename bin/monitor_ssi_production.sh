#!/usr/bin/env bash
# Monitor SSI production array jobs for ssi_workflow.
# Jobs are sorted by phi value.

ROOT="__SSI_ROOT__"
source "${ROOT}/config/workflow.env"

PDBID="${1:-}"
if [ -z "$PDBID" ]; then
  # Try to find PDBID from most recent jobid file
  PDBID=$(ls -t "${ROOT}/runs"/*/ssi_production.jobid 2>/dev/null | head -1 | awk -F'/' '{print $(NF-1)}')
fi
if [ -z "$PDBID" ]; then
  echo "Usage: $0 PDBID"
  exit 1
fi
PDBID=$(echo "$PDBID" | tr '[:lower:]' '[:upper:]')

JOBID_FILE="${ROOT}/runs/${PDBID}/ssi_production.jobid"
ARRAY_JOBID=""
if [ -f "$JOBID_FILE" ]; then
  ARRAY_JOBID=$(cat "$JOBID_FILE")
fi

PHI_ARRAY=(${PHI_LIST})

echo "=========================================================================="
echo "  SSI Production Monitor — ${PDBID}"
[ -n "$ARRAY_JOBID" ] && echo "  Array Job ID: ${ARRAY_JOBID}"
echo "  $(date)"
echo "=========================================================================="
echo

printf "%-6s %-12s %-18s %-10s %-10s %-10s %-12s %-12s %-10s %-10s %-12s %-10s\n" \
  "IDX" "JOBID" "PHIDIR" "STATE" "WALL" "STEP" "TIME(ps)" "TEMP(K)" "N" "Ntilde" "Bias" "DONE?"
printf "%-6s %-12s %-18s %-10s %-10s %-10s %-12s %-12s %-10s %-10s %-12s %-10s\n" \
  "------" "------------" "------------------" "----------" "----------" "----------" "------------" "------------" "----------" "----------" "------------" "----------"

# Build a map of array_index -> slurm state/jobid from squeue
declare -A SLURM_STATE
declare -A SLURM_WALL
declare -A SLURM_JOBID_FOR_IDX

if [ -n "$ARRAY_JOBID" ]; then
  while IFS='|' read -r FULL_JID NAME STATE WALL; do
    IDX="${FULL_JID##*_}"
    SLURM_STATE[$IDX]="$STATE"
    SLURM_WALL[$IDX]="$WALL"
    SLURM_JOBID_FOR_IDX[$IDX]="$FULL_JID"
  done < <(squeue -u "$USER" -h -o "%i|%j|%T|%M" | grep "^${ARRAY_JOBID}_")
fi

# Iterate phi values in sorted order (they are already 0..12 in PHI_LIST)
for IDX in "${!PHI_ARRAY[@]}"; do
  PHI="${PHI_ARRAY[$IDX]}"

  WHOLE=$(python3 -c "phi=float('${PHI}'); print(int(phi))")
  FRAC=$(python3 -c "phi=float('${PHI}'); w=int(phi); print(int(round((phi-w)*10)))")
  PHIDIR_NAME=$(printf "phi_%02dp%d" "$WHOLE" "$FRAC")

  WORK_DIR="${ROOT}/runs/${PDBID}/ssi_production/${PHIDIR_NAME}"

  STATE="${SLURM_STATE[$IDX]:-}"
  WALL="${SLURM_WALL[$IDX]:--}"
  FULL_JID="${SLURM_JOBID_FOR_IDX[$IDX]:--}"

  STEP="--"; TIMEPS="--"; TEMP="--"; PRESS="--"
  NVAL="--"; NTILDE="--"; BIAS="--"
  DONE="no"

  LOGFILE="${WORK_DIR}/ssi_phi.log"
  COLVAR="${WORK_DIR}/COLVAR"

  if [ -f "$LOGFILE" ]; then
    read -r STEP TIMEPS TEMP PRESS < <(
      awk '
        /^[[:space:]]*Step[[:space:]]+Time/ { getline; step=$1; time=$2 }
        /Temperature/ && /Pressure/ { getline; temp=$4; press=$5 }
        END {
          if(step=="") step="--"
          if(time=="") time="--"
          if(temp=="") temp="--"
          if(press=="") press="--"
          print step, time, temp, press
        }
      ' "$LOGFILE"
    )
  fi

  if [ -f "$COLVAR" ]; then
    read -r NVAL NTILDE BIAS < <(
      awk '$1 !~ /^#/ { n=$2; nt=$3; b=$4 }
        END {
          if(n=="") n="--"
          if(nt=="") nt="--"
          if(b=="") b="--"
          print n, nt, b
        }' "$COLVAR"
    )
  fi

  # Done = gro + cpt exist and job not running
  if [ -f "${WORK_DIR}/ssi_phi.gro" ] && [ -f "${WORK_DIR}/ssi_phi.cpt" ] && [ -z "$STATE" ]; then
    DONE="DONE"
  elif [ -n "$STATE" ]; then
    DONE="$STATE"
  elif [ -f "${WORK_DIR}/ssi_phi.tpr" ]; then
    DONE="PENDING?"
  fi

  printf "%-6s %-12s %-18s %-10s %-10s %-10s %-12s %-12s %-10s %-10s %-12s %-10s\n" \
    "$IDX" "$FULL_JID" "$PHIDIR_NAME" "$STATE" "$WALL" \
    "$STEP" "$TIMEPS" "$TEMP" "$NVAL" "$NTILDE" "$BIAS" "$DONE"
done

echo
echo "=== COLVAR N-values per phi (last line) ==========================="
for IDX in "${!PHI_ARRAY[@]}"; do
  PHI="${PHI_ARRAY[$IDX]}"
  WHOLE=$(python3 -c "phi=float('${PHI}'); print(int(phi))")
  FRAC=$(python3 -c "phi=float('${PHI}'); w=int(phi); print(int(round((phi-w)*10)))")
  PHIDIR_NAME=$(printf "phi_%02dp%d" "$WHOLE" "$FRAC")
  COLVAR="${ROOT}/runs/${PDBID}/ssi_production/${PHIDIR_NAME}/COLVAR"
  if [ -f "$COLVAR" ]; then
    LAST=$(awk '$1 !~ /^#/ {line=$0} END {print line}' "$COLVAR")
    printf "  phi=%-5s  %s\n" "$PHI" "$LAST"
  else
    printf "  phi=%-5s  (no COLVAR yet)\n" "$PHI"
  fi
done

echo
echo "=== Progress (stderr) ============================================="
for IDX in "${!PHI_ARRAY[@]}"; do
  FULL_JID="${SLURM_JOBID_FOR_IDX[$IDX]:-}"
  [ -z "$FULL_JID" ] || [ "$FULL_JID" = "-" ] && continue
  ERRFILE="${ROOT}/logs/04_ssi_production.${FULL_JID}.err"
  if [ -f "$ERRFILE" ]; then
    PROG=$(grep -E "step [0-9]|will finish|Performance" "$ERRFILE" | tail -1)
    printf "  %-12s phi=%-5s  %s\n" "$FULL_JID" "${PHI_ARRAY[$IDX]}" "${PROG:-(no progress line yet)}"
  fi
done
