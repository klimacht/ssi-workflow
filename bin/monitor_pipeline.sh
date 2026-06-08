#!/usr/bin/env bash
# Smart pipeline monitor for one or more PDBIDs.
#
# Usage:
#   bash bin/monitor_pipeline.sh                  # all PDBIDs under runs/
#   bash bin/monitor_pipeline.sh 1DPX 3LDJ        # specific ones
#   watch -n 30 bash bin/monitor_pipeline.sh 3LDJ 1RBB

ROOT="__SSI_ROOT__"
source "${ROOT}/config/workflow.env" 2>/dev/null

PHI_ARRAY=(${PHI_LIST})
NPHI=${#PHI_ARRAY[@]}
PROD_NSTEPS=1500000   # 3 ns at 2 fs — for % progress
PROD_TARGET_PS=3000

# ---- args: PDBIDs ----
if [ $# -eq 0 ]; then
  PDBIDS=()
  for d in "${ROOT}/runs"/*/; do
    [ -d "$d" ] && PDBIDS+=( "$(basename "$d")" )
  done
else
  PDBIDS=()
  for a in "$@"; do PDBIDS+=( "$(echo "$a" | tr '[:lower:]' '[:upper:]')" ); done
fi

# ---- snapshot squeue once ----
declare -A SQ_STATE
declare -A SQ_TIME
declare -A SQ_END
while IFS='|' read -r jid state time endt; do
  SQ_STATE[$jid]="$state"
  SQ_TIME[$jid]="$time"
  SQ_END[$jid]="$endt"
done < <(squeue -u "$USER" -h -o "%i|%T|%M|%e" 2>/dev/null)

# Precompute phi dir names once (avoids per-call python in watch loops)
PHIDIR_NAMES=()
for ((i=0;i<NPHI;i++)); do
  phi="${PHI_ARRAY[$i]}"
  w=${phi%.*}
  frac="${phi#*.}"; [ "$frac" = "$phi" ] && frac=0
  # handle integer phi like "5" -> frac 0; "5.5" -> frac 5
  PHIDIR_NAMES+=( "$(printf 'phi_%02dp%d' "$w" "${frac:0:1}")" )
done

phidir() { echo "${PHIDIR_NAMES[$1]}"; }

# helper: last COLVAR time for a phi dir
colvar_time() {
  local f="$1"
  [ -f "$f" ] || { echo ""; return; }
  awk '!/^#/{t=$1} END{print t}' "$f" 2>/dev/null
}

# helper: slurm state for a jobid (handles array base)
job_state() {
  local jid="$1"
  if [ -n "${SQ_STATE[$jid]+_}" ]; then echo "${SQ_STATE[$jid]}"; return; fi
  local base="${jid%%_*}"
  # any array task active?
  for k in "${!SQ_STATE[@]}"; do
    [[ "$k" == ${base}_* ]] && { echo "${SQ_STATE[$k]}"; return; }
  done
  echo "GONE"
}

NOW=$(date '+%H:%M:%S')
echo "=========================================================================="
echo "  SSI Pipeline Monitor  |  ${NOW}  |  ${#PDBIDS[@]} protein(s)"
echo "=========================================================================="

stages=(prepare equilibrate unbiased ssi_preflight ssi_production analysis)
labels=(prep equil unbias preflt prod analy)

for PDBID in "${PDBIDS[@]}"; do
  RUN_DIR="${ROOT}/runs/${PDBID}"
  [ -d "$RUN_DIR" ] || continue

  # ---- stage state line ----
  line=""
  prod_done=0
  for si in "${!stages[@]}"; do
    stage="${stages[$si]}"
    lbl="${labels[$si]}"
    jf="${RUN_DIR}/${stage}.jobid"

    if [ ! -f "$jf" ]; then
      line+=$(printf " %-7s:%-9s" "$lbl" "-")
      continue
    fi
    jid=$(cat "$jf")
    st=$(job_state "$jid")

    # Refine "GONE" using on-disk evidence
    if [ "$st" = "GONE" ]; then
      case "$stage" in
        ssi_production)
          # count complete phi windows
          done_n=0
          for ((i=0;i<NPHI;i++)); do
            pd=$(phidir $i)
            t=$(colvar_time "${RUN_DIR}/ssi_production/${pd}/COLVAR")
            [ -n "$t" ] && awk "BEGIN{exit !($t>=$PROD_TARGET_PS)}" && done_n=$((done_n+1))
          done
          if [ "$done_n" -eq "$NPHI" ]; then st="DONE✓"; prod_done=1
          else st="${done_n}/${NPHI}"; fi
          ;;
        analysis)
          [ -f "${RUN_DIR}/dewetting"/*.pdb ] 2>/dev/null && st="DONE✓" || \
          { ls "${RUN_DIR}/dewetting"/*.pdb >/dev/null 2>&1 && st="DONE✓" || st="?"; }
          ;;
        *)
          # stage-specific output file existence
          case "$stage" in
            prepare)      [ -f "${RUN_DIR}/prep/em.gro" ] && st="DONE✓" || st="FAIL?";;
            equilibrate)  [ -f "${RUN_DIR}/equil/npt.gro" ] && st="DONE✓" || st="FAIL?";;
            unbiased)     [ -f "${RUN_DIR}/unbiased/prod_unbiased.gro" ] && st="DONE✓" || st="FAIL?";;
            ssi_preflight)[ -f "${RUN_DIR}/ssi_preflight/phi_00p0/ssi_phi.tpr" ] && st="DONE✓" || st="FAIL?";;
          esac
          ;;
      esac
    fi
    line+=$(printf " %-7s:%-9s" "$lbl" "$st")
  done

  echo
  printf "▶ %-6s%s\n" "$PDBID" "$line"

  # ---- production phi detail (only if not fully done and dir exists) ----
  if [ -d "${RUN_DIR}/ssi_production" ] && [ "$prod_done" -eq 0 ]; then
    detail=""
    for ((i=0;i<NPHI;i++)); do
      pd=$(phidir $i)
      t=$(colvar_time "${RUN_DIR}/ssi_production/${pd}/COLVAR")
      if [ -z "$t" ]; then
        sym="·"
      elif awk "BEGIN{exit !($t>=$PROD_TARGET_PS)}"; then
        sym="█"
      else
        pct=$(awk "BEGIN{printf \"%d\", $t/$PROD_TARGET_PS*100}")
        if   [ "$pct" -ge 75 ]; then sym="▓"
        elif [ "$pct" -ge 50 ]; then sym="▒"
        elif [ "$pct" -ge 25 ]; then sym="░"
        else sym="."; fi
      fi
      detail+="$sym"
    done
    echo "         phi[0..12]: ${detail}   (· none  . <25%  ░▒▓ part  █ done)"
  fi

  # ---- show error tail only for FAIL? stages ----
  for si in "${!stages[@]}"; do
    stage="${stages[$si]}"
    jf="${RUN_DIR}/${stage}.jobid"
    [ -f "$jf" ] || continue
    jid=$(cat "$jf")
    st=$(job_state "$jid")
    [ "$st" != "GONE" ] && continue
    # re-derive fail
    fail=0
    case "$stage" in
      prepare)      [ -f "${RUN_DIR}/prep/em.gro" ] || fail=1;;
      equilibrate)  [ -f "${RUN_DIR}/equil/npt.gro" ] || fail=1;;
      unbiased)     [ -f "${RUN_DIR}/unbiased/prod_unbiased.gro" ] || fail=1;;
      ssi_preflight)[ -f "${RUN_DIR}/ssi_preflight/phi_00p0/ssi_phi.tpr" ] || fail=1;;
    esac
    [ "$fail" -eq 1 ] || continue
    err=$(ls "${ROOT}/logs/"*"${jid}"*.err 2>/dev/null | head -1)
    out=$(ls "${ROOT}/logs/"*"${jid}"*.out 2>/dev/null | head -1)
    echo "    ⚠ ${stage} (${jid}) appears FAILED:"
    if [ -f "$err" ] && [ -s "$err" ]; then
      grep -iE "fatal|error|abort|not found|missing" "$err" | tail -3 | sed 's/^/        /'
    elif [ -f "$out" ]; then
      grep -iE "ERROR|missing|not found" "$out" | tail -3 | sed 's/^/        /'
    fi
  done
done

# ---- active jobs summary with ETA ----
echo
echo "--- Active jobs (with estimated end) ---"
any_active=0
while IFS='|' read -r jid name state time endt; do
  [ -z "$jid" ] && continue
  any_active=1
  printf "  %-12s %-12s %-9s %-8s end:%s\n" "$jid" "$name" "$state" "$time" "$endt"
done < <(squeue -u "$USER" -h -o "%i|%j|%T|%M|%e" 2>/dev/null | grep -iE "ssi|prep|equil|unbias|prod|anal")
[ "$any_active" -eq 0 ] && echo "  (none)"

echo "=========================================================================="
