#!/usr/bin/env bash
# Status overview of the DOPA study: RESP QM/fit, setup, production, analysis,
# with convergence numbers per replica and gnuplot terminal plots.
#
# Usage:
#   bash dopa/bin/dopa_monitor.sh                 all systems, RMSD plot per system
#   bash dopa/bin/dopa_monitor.sh wt              one system, all convergence plots
#   bash dopa/bin/dopa_monitor.sh --no-plots      numbers only
#   watch -n 120 bash dopa/bin/dopa_monitor.sh

ROOT="__SSI_ROOT__"
source "${ROOT}/dopa/bin/dopa_env.sh" "${ROOT}" 2>/dev/null

DETAIL=""; PLOTS=1
for a in "$@"; do
  case "$a" in
    --no-plots) PLOTS=0 ;;
    *)          DETAIL="$a" ;;
  esac
done

# python (MDAnalysis) + gnuplot for the convergence part
type module >/dev/null 2>&1 || source /etc/profile.d/modules.sh 2>/dev/null
module load "${MODULE_PYTHON}" "${MODULE_GNUPLOT:-gnuplot/6.0}" >/dev/null 2>&1
CONV=0; command -v python3 >/dev/null 2>&1 && CONV=1
PLOT_MODE="none"
if [ "$PLOTS" -eq 1 ] && command -v gnuplot >/dev/null 2>&1; then
  [ -n "$DETAIL" ] && PLOT_MODE="full" || PLOT_MODE="compact"
fi

QM_CONFS=(dah_alpha_a dah_alpha_b dah_beta_a dah_beta_b tyr_alpha tyr_beta)
QM_DIR="${DOPA_RUNS}/qm"
LOGS="${ROOT}/logs"
BAR_W=24
NOW_S=$(date +%s)

# ---- snapshot squeue once (keys: job name, or name#index for array tasks) ----
declare -A J_STATE J_TIME J_ID
while IFS='|' read -r jid name state time; do
  [[ "$name" == dopa_* ]] || continue
  if [[ "$jid" == *_\[* ]]; then                 # pending array range, e.g. 123_[1-5%4]
    range="${jid#*_[}"; range="${range%]}"; range="${range%%%*}"
    IFS=',' read -r -a parts <<< "$range"
    for p in "${parts[@]}"; do
      a="${p%-*}"; b="${p#*-}"
      for ((k=a; k<=b; k++)); do
        J_STATE["${name}#${k}"]="$state"; J_TIME["${name}#${k}"]="$time"; J_ID["${name}#${k}"]="${jid%%_*}_${k}"
      done
    done
  elif [[ "$jid" == *_* ]]; then                 # single array task
    k="${jid##*_}"
    J_STATE["${name}#${k}"]="$state"; J_TIME["${name}#${k}"]="$time"; J_ID["${name}#${k}"]="$jid"
  else
    J_STATE["$name"]="$state"; J_TIME["$name"]="$time"; J_ID["$name"]="$jid"
  fi
done < <(squeue -u "$USER" -h -o "%i|%j|%T|%M" 2>/dev/null)

bar() {   # bar FRACTION -> ████····
  awk -v f="$1" -v w="$BAR_W" 'BEGIN{ if (f < 0) f = 0; if (f > 1) f = 1; n = int(f * w + 0.5)
    s = ""; for (i = 0; i < w; i++) s = s (i < n ? "█" : "·"); print s }'
}

fmt_dur() {   # fmt_dur SECONDS -> "1d 05h" / "3h 12m" / "8m"
  awk -v s="$1" 'BEGIN{ s = int(s); d = int(s / 86400); h = int((s % 86400) / 3600); m = int((s % 3600) / 60)
    if (d > 0) printf "%dd %02dh", d, h; else if (h > 0) printf "%dh %02dm", h, m; else printf "%dm", m }'
}

last_log() {   # last_log PATTERN -> newest matching file in logs/
  ls -t ${LOGS}/$1 2>/dev/null | head -1
}

show_errors() {   # show_errors LOGFILE...
  local f
  for f in "$@"; do
    [ -s "$f" ] || continue
    grep -iE "error|fatal|abort|traceback|not found|missing|killed|cancelled" "$f" | tail -3 | sed 's/^/        /'
  done
}

NOW=$(date '+%H:%M:%S')
echo "=========================================================================="
echo "  DOPA MD Monitor  |  ${NOW}  |  Y${DOPA_MUT} -> DOPA, domain ${DOPA_FIRST}-${DOPA_LAST}"
echo "=========================================================================="

# ---------------------------------------------------------------------------
# RESP: QM per conformer + fit (overview only)
# ---------------------------------------------------------------------------
if [ -z "$DETAIL" ]; then
echo
echo "▶ RESP charges (HF/6-31G*, PySCF)"
for i in "${!QM_CONFS[@]}"; do
  c="${QM_CONFS[$i]}"
  d="${QM_DIR}/${c}"
  st="${J_STATE["dopa_qm#${i}"]:-}"
  info=""
  if [ -f "${d}/${c}_geomeTRIC_optim.xyz" ]; then   # comment lines: "Iteration N Energy E"
    info=$(grep "^Iteration" "${d}/${c}_geomeTRIC_optim.xyz" | tail -1 | awk '{printf "opt step %s, E = %s Eh", $2, $4}')
  fi
  if [ -s "${d}/${c}_vpot.out" ]; then
    status="DONE✓"
    info="$(head -n 2 "${d}/${c}.xyz" 2>/dev/null | tail -1 | sed -n 's/.*E=\([-0-9.]*\).*/E = \1 Eh/p'), $(head -n 1 "${d}/${c}_vpot.out") ESP points"
  elif [ -n "$st" ]; then
    status="$st"
    [ "$st" = "RUNNING" ] && status="RUN ${J_TIME["dopa_qm#${i}"]}"
  elif [ -e "${d}/${c}_pyscf.log" ]; then
    status="FAIL?"
  else
    status="-"
  fi
  printf "    %-12s %-12s %s\n" "$c" "$status" "$info"
  if [ "$status" = "FAIL?" ]; then
    show_errors "$(last_log "dopa_qm.*_${i}.err")" "$(last_log "dopa_qm.*_${i}.out")"
  fi
done

REPORT="${ROOT}/dopa/params/results/resp_report.txt"
fit_st="${J_STATE[dopa_fit]:-}"
if [ -s "$REPORT" ]; then
  maxd=$(grep "max |diff|" "$REPORT" | awk '{print $NF}')
  qoh=$(awk '$1 == "OH" && $2 == "OH" {print $3}' "$REPORT" | tail -1)
  qoe=$(awk '$1 == "OE2" {print $3}' "$REPORT" | tail -1)
  printf "    %-12s %-12s Tyr refit vs ff99SB: max |dq| = %s;  DAH q(OH) = %s, q(OE2) = %s\n" \
    "fit" "DONE✓" "$maxd" "$qoh" "$qoe"
elif [ -n "$fit_st" ]; then
  printf "    %-12s %-12s\n" "fit" "$fit_st"
else
  printf "    %-12s %-12s\n" "fit" "-"
fi
fi

# ---------------------------------------------------------------------------
# Systems: setup + production + convergence
# ---------------------------------------------------------------------------
for SYS in ${DETAIL:-${DOPA_SYSTEMS}}; do
  SD="${DOPA_RUNS}/${SYS}"
  echo
  if [ ! -d "$SD" ]; then
    printf "▶ %-8s setup:-\n" "$SYS"
    continue
  fi

  reps=()
  for r in $(seq 1 "${DOPA_NREP}"); do reps+=("${SD}/rep${r}"); done

  # setup stage from files on disk
  n_tpr=0
  for rd in "${reps[@]}"; do [ -s "${rd}/md.tpr" ] && n_tpr=$((n_tpr + 1)); done
  if   [ "$n_tpr" -eq "${#reps[@]}" ];  then stage="DONE✓"
  elif [ -s "${reps[0]}/npt.log" ];    then stage="NPT"
  elif [ -s "${reps[0]}/nvt.log" ];    then stage="NVT"
  elif [ -s "${SD}/em.gro" ];          then stage="EM✓"
  elif [ -s "${SD}/em.log" ];          then stage="EM"
  elif [ -s "${SD}/solv_ions.gro" ];   then stage="solvated"
  elif [ -s "${SD}/processed.gro" ];   then stage="topology"
  elif [ -s "${SD}/start.pdb" ];       then stage="built"
  else stage="-"; fi

  sst="${J_STATE["dopa_setup_${SYS}"]:-}"
  if   [ "$sst" = "RUNNING" ]; then setup="RUN:${stage}"
  elif [ -n "$sst" ];          then setup="$sst"
  elif [ "$stage" != "DONE✓" ] && [ "$stage" != "-" ]; then setup="FAIL?(${stage})"
  else setup="$stage"; fi

  atoms=""
  [ -s "${SD}/em.gro" ] && atoms="$(sed -n 2p "${SD}/em.gro" | tr -d ' ') atoms"

  mst="${J_STATE["dopa_md_${SYS}"]:-}"
  restarts=$(cat "${SD}/resubmit.count" 2>/dev/null || echo 0)
  prod="-"
  if   [ "$mst" = "RUNNING" ]; then prod="RUN ${J_TIME["dopa_md_${SYS}"]}"
  elif [ -n "$mst" ];          then prod="$mst"; fi

  printf "▶ %-8s setup:%-14s %-14s prod:%-14s restarts:%s\n" "$SYS" "$setup" "$atoms" "$prod" "$restarts"
  if [[ "$setup" == FAIL* ]]; then
    show_errors "$(last_log "dopa_setup_${SYS}.*.out")" "$(last_log "dopa_setup_${SYS}.*.err")"
  fi

  # benchmarks (rep*/bench_<N>n.log, N nodes; plain bench.log = 1 node)
  bench_line=""
  for b in $(ls "${SD}"/rep1/bench*.log 2>/dev/null | xargs -n1 basename 2>/dev/null); do
    perf=$(grep -h "^Performance:" "${SD}"/rep*/"$b" 2>/dev/null | awk '{s += $2; n++} END{if (n) printf "%.1f", s / n}')
    nn=$(echo "$b" | sed -n 's/^bench_\([0-9]*\)n\.log$/\1/p'); nn="${nn:-1}"
    bench_line+="  ${nn} node(s): ${perf:-running}"
  done
  bst="${J_STATE["dopa_bench_${SYS}"]:-}"
  [ -n "$bst" ] && bench_line+="  (job ${bst})"
  [ -n "$bench_line" ] && echo "    bench ns/day per replica:${bench_line}"

  if [ "$n_tpr" -gt 0 ]; then
  target_ps=$(awk -F= '/^nsteps/{n=$2} /^dt/{dt=$2} END{printf "%.0f", n * dt}' "${reps[0]}/mdout_md.mdp" 2>/dev/null)
  [ -n "$target_ps" ] && [ "$target_ps" -gt 0 ] 2>/dev/null || target_ps=300000

  sum_ps=0; sum_rate=0; n_rate=0; n_done=0
  for rd in "${reps[@]}"; do
    name=$(basename "$rd")
    if [ ! -s "${rd}/md.log" ]; then
      printf "    %-5s [%s]   not started\n" "$name" "$(bar 0)"
      continue
    fi
    read -r first last seg_first start perf < <(awk '
      /^Started mdrun on rank 0/ { sub(/^Started mdrun on rank 0 /, ""); gsub(/ /, "_"); start = $0; seg = "" }
      /^ +Step +Time/ { getline; if (first == "") first = $2; if (seg == "") seg = $2; last = $2 }
      /^Performance:/ { perf = $2 }
      END { print (first == "" ? 0 : first), (last == "" ? 0 : last), (seg == "" ? 0 : seg), (start == "" ? "-" : start), (perf == "" ? "-" : perf) }
    ' "${rd}/md.log")
    done_ps=$(awk -v a="$first" -v b="$last" 'BEGIN{print b - a}')
    sum_ps=$(awk -v s="$sum_ps" -v d="$done_ps" 'BEGIN{print s + d}')
    frac=$(awk -v d="$done_ps" -v t="$target_ps" 'BEGIN{print d / t}')

    rate=""
    if [ -s "${rd}/md.gro" ]; then
      n_done=$((n_done + 1))
    elif [ "$mst" = "RUNNING" ] && [ "$start" != "-" ]; then
      t0=$(date -d "${start//_/ }" +%s 2>/dev/null)
      if [ -n "$t0" ] && [ $((NOW_S - t0)) -gt 600 ]; then
        rate=$(awk -v a="$seg_first" -v b="$last" -v w=$((NOW_S - t0)) 'BEGIN{ r = (b - a) / 1000 / (w / 86400); if (r > 0) printf "%.1f", r }')
      fi
    fi
    [ -z "$rate" ] && [ "$perf" != "-" ] && rate="$perf"
    if [ -n "$rate" ]; then sum_rate=$(awk -v s="$sum_rate" -v r="$rate" 'BEGIN{print s + r}'); n_rate=$((n_rate + 1)); fi

    tag=""; [ -s "${rd}/md.gro" ] && tag="  done✓"
    printf "    %-5s [%s] %7.1f / %.0f ns%s%s\n" "$name" "$(bar "$frac")" \
      "$(awk -v d="$done_ps" 'BEGIN{print d / 1000}')" "$(awk -v t="$target_ps" 'BEGIN{print t / 1000}')" \
      "${rate:+   ${rate} ns/day}" "$tag"
  done

  total_ps=$(awk -v t="$target_ps" -v n="${#reps[@]}" 'BEGIN{print t * n}')
  line=$(awk -v s="$sum_ps" -v t="$total_ps" 'BEGIN{printf "%.1f / %.0f ns (%d%%)", s / 1000, t / 1000, 100 * s / t}')
  if [ "$n_rate" -gt 0 ] && [ "$n_done" -lt "${#reps[@]}" ]; then
    left=$(awk -v s="$sum_ps" -v t="$total_ps" -v r="$sum_rate" -v n="$n_rate" -v nr="${#reps[@]}" \
      'BEGIN{ per = r / n; rem = (t - s) / 1000 / nr; if (per > 0) print rem / per * 86400 }')
    [ -n "$left" ] && line+="   ≈ $(fmt_dur "$left") compute left (+ queue)"
  fi
  echo "    total $line"
  if [ -z "$mst" ] && [ "$n_done" -lt "${#reps[@]}" ] && [ -s "${reps[0]}/md.log" ]; then
    echo "    ⚠ production not in queue and not finished:"
    show_errors "$(last_log "dopa_md_${SYS}.*.out")" "$(last_log "dopa_md_${SYS}.*.err")"
  fi
  fi

  # convergence numbers (+ plots): energies from the logs, RMSD/Rg/ring distance from md.xtc
  # site of this system (y112_... -> 112; wt/dah_* -> DOPA_MUT); partner ring: Y112, or Y65
  SITE="${DOPA_MUT}"
  [[ "$SYS" =~ ^y([0-9]+)_ ]] && SITE="${BASH_REMATCH[1]}"
  PARTNER=112; [ "$SITE" != "65" ] && PARTNER=65
  [[ "$SYS" =~ ^x([0-9]+)_([0-9]+)$ ]] && SITE="${BASH_REMATCH[1]}" && PARTNER="${BASH_REMATCH[2]}"
  if [ "$CONV" -eq 1 ] && ls "${SD}"/rep*/nvt.log >/dev/null 2>&1; then
    python3 "${ROOT}/dopa/src/monitor_conv.py" --sysdir "$SD" --nrep "${DOPA_NREP}" \
      --first "${DOPA_FIRST}" --last "${DOPA_LAST}" --mut "${SITE}" --partner "${PARTNER}" \
      --plot "$PLOT_MODE" 2>&1 | grep -v "^\s*$"
  fi
done

# ---------------------------------------------------------------------------
# Analysis (overview only), one table per site
# ---------------------------------------------------------------------------
ast="${J_STATE[dopa_analysis]:-}"
n_summary=0
for SUMMARY in "${DOPA_RUNS}"/analysis/site*/summary.tsv; do
  [ -z "$DETAIL" ] && [ -s "$SUMMARY" ] || continue
  n_summary=$((n_summary + 1))
  systems=$(awk -F'\t' 'NR > 1 && !($2 in s) { s[$2] = 1; printf "%s ", $2 }' "$SUMMARY")
  echo
  echo "▶ Analysis $(basename "$(dirname "$SUMMARY")")  (${SUMMARY#${ROOT}/}, $(date -r "$SUMMARY" '+%d.%m. %H:%M'))${ast:+  job: $ast}"
  awk -F'\t' -v systems="${systems}" '
    NR == 1 { next }
    { key = $1 SUBSEP $2; val[key] = sprintf("%.3f±%.3f", $4, $5); if (!($1 in seen)) { seen[$1] = 1; order[++n] = $1 } }
    END {
      ns = split(systems, sys, " ")
      printf "    %-30s", "metric"; for (j = 1; j <= ns; j++) printf " %-16s", sys[j]; printf "\n"
      for (i = 1; i <= n; i++) {
        printf "    %-30s", order[i]
        for (j = 1; j <= ns; j++) { k = order[i] SUBSEP sys[j]; printf " %-16s", (k in val ? val[k] : "-") }
        printf "\n"
      }
    }' "$SUMMARY"
done
# exposure of the reactive Tyr atoms vs. the paper, Tyr-Tyr pairs that come close
for TVP in "${DOPA_RUNS}"/analysis/tyr_exposure/*/tyr_vs_paper.tsv; do
  [ -z "$DETAIL" ] && [ -s "$TVP" ] || continue
  n_summary=$((n_summary + 1))
  echo
  echo "▶ Tyr exposure $(basename "$(dirname "$TVP")")  (${TVP#${ROOT}/}, $(date -r "$TVP" '+%d.%m. %H:%M'))"
  echo "    ASA (Å²) paper → MD mean (fraction of time > 1 Å²); ND/HD (%) at ONOO⁻:Y 5:1"
  awk -F'\t' '
    function f1(x) { return (x == "NA" ? "NA" : sprintf("%.1f", x)) }
    function cell(p, m, a) { return sprintf("%s→%s (%.2f)", f1(p), f1(m), a) }
    BEGIN { printf "    %-5s %5s %5s  %-22s %-22s %-22s\n", "Tyr", "ND", "HD", "O", "ortho CE1+CE2", "meta CD1+CD2" }
    /^#/ || $1 == "tyr" { next }
    { printf "    Y%-4s %5s %5s  %-22s %-22s %-22s\n", $1, f1($3), f1($5),
        cell($9, $10, $11), cell($12, $13, $14), cell($15, $16, $17) }' "$TVP"
  PAIRS="$(dirname "$TVP")/pairs.tsv"
  [ -s "$PAIRS" ] && awk -F'\t' '
    /^#/ || $1 == "pair" { next }
    $9 + 0 < 8 { p = $1; sub("-", "–Y", p)
      printf "    Y%-11s ortho C···C %4.1f Å, ≤ 4.5 Å in %3.0f %% of frames  (cross-link: %s, NMR contact: %s)\n",
             p, $9, 100 * $12, $2, $4 }' "$PAIRS"
done
if [ -z "$DETAIL" ] && [ "$n_summary" -eq 0 ]; then
  echo
  echo "▶ Analysis  ${ast:--}"
fi

# ---------------------------------------------------------------------------
# Active jobs
# ---------------------------------------------------------------------------
echo
echo "--- Active jobs ---"
any=0
while IFS='|' read -r jid name state time endt reason; do
  [[ "$name" == dopa_* ]] || continue
  any=1
  printf "  %-16s %-18s %-9s %-10s wall<=%-19s %s\n" "$jid" "$name" "$state" "$time" "$endt" "$reason"
done < <(squeue -u "$USER" -h -o "%i|%j|%T|%M|%e|%R" 2>/dev/null)
[ "$any" -eq 1 ] || echo "  (none)"
echo "=========================================================================="
