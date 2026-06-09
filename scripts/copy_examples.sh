#!/usr/bin/env bash
# Copy completed run results into examples/ for distribution.
# Drops trajectories, keeps maps + per-atom data + N_vs_phi.tsv.
set -euo pipefail

SRC=/ptmp/tklimach/ssi_workflow/runs
DST=/ptmp/tklimach/ssi-workflow/examples
PROTEINS=(1DPX 1HRC 1FS3 1RBB 3LDJ 4IBA PHLP5)

mkdir -p "$DST"

for PDBID in "${PROTEINS[@]}"; do
  RUN="$SRC/$PDBID"
  [ -d "$RUN" ] || { echo "  skip $PDBID (no run dir)"; continue; }

  out="$DST/$PDBID"
  mkdir -p "$out/analysis" "$out/dewetting"

  # cleaned input structure
  if [ -f "$RUN/pdb/${PDBID}.pdb" ]; then
    cp -p "$RUN/pdb/${PDBID}.pdb" "$out/${PDBID}_input.pdb"
  fi


  # analysis tables
  for f in N_vs_phi.tsv n0_per_atom.npy; do
    [ -f "$RUN/analysis/$f" ] && cp -p "$RUN/analysis/$f" "$out/analysis/$f"
  done

  # dewetting maps
  for f in "$RUN"/dewetting/dewetting_phi_*.pdb \
           "$RUN"/dewetting/dewetting_phi_*.dat \
           "$RUN"/dewetting/dewetting_phistar.pdb \
           "$RUN"/dewetting/dewetting_phistar.dat; do
    [ -f "$f" ] && cp -p "$f" "$out/dewetting/$(basename "$f")"
  done

  echo "  copied $PDBID"
done

# Top-level comparison table
[ -f /ptmp/tklimach/ssi_workflow/analysis_comparison.tsv ] && \
  cp -p /ptmp/tklimach/ssi_workflow/analysis_comparison.tsv "$DST/"

echo
echo "Done. Sizes:"
du -sh "$DST"/* 2>/dev/null | sort -k2
