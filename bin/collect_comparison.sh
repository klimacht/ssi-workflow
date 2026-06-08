#!/usr/bin/env bash
# Collect a cross-protein comparison table from completed analyses.
set -euo pipefail
ROOT="__SSI_ROOT__"
source "${ROOT}/config/workflow.env"

PROTEINS=("$@")
if [ ${#PROTEINS[@]} -eq 0 ]; then
  PROTEINS=(1DPX 1HRC 3LDJ 1FS3 1RBB 4IBA PHLP5)
fi

OUT="${ROOT}/analysis_comparison.tsv"
python3 - "$OUT" "${PROTEINS[@]}" << 'PYEOF'
import sys, csv
from pathlib import Path
import numpy as np

out = sys.argv[1]
proteins = sys.argv[2:]
root = Path("__SSI_ROOT__")

rows = []
for pdbid in proteins:
    adir = root / "runs" / pdbid / "analysis"
    tsv = adir / "N_vs_phi.tsv"
    if not tsv.is_file():
        print(f"  {pdbid}: no N_vs_phi.tsv, skipping")
        continue
    n0 = phistar_g = var = nstar = None
    with open(tsv) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if float(r["phi"]) == 0.0:
                n0 = float(r["N_mean"])
            if r.get("phi_star","0") == "1":
                phistar_g = float(r["phi"]); var = float(r["N_var"]); nstar = float(r["N_mean"])
    gamma = var / n0 if (var and n0) else float("nan")

    # patch contrast from per-atom phistar .dat (column 9 = phistar, may be centered)
    contrast = float("nan"); n_surf = 0
    dat = root / "runs" / pdbid / "dewetting" / "dewetting_phistar.dat"
    if dat.is_file():
        vals = []
        with open(dat) as f:
            for line in f:
                if line.startswith("#"): continue
                p = line.split()
                if len(p) > 10 and p[10] == "1" and p[8] != "nan":
                    vals.append(float(p[8]))
        if vals:
            v = np.array(vals); n_surf = len(v)
            contrast = float(np.percentile(v,95) - np.percentile(v,5))

    rows.append((pdbid, n0, phistar_g, gamma, contrast, n_surf))

rows.sort(key=lambda r: (-(r[3] if r[3]==r[3] else -1)))  # by gamma desc

hdr = ["protein","N0","phi_star","gamma","patch_contrast_kJmol","n_surface"]
with open(out,"w") as f:
    f.write("\t".join(hdr)+"\n")
    for r in rows:
        f.write(f"{r[0]}\t{r[1]:.1f}\t{r[2]:.1f}\t{r[3]:.3f}\t{r[4]:.2f}\t{r[5]}\n")

print(f"{'protein':8} {'N0':>7} {'phi*':>5} {'gamma':>6} {'contrast':>9} {'n_surf':>7}")
for r in rows:
    print(f"{r[0]:8} {r[1]:7.0f} {r[2]:5.1f} {r[3]:6.3f} {r[4]:9.2f} {r[5]:7d}")
print(f"\nWrote {out}")
PYEOF
