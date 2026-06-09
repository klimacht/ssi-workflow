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
module load python-waterboa/2025.06
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
                if len(p) > 9 and p[9] == "1" and p[8] != "nan":
                    vals.append(float(p[8]))
        if vals:
            v = np.array(vals); n_surf = len(v)
            contrast = float(np.percentile(v,95) - np.percentile(v,5))

    # global half-dewetting phi: prefer sidecar .center file (computed by
    # combine_dewetting_phistar_map with the same cubic interp), else
    # recompute from N_vs_phi here.
    phi_star_m = float("nan")
    center_file = root / "runs" / pdbid / "dewetting" / "dewetting_phistar.center"
    if center_file.is_file():
        meta = {}
        with open(center_file) as f:
            for line in f:
                k, v = line.strip().split("\t", 1)
                meta[k] = v
        if meta.get("center_mode") == "global_half":
            try: phi_star_m = float(meta.get("center_value", "nan"))
            except ValueError: pass
    if phi_star_m != phi_star_m:  # NaN -> fallback: linear interp from N_vs_phi
        phis, Ns = [], []
        with open(tsv) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                phis.append(float(r["phi"])); Ns.append(float(r["N_mean"]))
        if Ns:
            target = 0.5 * Ns[0]
            for k in range(len(phis)-1):
                if Ns[k] >= target >= Ns[k+1] and Ns[k] != Ns[k+1]:
                    t = (Ns[k] - target) / (Ns[k] - Ns[k+1])
                    phi_star_m = phis[k] + t * (phis[k+1] - phis[k])
                    break

    rows.append((pdbid, n0, phistar_g, gamma, phi_star_m, contrast, n_surf))

rows.sort(key=lambda r: (-(r[3] if r[3]==r[3] else -1)))  # by gamma desc

hdr = ["protein","N0","phi_star","gamma","phi_star_m","patch_contrast_kJmol","n_surface"]
with open(out,"w") as f:
    f.write("\t".join(hdr)+"\n")
    for r in rows:
        f.write(f"{r[0]}\t{r[1]:.1f}\t{r[2]:.1f}\t{r[3]:.3f}\t{r[4]:.3f}\t{r[5]:.2f}\t{r[6]}\n")

print(f"{'protein':8} {'N0':>7} {'phi*':>5} {'gamma':>6} {'phi_star_m':>10} {'contrast':>9} {'n_surf':>7}")
for r in rows:
    print(f"{r[0]:8} {r[1]:7.0f} {r[2]:5.1f} {r[3]:6.3f} {r[4]:9.3f} {r[5]:9.2f} {r[6]:7d}")
print(f"\nWrote {out}")
PYEOF
