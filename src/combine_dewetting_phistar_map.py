#!/usr/bin/env python3
"""
combine_dewetting_phistar_map.py

Per-atom dewetting map that uses ALL phi windows instead of a single phi*.

For each surface atom i we have the hydration curve <n_i>_phi over all phi
values. We define the atom-wise dewetting transition point phi*_i as the phi
at which the atom has lost half of its unbiased hydration:

    <n_i>_{phi*_i} = 0.5 * <n_i>_0

phi*_i is obtained by linear interpolation between the two bracketing windows,
giving a continuous value (e.g. 4.7 kJ/mol) rather than a grid point.

Interpretation:
    low  phi*_i  -> atom dewets early  -> hydrophobic  (red)
    high phi*_i  -> atom holds water   -> hydrophilic  (blue)

This uses every simulation that was run, is continuous, and is far less
sensitive to noise in any single window than the single-phi* eta map.

Outputs (in runs/PDBID/dewetting/):
    dewetting_phistar.dat   per-atom table (n0, phistar_i, all <n_i>_phi)
    dewetting_phistar.pdb   B-factor = phi*_i for surface atoms, -1 otherwise

Usage:
    python3 combine_dewetting_phistar_map.py --pdbid 1DPX
        [--root __SSI_ROOT__] [--surface-threshold 4.0]
        [--phi-list "0 1 2 3 4 5 6 7 8 9 10 11 12"]
"""
import argparse
from pathlib import Path

import numpy as np
import MDAnalysis as mda


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdbid", required=True)
    ap.add_argument("--root", default="__SSI_ROOT__")
    ap.add_argument("--surface-threshold", type=float, default=4.0,
                    help="Atoms with n0 > threshold are surface-exposed (default 4.0)")
    ap.add_argument("--phi-list", default="0 1 2 3 4 5 6 7 8 9 10 11 12",
                    help="Space-separated phi values to include")
    ap.add_argument("--frac", type=float, default=0.5,
                    help="Hydration fraction defining the transition (default 0.5)")
    ap.add_argument("--center", action="store_true", default=True,
                    help="Subtract a reference so B-factor is delta around 0 (default on)")
    ap.add_argument("--no-center", dest="center", action="store_false",
                    help="Write raw phi*_i in kJ/mol instead of centered delta")
    ap.add_argument("--center-mode", choices=["phistar_global","median"], default="phistar_global",
                    help="Centering reference: phistar_global (susceptibility max) or median of phi*_i")
    return ap.parse_args()


def phi_to_dirname(phi):
    whole = int(phi)
    frac = int(round((phi - whole) * 10))
    return f"phi_{whole:02d}p{frac}"


def interp_phistar(phis, n_curve, n0, frac):
    """
    Find the phi at which n_curve crosses frac*n0, by linear interpolation.

    phis    : 1D array of phi values (sorted ascending)
    n_curve : 1D array of <n_i>_phi, same length
    n0      : unbiased hydration for this atom
    frac    : target fraction (0.5 = half-dehydration)

    Returns phi*_i, or np.nan if the curve never crosses the target
    (e.g. atom already below target at phi=0, or never drops that far).
    """
    target = frac * n0
    # Curve should be (mostly) decreasing. Find first window where it drops
    # from >= target to < target.
    for k in range(len(phis) - 1):
        a, b = n_curve[k], n_curve[k + 1]
        if a >= target >= b and a != b:
            # linear interpolation between phi[k] and phi[k+1]
            t = (a - target) / (a - b)
            return phis[k] + t * (phis[k + 1] - phis[k])
    # No clean crossing
    if n_curve[0] < target:
        return phis[0]            # already dehydrated at phi=0
    return np.nan                 # never reaches target within phi range


def main():
    args = parse_args()
    root = Path(args.root)
    pdbid = args.pdbid.upper()

    run_dir = root / "runs" / pdbid
    ana_dir = run_dir / "analysis"
    out_dir = run_dir / "dewetting"
    out_dir.mkdir(parents=True, exist_ok=True)

    gro = run_dir / "unbiased" / "prod_unbiased.gro"
    n0_file = ana_dir / "n0_per_atom.npy"

    if not n0_file.is_file():
        raise FileNotFoundError(f"Missing {n0_file} (run compute_n0.py first)")
    if not gro.is_file():
        raise FileNotFoundError(f"Missing {gro}")

    n0 = np.load(str(n0_file))
    phis = [float(x) for x in args.phi_list.split()]

    # Load all available nS_phi_*.npy windows.
    n_curves = []          # list of (phi, array)
    for phi in phis:
        dirname = phi_to_dirname(phi)
        f = ana_dir / f"nS_{dirname}.npy"
        if not f.is_file():
            print(f"  WARNING: missing {f.name} — skipping phi={phi}")
            continue
        arr = np.load(str(f))
        if len(arr) != len(n0):
            raise RuntimeError(
                f"Length mismatch in {f.name}: {len(arr)} vs n0 {len(n0)}")
        n_curves.append((phi, arr))

    if len(n_curves) < 3:
        raise RuntimeError(
            f"Only {len(n_curves)} phi windows found; need at least 3 for a "
            f"meaningful per-atom transition. Run compute_nS.py for all phi "
            f"first (see note below).")

    # phi=0 must be present as the reference top of each curve.
    phi_vals = np.array([p for p, _ in n_curves], dtype=float)
    order = np.argsort(phi_vals)
    phi_vals = phi_vals[order]
    stack = np.vstack([n_curves[i][1] for i in order])   # shape (n_phi, n_atoms)

    n_atoms = len(n0)
    surface = n0 > args.surface_threshold

    phistar = np.full(n_atoms, np.nan, dtype=float)
    for i in range(n_atoms):
        if not surface[i]:
            continue
        phistar[i] = interp_phistar(phi_vals, stack[:, i], n0[i], args.frac)

    # Global phi* (collective susceptibility maximum) for centering.
    phi_star_global = None
    tsv = ana_dir / "N_vs_phi.tsv"
    if tsv.is_file():
        import csv as _csv
        with open(tsv) as _f:
            for _row in _csv.DictReader(_f, delimiter="\t"):
                if _row.get("phi_star", "0") == "1":
                    phi_star_global = float(_row["phi"])
                    break

    raw_phistar = phistar.copy()
    valid_init = phistar[surface & ~np.isnan(phistar)]
    median_phistar = float(np.median(valid_init)) if len(valid_init) > 0 else None
    if args.center:
        if args.center_mode == "median" and median_phistar is not None:
            center_value = median_phistar
            center_label = f"median(phi*_i)={median_phistar:.3f}"
        elif phi_star_global is not None:
            center_value = phi_star_global
            center_label = f"phi*_global={phi_star_global:.3f}"
        else:
            center_value = None; center_label = "none"
        if center_value is not None:
            phistar = phistar - center_value

    valid = phistar[surface & ~np.isnan(phistar)]

    if args.center and (phi_star_global is not None or median_phistar is not None):
        mode = f"delta (phi*_i - {center_label})"
    else:
        mode = "raw phi*_i"
    print("=== Per-atom dewetting transition map ===")
    print(f"pdbid                 {pdbid}")
    print(f"phi windows used      {[float(x) for x in phi_vals]}")
    print(f"surface_threshold     {args.surface_threshold}")
    print(f"transition fraction   {args.frac}")
    print(f"global phi*           {phi_star_global}  (collective susceptibility max)")
    print(f"B-factor quantity     {mode}")
    print(f"n_atoms_total         {n_atoms}")
    print(f"n_surface_atoms       {int(surface.sum())}")
    if len(valid) > 0:
        lo_label = "more hydrophobic" if args.center else "most hydrophobic"
        hi_label = "more hydrophilic" if args.center else "most hydrophilic"
        print(f"value min             {valid.min():+.3f}  ({lo_label})")
        print(f"value p05             {np.percentile(valid, 5):+.3f}")
        print(f"value p25             {np.percentile(valid, 25):+.3f}")
        print(f"value median          {np.median(valid):+.3f}")
        print(f"value p75             {np.percentile(valid, 75):+.3f}")
        print(f"value p95             {np.percentile(valid, 95):+.3f}")
        print(f"value max             {valid.max():+.3f}  ({hi_label})")
        # Patch-contrast metric: how strongly differentiated is the surface?
        spread_p = float(np.percentile(valid, 95) - np.percentile(valid, 5))
        spread_sd = float(np.std(valid))
        print(f"patch contrast (p95-p05) {spread_p:.3f} kJ/mol")
        print(f"patch contrast (std)     {spread_sd:.3f} kJ/mol")
        n_nan = int(surface.sum() - len(valid))
        print(f"surface atoms w/o crossing: {n_nan} "
              f"(never reached {args.frac:.0%} dehydration)")

    u = mda.Universe(str(gro))
    protein_heavy = u.select_atoms("protein and not name H*")
    if len(protein_heavy) != n_atoms:
        raise RuntimeError(
            f"Atom count mismatch: gro {len(protein_heavy)} vs arrays {n_atoms}")

    # --- table ---
    dat = out_dir / "dewetting_phistar.dat"
    phi_hdr = " ".join(f"n_phi{p:g}" for p in phi_vals)
    with open(dat, "w") as fh:
        fh.write(f"# atom_index_1based resid resname atomname x_nm y_nm z_nm "
                 f"n0 phistar surface {phi_hdr}\n")
        for i, atom in enumerate(protein_heavy):
            x, y, z = atom.position / 10.0
            ps = "nan" if np.isnan(phistar[i]) else f"{phistar[i]:.4f}"
            curve = " ".join(f"{stack[k, i]:.4f}" for k in range(len(phi_vals)))
            fh.write(f"{atom.index+1:6d} {atom.resid:6d} {atom.resname:6s} "
                     f"{atom.name:6s} {x:8.3f} {y:8.3f} {z:8.3f} "
                     f"{n0[i]:8.4f} {ps:>8s} {int(surface[i])} {curve}\n")
    print(f"\nWrote {dat}")

    # --- pdb: B-factor = (centered) phi*_i for surface atoms ---
    # Buried atoms get a sentinel well below the data range so they are easy
    # to hide ("beta < BURIED+1") and never collide with real (possibly
    # negative) centered values.
    # Atoms that never cross the threshold are the most hydrophilic; fill them
    # with the maximum observed surface value so they color at that end.
    if len(valid) > 0:
        fill_hydrophilic = float(valid.max())
        buried_sentinel = float(np.floor(valid.min() - 5.0))
    else:
        fill_hydrophilic = 0.0
        buried_sentinel = -99.0

    pdb = out_dir / "dewetting_phistar.pdb"
    with open(pdb, "w") as fh:
        for i, atom in enumerate(protein_heavy):
            x, y, z = atom.position
            if surface[i]:
                b = phistar[i] if not np.isnan(phistar[i]) else fill_hydrophilic
            else:
                b = buried_sentinel
            fh.write(
                "ATOM  {serial:5d} {name:^4s} {resname:>3s} A{resid:4d}    "
                "{x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{b:6.2f}          {elem:>2s}\n".format(
                    serial=i + 1, name=atom.name[:4], resname=atom.resname[:3],
                    resid=int(atom.resid), x=float(x), y=float(y), z=float(z),
                    occ=1.0, b=float(b),
                    elem=atom.name[0] if atom.name else "C"))
        fh.write("END\n")
    print(f"Wrote {pdb}")

    if args.center and (phi_star_global is not None or median_phistar is not None):
        print(f"\nB-factor = {center_label} (kJ/mol), centered on 0:")
        print(f"  negative (RED)   = dewets before the collective = more hydrophobic")
        print(f"  ~0     (WHITE)   = behaves like the protein average")
        print(f"  positive (BLUE)  = holds water longer = more hydrophilic")
        if len(valid) > 0:
            m = max(abs(valid.min()), abs(valid.max()))
            print(f"  symmetric color range: -{m:.1f} to +{m:.1f}")
            print(f"  buried-atom sentinel:  {buried_sentinel:.1f} (hide with 'beta > {buried_sentinel+1:.0f}')")
        print(f"\nVMD:  color scale method BWR   (B=neg→blue? use RWB if reversed)")
        print(f"      mol modselect 0 top \"beta > {buried_sentinel+1:.0f}\"")
        print(f"      mol modcolor 0 top Beta ; mol scaleminmax top 0 -{m:.1f} {m:.1f}")
    else:
        print(f"\nB-factor = raw phi*_i (kJ/mol). low=hydrophobic(RED), high=hydrophilic(BLUE)")
        if len(valid) > 0:
            print(f"  suggested color range: {valid.min():.1f} to {valid.max():.1f}")
        print(f"  buried-atom sentinel: {buried_sentinel:.1f}")


if __name__ == "__main__":
    main()
