#!/usr/bin/env python3
"""
Combine n0 and nS arrays into a per-atom dewetting map (eta = nS / n0).
Writes a .dat table and a .pdb with eta in the B-factor column for ChimeraX.

Usage:
    python3 combine_dewetting_map.py --pdbid 1DPX --phi phi_06p0
                                     [--root __SSI_ROOT__]
                                     [--surface-threshold 4.0]
"""
import argparse
from pathlib import Path

import numpy as np
import MDAnalysis as mda


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdbid", required=True)
    ap.add_argument("--phi", required=True,
                    help="phi directory name used for nS, e.g. phi_06p0")
    ap.add_argument("--root", default="__SSI_ROOT__")
    ap.add_argument("--surface-threshold", type=float, default=4.0,
                    help="Atoms with n0 > threshold are surface-exposed (default: 4.0)")
    return ap.parse_args()


def main():
    args    = parse_args()
    root    = Path(args.root)
    pdbid   = args.pdbid.upper()
    phi     = args.phi

    run_dir = root / "runs" / pdbid
    ana_dir = run_dir / "analysis"
    out_dir = run_dir / "dewetting"
    out_dir.mkdir(parents=True, exist_ok=True)

    n0_file = ana_dir / "n0_per_atom.npy"
    nS_file = ana_dir / f"nS_{phi}.npy"
    gro     = run_dir / "unbiased" / "prod_unbiased.gro"

    for f in [n0_file, nS_file, gro]:
        if not f.is_file():
            raise FileNotFoundError(f"Missing: {f}")

    n0 = np.load(str(n0_file))
    nS = np.load(str(nS_file))

    if len(n0) != len(nS):
        raise RuntimeError(f"Length mismatch: n0={len(n0)}, nS={len(nS)}")

    u = mda.Universe(str(gro))
    protein_heavy = u.select_atoms("protein and not name H*")

    if len(protein_heavy) != len(n0):
        raise RuntimeError(
            f"Heavy atom count mismatch: gro={len(protein_heavy)}, arrays={len(n0)}"
        )

    surface = n0 > args.surface_threshold
    eta     = np.full(len(n0), np.nan, dtype=float)
    eta[surface] = nS[surface] / n0[surface]

    valid = eta[~np.isnan(eta)]

    print("=== Dewetting map summary ===")
    print(f"pdbid                   {pdbid}")
    print(f"phi                     {phi}")
    print(f"surface_threshold       {args.surface_threshold}")
    print(f"n_atoms_total           {len(n0)}")
    print(f"n_surface_atoms         {int(surface.sum())}")
    if len(valid) > 0:
        print(f"eta_min                 {np.min(valid):.4f}")
        print(f"eta_p25                 {np.percentile(valid, 25):.4f}")
        print(f"eta_median              {np.median(valid):.4f}")
        print(f"eta_p75                 {np.percentile(valid, 75):.4f}")
        print(f"eta_max                 {np.max(valid):.4f}")
        print(f"n_eta_lt_0.25           {int(np.sum(valid < 0.25))}")
        print(f"n_eta_lt_0.50           {int(np.sum(valid < 0.50))}")
        print(f"n_eta_lt_0.75           {int(np.sum(valid < 0.75))}")

    # --- Write .dat table ---
    dat_file = out_dir / f"dewetting_{phi}.dat"
    with open(dat_file, "w") as fh:
        fh.write("# atom_index_1based resid resname atomname x_nm y_nm z_nm "
                 "n0 nS eta surface\n")
        for i, atom in enumerate(protein_heavy):
            x, y, z = atom.position / 10.0
            eta_str = "nan" if np.isnan(eta[i]) else f"{eta[i]:.6f}"
            fh.write(
                f"{atom.index+1:6d} {atom.resid:6d} {atom.resname:6s} "
                f"{atom.name:6s} "
                f"{x:10.5f} {y:10.5f} {z:10.5f} "
                f"{n0[i]:10.6f} {nS[i]:10.6f} {eta_str:10s} "
                f"{int(surface[i])}\n"
            )
    print(f"\nWrote {dat_file}")

    # --- Write .pdb for ChimeraX ---
    # B-factor = eta for surface atoms, -1 for buried/non-surface
    pdb_file = out_dir / f"dewetting_{phi}.pdb"
    with open(pdb_file, "w") as fh:
        for i, atom in enumerate(protein_heavy):
            x, y, z = atom.position          # already in Angstrom
            b = float(eta[i]) if (surface[i] and not np.isnan(eta[i])) else -1.0
            fh.write(
                "ATOM  {serial:5d} {name:^4s} {resname:>3s} A{resid:4d}    "
                "{x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{b:6.2f}          {elem:>2s}\n".format(
                    serial=i + 1,
                    name=atom.name[:4],
                    resname=atom.resname[:3],
                    resid=int(atom.resid),
                    x=float(x), y=float(y), z=float(z),
                    occ=1.0, b=b,
                    elem=atom.name[0] if atom.name else "C",
                )
            )
        fh.write("END\n")
    print(f"Wrote {pdb_file}")
    print("B-factor: eta for surface atoms, -1 for non-surface.")


if __name__ == "__main__":
    main()
