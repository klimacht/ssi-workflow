#!/usr/bin/env python3
"""
Compute <n_i>_0 (unbiased per-atom hydration) from the unbiased production run.

Usage:
    python3 compute_n0.py --pdbid 1DPX [--root __SSI_ROOT__]
                          [--t-start 1000] [--radius 0.6]
"""
import argparse
import time
from pathlib import Path

import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import capped_distance


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdbid", required=True)
    ap.add_argument("--root", default="__SSI_ROOT__")
    ap.add_argument("--t-start", type=float, default=1000.0,
                    help="Discard frames before this time in ps (default: 1000)")
    ap.add_argument("--radius", type=float, default=0.6,
                    help="Sphere radius in nm (default: 0.6)")
    return ap.parse_args()


def main():
    args = parse_args()
    root     = Path(args.root)
    pdbid    = args.pdbid.upper()
    run_dir  = root / "runs" / pdbid
    ana_dir  = run_dir / "analysis"
    ana_dir.mkdir(parents=True, exist_ok=True)

    gro = run_dir / "unbiased" / "prod_unbiased.gro"
    xtc = run_dir / "unbiased" / "prod_unbiased.xtc"
    out_npy = ana_dir / "n0_per_atom.npy"
    out_dat = ana_dir / "n0_per_atom.dat"

    for f in [gro, xtc]:
        if not f.is_file():
            raise FileNotFoundError(f"Missing: {f}")

    r_ang = args.radius * 10.0

    print(f"Loading universe: {gro}")
    u = mda.Universe(str(gro), str(xtc))

    protein_heavy = u.select_atoms("protein and not name H*")
    water_O       = u.select_atoms("resname SOL and name OW")

    n_atoms  = len(protein_heavy)
    counts   = np.zeros(n_atoms, dtype=np.float64)
    n_frames = 0
    t_start  = time.time()

    print(f"Protein heavy atoms : {n_atoms}")
    print(f"Water oxygens       : {len(water_O)}")
    print(f"Total frames        : {len(u.trajectory)}")
    print(f"t_start             : {args.t_start} ps")
    print(f"radius              : {args.radius} nm")

    # Protein is restrained — use first frame positions as fixed centers
    u.trajectory[0]
    protein_pos = protein_heavy.positions.copy()
    prot_min    = protein_pos.min(axis=0) - r_ang
    prot_max    = protein_pos.max(axis=0) + r_ang

    for ts in u.trajectory:
        if ts.time < args.t_start:
            continue

        water_pos  = water_O.positions
        mask       = np.all(
            (water_pos >= prot_min) & (water_pos <= prot_max), axis=1
        )
        water_near = water_pos[mask]

        if len(water_near) > 0:
            result = capped_distance(
                protein_pos, water_near,
                max_cutoff=r_ang,
                box=ts.dimensions,
                return_distances=False,
            )
            pairs = result[0] if isinstance(result, tuple) else result
            if len(pairs) > 0:
                np.add.at(counts, pairs[:, 0], 1)

        n_frames += 1
        if n_frames % 200 == 0:
            elapsed   = time.time() - t_start
            rate      = n_frames / elapsed
            remaining = (len(u.trajectory) - n_frames) / max(rate, 1e-9)
            print(f"  frame {ts.frame:5d}  t={ts.time:.0f} ps  "
                  f"{n_frames} analyzed  ~{remaining/60:.1f} min left", flush=True)

    if n_frames == 0:
        raise RuntimeError(f"No frames found after t_start={args.t_start} ps")

    n0 = counts / n_frames

    np.save(str(out_npy), n0)
    print(f"\nSaved: {out_npy}")

    with open(out_dat, "w") as fh:
        fh.write("# <n_i>_0 per protein heavy atom\n")
        fh.write(f"# frames={n_frames}  t_start={args.t_start} ps  r={args.radius} nm\n")
        fh.write("# atom_index_1based resid resname atomname n0\n")
        u.trajectory[-1]
        for i, atom in enumerate(protein_heavy):
            fh.write(f"{atom.index+1:6d} {atom.resid:6d} {atom.resname:6s} "
                     f"{atom.name:6s} {n0[i]:10.6f}\n")

    print(f"Saved: {out_dat}")
    print(f"\nSum <n_i>_0 = {n0.sum():.1f}  (= N0)")
    print(f"Atoms with n0 > 4  : {(n0 > 4).sum()} / {n_atoms}  (surface atoms)")
    print(f"Frames analyzed    : {n_frames}")
    print(f"Time               : {(time.time()-t_start)/60:.1f} min")


if __name__ == "__main__":
    main()
