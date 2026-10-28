"""Frames for the QM/MM single points of qmmm_tyr.py: coordinates and force-field charges.

For each replica of one system, FRAMES frames evenly spread over t >= skip are taken from the
centred protein trajectory of 14_analysis.slurm (rep*/md_center.xtc). Atom names and the
original residue numbers come from rep*/protein_ref.gro, charges and elements from
rep*/md_protein.tpr (same atom order; the tpr numbers residues from 1). The trajectories hold
the protein only, so the MM environment is the protein; with full-system frames
(--top/--tpr/--xtc pointing at them) waters and ions within --solvent-cutoff of the Tyr
cluster are kept as well (they must be whole and centred around the protein).

    python qmmm_prep.py --runs runs/dopa --system wt --nrep 4 --frames 10 --skip-ns 50 \
        --sites 65 108 112 124 127 --out runs/dopa/qmmm/wt

Writes OUT/frames/<system>_rep<r>_<time>ps.npz and OUT/frames.tsv (one line per frame, the
index used by the array job).
"""
import argparse
import warnings
from pathlib import Path

import numpy as np
import MDAnalysis as mda

warnings.filterwarnings("ignore")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", required=True)
    ap.add_argument("--system", default="wt")
    ap.add_argument("--nrep", type=int, default=4)
    ap.add_argument("--frames", type=int, default=10, help="frames per replica")
    ap.add_argument("--skip-ns", type=float, default=50.0)
    ap.add_argument("--sites", type=int, nargs="+", default=[65, 108, 112, 124, 127])
    ap.add_argument("--top", default="protein_ref.gro", help="per replica: names and residue numbers")
    ap.add_argument("--tpr", default="md_protein.tpr", help="per replica: charges and elements")
    ap.add_argument("--xtc", default="md_center.xtc")
    ap.add_argument("--solvent-cutoff", type=float, default=20.0,
                    help="keep non-protein atoms within this distance (A) of a site ring")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    out = Path(a.out)
    (out / "frames").mkdir(parents=True, exist_ok=True)
    rows = ["index\tfile\tsystem\trep\ttime_ps\tn_atoms\tn_solvent"]
    k = 0
    for r in range(1, a.nrep + 1):
        rep = Path(a.runs) / a.system / f"rep{r}"
        if not (rep / a.xtc).exists():
            print(f"skip {a.system}/rep{r}: no {a.xtc}")
            continue
        u = mda.Universe(str(rep / a.top), str(rep / a.xtc))
        tpr = mda.Universe(str(rep / a.tpr))
        if len(tpr.atoms) != len(u.atoms):
            raise SystemExit(f"{rep}: {a.tpr} has {len(tpr.atoms)} atoms, {a.top} {len(u.atoms)}")
        charges = tpr.atoms.charges.astype(float)
        elements = np.array([e.capitalize() for e in tpr.atoms.elements])
        protein = u.select_atoms("protein").indices
        rings = u.select_atoms("resid " + " ".join(map(str, a.sites)) + " and name CG CD1 CD2 CE1 CE2 CZ")
        if len(rings) != 6 * len(a.sites):
            raise SystemExit(f"{rep}: expected {6 * len(a.sites)} ring atoms, found {len(rings)}")

        times = np.array([ts.time for ts in u.trajectory])
        pick = np.flatnonzero(times >= a.skip_ns * 1000.0)
        if len(pick) < a.frames:
            raise SystemExit(f"{rep}: only {len(pick)} frames after {a.skip_ns} ns")
        pick = pick[np.linspace(0, len(pick) - 1, a.frames).round().astype(int)]

        for i in pick:
            u.trajectory[i]
            keep = np.zeros(len(u.atoms), bool)
            keep[protein] = True
            other = np.setdiff1d(np.arange(len(u.atoms)), protein)
            if len(other):
                d = np.linalg.norm(u.atoms.positions[other][:, None] - rings.positions[None], axis=-1).min(1)
                # whole residues (water, ion) within the cutoff
                near = set(u.atoms[other[d <= a.solvent_cutoff]].residues.ix)
                keep[[at.ix for res in u.residues if res.ix in near for at in res.atoms]] = True
            sel = np.flatnonzero(keep)
            t = u.trajectory.time
            name = f"{a.system}_rep{r}_{t:09.0f}ps.npz"
            np.savez_compressed(
                out / "frames" / name,
                xyz=u.atoms.positions[sel].astype(np.float64),
                names=u.atoms.names[sel].astype("U4"),
                resnames=u.atoms.resnames[sel].astype("U4"),
                resids=u.atoms.resids[sel].astype(int),
                elements=elements[sel].astype("U2"),
                charges=charges[sel],
                system=a.system, rep=r, time_ps=t, sites=np.array(a.sites))
            rows.append(f"{k}\t{name}\t{a.system}\t{r}\t{t:.0f}\t{len(sel)}\t{len(sel) - len(protein)}")
            k += 1
        print(f"{a.system}/rep{r}: {len(pick)} frames, {times[pick[0]] / 1000:.1f}-{times[pick[-1]] / 1000:.1f} ns")

    (out / "frames.tsv").write_text("\n".join(rows) + "\n", encoding="utf-8")
    print(f"wrote {k} frames to {out / 'frames'}")


if __name__ == "__main__":
    main()
