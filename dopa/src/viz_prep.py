#!/usr/bin/env python3
"""Structures and movie frames for the PyMOL renderings (render_dopa_pymol.py).

Reads the replicas (RUNS/<system>/rep<i>/protein_ref.gro + md_center.xtc, written by
14_analysis.slurm) and the analysis results (ANALYSIS/site*/per_residue.tsv,
ANALYSIS/tyr_exposure/wt/exposure.tsv and frames_rep*.tsv) and writes to OUT:

  rep_<system>.pdb        representative frame: the frame closest to the mean structure
                          (core Calpha + Tyr rings of the cluster) over all replicas after
                          --skip-ns, fitted on the core Calpha to the WT reference (rep1 start);
                          B = Calpha RMSF of the residue (A), occupancy = RMSF - RMSF(wt)
  exposure_paper.pdb      rep_wt; reactive Tyr atoms: B = ASA of Table 1 (A^2)
  exposure_md.pdb         rep_wt; reactive Tyr atoms: B = MD mean ASA, occupancy = fraction
                          of time accessible (> 1 A^2)
  movie_<name>.pdb/.xtc   one WT replica window, fitted, coordinates averaged over --smooth-ns
  movie_<name>.tsv        per movie frame: time, ASA of the reactive atoms and Tyr-Tyr distances
                          (tyr_exposure frames, same smoothing), closest ortho C pair of --pair
  viz_info.json           systems, modified residues, movie windows (read by render and compose)

Movies (WT, chosen automatically unless --exposure-window / --pair-window):
  exposure  the largest rise of the phenol-O ASA of --exposure-tyr (default 108: buried in
            Table 1, strongly nitrated), placed in the first third of the window
  pair      the window with the widest range of the ortho...ortho distance of --pair that
            reaches the coupling cutoff (default 65-127, cross-linked in the paper)
"""
import argparse
import json
from pathlib import Path

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.align import rotation_matrix

import tyr_exposure as T

CLUSTER = (65, 108, 112, 124, 127)
RING = "CG CD1 CD2 CE1 CE2 CZ"
MODIFIED = ("DAH", "TOP", "DTY")
ORTHO = ("CE1", "CE2")


def universe(runs, system, r):
    d = Path(runs) / system / f"rep{r}"
    return mda.Universe(str(d / "protein_ref.gro"), str(d / "md_center.xtc"))


def fit(mobile, ref_centered):
    """(rotation, centre) that put `mobile` onto the centred reference: x -> (x - c) @ R.T."""
    c = mobile.mean(0)
    R, _ = rotation_matrix(mobile - c, ref_centered)
    return R, c


def read_rmsf(analysis):
    """{system: {resid: RMSF}} from all site*/per_residue.tsv (first value per system wins)."""
    rmsf = {}
    for f in sorted(Path(analysis).glob("site*/per_residue.tsv")):
        for row in T.read_tsv(f):
            sys_ = rmsf.setdefault(row["system"], {})
            sys_.setdefault(int(row["resid"]), float(row["rmsf_A_mean"]))
    return rmsf


def write_pdb(u, path, positions, bfac, occ):
    u.atoms.positions = positions
    u.atoms.tempfactors = np.nan_to_num(bfac)
    u.atoms.occupancies = np.nan_to_num(occ)
    u.atoms.write(str(path))


def prepare_universe(u):
    n = len(u.atoms)
    u.add_TopologyAttr("tempfactors", np.zeros(n))
    u.add_TopologyAttr("occupancies", np.zeros(n))
    for attr, value in (("chainIDs", ["A"] * n), ("altLocs", [""] * n),
                        ("elements", [nm.lstrip("0123456789")[0] for nm in u.atoms.names])):
        if not hasattr(u.atoms, attr):
            u.add_TopologyAttr(attr, value)
    return u


def representative(args, system, ref_core):
    sel = f"(name CA and resid {args.core}) or (resid {' '.join(map(str, CLUSTER))} and name {RING})"
    X, where = [], []
    for r in range(1, args.nrep + 1):
        if not (Path(args.runs) / system / f"rep{r}" / "md_center.xtc").exists():
            continue
        u = universe(args.runs, system, r)
        core = u.select_atoms(f"name CA and resid {args.core}")
        med = u.select_atoms(sel)
        step = max(1, round(args.every_ns * 1000.0 / u.trajectory.dt))
        t_first = u.trajectory[0].time
        for ts in u.trajectory[::step]:
            t_ns = (ts.time - t_first) / 1000.0
            if t_ns < args.skip_ns:
                continue
            R, c = fit(core.positions, ref_core)
            X.append((med.positions - c) @ R.T)
            where.append((r, ts.frame, t_ns))
    if not X:
        return None
    X = np.array(X)
    dev = np.sqrt(((X - X.mean(0)) ** 2).sum(-1).mean(-1))
    k = int(dev.argmin())
    r, frame, t_ns = where[k]
    u = prepare_universe(universe(args.runs, system, r))
    u.trajectory[frame]
    core = u.select_atoms(f"name CA and resid {args.core}")
    R, c = fit(core.positions, ref_core)
    return u, (u.atoms.positions - c) @ R.T, {"rep": r, "frame": frame, "time_ns": round(t_ns, 3),
                                             "rmsd_to_mean_A": round(float(dev[k]), 3), "n_frames": len(X)}


def per_atom(u, values_by_resid):
    return np.array([values_by_resid.get(int(a.resid), np.nan) for a in u.atoms])


def smooth(y, n):
    """Centred running mean over n points (edges: shorter windows)."""
    if n <= 1:
        return y
    c = np.cumsum(np.insert(np.asarray(y, float), 0, 0.0, axis=0), axis=0)
    h = n // 2
    i = np.arange(len(y))
    lo, hi = np.clip(i - h, 0, len(y)), np.clip(i + h + 1, 0, len(y))
    return (c[hi] - c[lo]) / (hi - lo)[:, None] if np.ndim(y) > 1 else (c[hi] - c[lo]) / (hi - lo)


def wt_frames(args):
    out = []
    for r in range(1, args.nrep + 1):
        f = Path(args.analysis) / "tyr_exposure" / "wt" / f"frames_rep{r}.tsv"
        if f.exists():
            times, cols, data = T.read_frames(f)
            out.append((r, times, cols, data))
    if not out:
        raise SystemExit("no tyr_exposure/wt/frames_rep*.tsv: run tyr_exposure.py md first")
    return out


def pick_exposure(frames, args):
    col = f"Y{args.exposure_tyr}:O"
    best = None
    for r, t, cols, data in frames:
        dt = t[1] - t[0]
        y = smooth(data[:, cols.index(col)], max(1, round(1.0 / dt)))
        lag = max(1, round(5.0 / dt))
        ok = np.where(t[:-lag] >= args.skip_ns)[0]
        if len(ok) == 0:
            continue
        rise = y[ok + lag] - y[ok]
        k = ok[int(rise.argmax())]
        if best is None or rise.max() > best[0]:
            best = (float(rise.max()), r, float(t[k]), float(t[-1]))
    _, r, t_event, t_end = best
    t0 = min(max(args.skip_ns, t_event - args.movie_ns / 3.0), t_end - args.movie_ns)
    return r, t0


def pick_pair(frames, args):
    """Window with the widest range of the ortho...ortho distance that reaches the cutoff."""
    a, b = args.pair.split("-")
    col = f"{a}-{b}:oo"
    best = None
    for r, t, cols, data in frames:
        dt = t[1] - t[0]
        y = smooth(data[:, cols.index(col)], max(1, round(0.5 / dt)))
        n = max(1, round(args.movie_ns / dt))
        step = max(1, round(1.0 / dt))
        for i in range(0, len(t) - n, step):
            if t[i] < args.skip_ns:
                continue
            w = y[i:i + n]
            score = (bool(w.min() <= T.CONTACT), float(w.max() - w.min()))
            if best is None or score > best[0]:
                best = (score, r, float(t[i]))
    return best[1], best[2]


def movie(args, name, kind, r, t0, ref_core, frames, info):
    u = prepare_universe(universe(args.runs, "wt", r))
    core = u.select_atoms(f"name CA and resid {args.core}")
    dt_ns = u.trajectory.dt / 1000.0
    t_first = u.trajectory[0].time / 1000.0
    n_out = int(round(args.movie_s * args.fps))
    t_out = t0 + np.arange(n_out) * args.movie_ns / n_out
    half = args.smooth_ns / 2.0
    i0 = max(0, int(np.floor((t_out[0] - half) / dt_ns)))         # times relative to the first frame
    i1 = min(len(u.trajectory) - 1, int(np.ceil((t_out[-1] + half) / dt_ns)))
    src_t, src_x = [], []
    for ts in u.trajectory[i0:i1 + 1]:
        R, c = fit(core.positions, ref_core)
        src_t.append(ts.time / 1000.0 - t_first)
        src_x.append((u.atoms.positions - c) @ R.T)
    src_t, src_x = np.array(src_t), np.array(src_x, dtype=np.float64)
    csum = np.concatenate([np.zeros((1,) + src_x.shape[1:]), np.cumsum(src_x, axis=0)])
    lo = np.searchsorted(src_t, t_out - half, side="left")
    hi = np.maximum(np.searchsorted(src_t, t_out + half, side="right"), lo + 1)
    xyz = (csum[hi] - csum[lo]) / (hi - lo)[:, None, None]

    out = Path(args.out)
    write_pdb(u, out / f"movie_{name}.pdb", xyz[0], np.zeros(len(u.atoms)), np.ones(len(u.atoms)))
    with mda.Writer(str(out / f"movie_{name}.xtc"), len(u.atoms)) as w:
        for k in range(n_out):
            u.atoms.positions = xyz[k]
            w.write(u.atoms)

    # per-frame series from tyr_exposure (same box smoothing), closest ortho pair from the coordinates
    _, t, cols, data = next(f for f in frames if f[0] == r)
    dt = t[1] - t[0]
    ys = smooth(data, max(1, round(args.smooth_ns / dt)))
    series = {c: np.interp(t_out, t, ys[:, j]) for j, c in enumerate(cols)}
    pa, pb = (int(x) for x in args.pair.split("-"))
    names_a = [u.select_atoms(f"resid {pa} and name {n}")[0].index for n in ORTHO]
    names_b = [u.select_atoms(f"resid {pb} and name {n}")[0].index for n in ORTHO]
    pairs = [(i, j) for i in names_a for j in names_b]
    d = np.array([[np.linalg.norm(x[i] - x[j]) for i, j in pairs] for x in xyz])
    best = d.argmin(1)
    with open(out / f"movie_{name}.tsv", "w") as fh:
        fh.write("time_ns\t" + "\t".join(cols) + "\tpair_a\tpair_b\tpair_d\n")
        for k in range(n_out):
            i, j = pairs[best[k]]
            fh.write(f"{t_out[k]:.3f}\t" + "\t".join(f"{series[c][k]:.3f}" for c in cols)
                     + f"\t{u.atoms[i].name}\t{u.atoms[j].name}\t{d[k, best[k]]:.3f}\n")
    info["movies"][name] = {"kind": kind, "system": "wt", "rep": r, "t0_ns": round(t0, 3),
                            "length_ns": args.movie_ns, "n_frames": n_out, "fps": args.fps,
                            "smooth_ns": args.smooth_ns, "exposure_tyr": args.exposure_tyr,
                            "pair": [pa, pb], "cutoff_A": T.CONTACT}
    print(f"movie {name}: wt rep{r}, {t0:.1f}-{t0 + args.movie_ns:.1f} ns, {n_out} frames "
          f"from {len(src_t)} trajectory frames")


def window(spec):
    r, t0 = spec.split(":")
    return int(r), float(t0)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", required=True)
    ap.add_argument("--analysis", required=True, help="directory with site*/ and tyr_exposure/")
    ap.add_argument("--systems", nargs="+", required=True)
    ap.add_argument("--nrep", type=int, required=True)
    ap.add_argument("--skip-ns", type=float, default=50.0)
    ap.add_argument("--core", default="37:140", help="residues for fitting")
    ap.add_argument("--every-ns", type=float, default=0.1, help="frame spacing for the representative frame")
    ap.add_argument("--exposure-tyr", type=int, default=108)
    ap.add_argument("--pair", default="65-127")
    ap.add_argument("--exposure-window", type=window, help="REP:T0_NS instead of the automatic choice")
    ap.add_argument("--pair-window", type=window, help="REP:T0_NS instead of the automatic choice")
    ap.add_argument("--movie-ns", type=float, default=25.0, help="simulated time per movie")
    ap.add_argument("--movie-s", type=float, default=20.0, help="movie length in seconds")
    ap.add_argument("--fps", type=int, default=25)
    ap.add_argument("--smooth-ns", type=float, default=0.25, help="coordinate averaging window")
    ap.add_argument("--no-movies", action="store_true")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    if "wt" not in args.systems:
        raise SystemExit("the wild type is the reference: include wt in --systems")
    ref_u = universe(args.runs, "wt", 1)
    ref_core = ref_u.select_atoms(f"name CA and resid {args.core}").positions
    ref_core = ref_core - ref_core.mean(0)
    rmsf = read_rmsf(args.analysis)
    info = {"skip_ns": args.skip_ns, "cluster": list(CLUSTER), "systems": {}, "movies": {}}

    reps = {}
    for system in args.systems:
        res = representative(args, system, ref_core)
        if res is None:
            print(f"skip {system}: no trajectories")
            continue
        u, xyz, meta = res
        mod = sorted({(int(a.resid), a.resname) for a in u.atoms if a.resname in MODIFIED})
        r_sys = per_atom(u, rmsf.get(system, {}))
        r_wt = per_atom(u, rmsf.get("wt", {}))
        write_pdb(u, out / f"rep_{system}.pdb", xyz, r_sys, r_sys - r_wt if system != "wt" else 0 * r_sys)
        info["systems"][system] = {"modified": [list(m) for m in mod], "rmsf": system in rmsf, **meta}
        reps[system] = (u, xyz)
        print(f"{system}: rep{meta['rep']} at {meta['time_ns']:.1f} ns, {meta['rmsd_to_mean_A']:.2f} A "
              f"from the mean of {meta['n_frames']} frames; modified {mod or '-'}")

    # exposure maps on the representative WT structure
    u, xyz = reps["wt"]
    rows = T.read_tsv(Path(args.analysis) / "tyr_exposure" / "wt" / "exposure.tsv")
    val = {(int(r["tyr"]), r["quantity"]): r for r in rows}
    paper, md, acc = (np.zeros(len(u.atoms)) for _ in range(3))
    for k, a in enumerate(u.atoms):
        r = val.get((int(a.resid), T.PAPER_NAME.get(a.name, "")))
        if a.resname == "TYR" and a.name in T.REACTIVE and r:
            paper[k], md[k], acc[k] = T.num(r["paper"]), T.num(r["md_mean"]), T.num(r["frac_accessible"])
    write_pdb(u, out / "exposure_paper.pdb", xyz, paper, np.zeros(len(u.atoms)))
    write_pdb(u, out / "exposure_md.pdb", xyz, md, acc)

    if not args.no_movies:
        frames = wt_frames(args)
        r, t0 = args.exposure_window or pick_exposure(frames, args)
        movie(args, "exposure", "exposure", r, t0, ref_core, frames, info)
        r, t0 = args.pair_window or pick_pair(frames, args)
        movie(args, "pair", "pair", r, t0, ref_core, frames, info)

    (out / "viz_info.json").write_text(json.dumps(info, indent=1))
    print(f"structures in {out}")


if __name__ == "__main__":
    main()
