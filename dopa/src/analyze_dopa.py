#!/usr/bin/env python3
"""Compare WT and DOPA variants of a protein domain over independent replicas.

Expects per replica (written by dopa/slurm/14_analysis.slurm):
  RUNS/<system>/rep<i>/protein_ref.gro        protein at the start of production
  RUNS/<system>/rep<i>/md_center.xtc          protein trajectory, whole + centered
  RUNS/<system>/rep<i>/sasa_res<MUT>.xvg      side-chain SASA of the mutation site

Per replica (after discarding --skip-ns):
  global   Calpha RMSD of the helical core, RMSF per residue, helicity (DSSP)
  site     chi1/chi2 of residue MUT, ring flips (|delta chi2| > 90 deg from start),
           ring-centroid distances to aromatic neighbours, side-chain contacts,
           hydrogen bonds of the hydroxyl groups, side-chain SASA
Systems are compared by mean +- SEM over replicas; p-values (Welch t-test on
replica means, n = number of replicas) are a rough guide only.

Outputs in --out: summary.tsv, per_residue.tsv, contacts.tsv, hbonds.tsv, *.png
"""
import argparse
from collections import defaultdict
from pathlib import Path

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.align import rotation_matrix
from MDAnalysis.lib.distances import distance_array, calc_dihedrals

RING = ("CG", "CD1", "CD2", "CE1", "CE2", "CZ")
SITE_OXYGENS = ("OH", "OE2", "OD1")       # Tyr / DOPA / TOPA hydroxyl oxygens
SITE_HYDROGENS = ("HH", "HOE2", "HOD1")
CONTACT_CUTOFF = 4.5      # A, heavy atoms
HB_DIST, HB_ANGLE = 3.5, 150.0


def read_xvg(path):
    rows = [l.split() for l in open(path) if l.strip() and l[0] not in "#@"]
    return np.array(rows, dtype=float)


def polar_hydrogens(u):
    """(H atom, donor atom) pairs: hydrogens within 1.2 A of an N or O in the first frame."""
    hs = u.select_atoms("name H*")
    heavy = u.select_atoms("name N* O*")
    d = distance_array(hs.positions, heavy.positions)
    pairs = []
    for i, h in enumerate(hs):
        j = int(d[i].argmin())
        if d[i, j] < 1.2:
            pairs.append((h.index, heavy[j].index))
    return pairs


def angle_deg(a, b, c):
    v1, v2 = a - b, c - b
    cosang = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    return np.degrees(np.arccos(np.clip(cosang, -1.0, 1.0)))


def analyze_replica(rep_dir, args):
    u = mda.Universe(str(rep_dir / "protein_ref.gro"), str(rep_dir / "md_center.xtc"))
    site = u.select_atoms(f"resid {args.mut}")
    resname = site.residues.resnames[0] if len(site) else "?"
    if resname not in ("TYR", "DAH", "TOP", "DTY"):
        raise SystemExit(f"{rep_dir}: residue {args.mut} is {resname}; residue numbering not preserved?")

    core_sel = f"name CA and resid {args.core}"
    core = u.select_atoms(core_sel)
    ca = u.select_atoms(f"name CA and resid {args.first}:{args.last}")
    ref_core = core.positions - core.positions.mean(0)

    chi1_atoms = [site.select_atoms(f"name {n}")[0] for n in ("N", "CA", "CB", "CG")]
    chi2_atoms = [site.select_atoms(f"name {n}")[0] for n in ("CA", "CB", "CG", "CD1")]
    ring_site = site.select_atoms("name " + " ".join(RING))
    partners = {}
    for p in args.partners:
        ring = u.select_atoms(f"resid {p} and name " + " ".join(RING))
        if len(ring) == 6:
            partners[p] = ring

    sc_heavy = site.select_atoms("not name N CA C O and not name H*")
    env_heavy = u.select_atoms(f"not resid {args.mut - 1}:{args.mut + 1} and not name H*")
    env_labels = np.array([f"{a.resname}{a.resid}" for a in env_heavy])
    chi2_start = np.degrees(calc_dihedrals(*[a.position for a in chi2_atoms]))   # frame 0

    hpairs = polar_hydrogens(u)
    site_idx = set(site.indices)
    site_donors = [(h, d) for h, d in hpairs if h in site_idx and u.atoms[h].name in SITE_HYDROGENS]
    other_donors = [(h, d) for h, d in hpairs if h not in site_idx]
    site_acc = site.select_atoms("name " + " ".join(SITE_OXYGENS))
    other_acc = u.select_atoms(f"(name O* or (resname HIS HID HIE HIP and name ND1 NE2)) and not resid {args.mut}")

    t_first = u.trajectory[0].time
    t0 = args.skip_ns * 1000.0
    times, rmsd, ca_xyz, chi1, chi2 = [], [], [], [], []
    cent = defaultdict(list)
    contact_counts = defaultdict(int)
    hb_counts = defaultdict(int)
    n_frames = 0
    frame_ids = []
    for ts in u.trajectory[::args.stride]:
        t_rel = ts.time - t_first
        if t_rel < t0:
            continue
        n_frames += 1
        frame_ids.append(ts.frame)
        times.append(t_rel / 1000.0)
        mob = core.positions
        com = mob.mean(0)
        R, r = rotation_matrix(mob - com, ref_core)
        rmsd.append(r)
        ca_xyz.append((ca.positions - com) @ R.T)

        chi1.append(np.degrees(calc_dihedrals(*[a.position for a in chi1_atoms])))
        chi2.append(np.degrees(calc_dihedrals(*[a.position for a in chi2_atoms])))
        c_site = ring_site.positions.mean(0)
        for p, ring in partners.items():
            cent[p].append(np.linalg.norm(ring.positions.mean(0) - c_site))

        dmin = distance_array(sc_heavy.positions, env_heavy.positions).min(axis=0)
        for label in np.unique(env_labels[dmin < CONTACT_CUTOFF]):
            contact_counts[str(label)] += 1

        pos = u.atoms.positions
        seen = set()
        for h, d in site_donors:                      # site OH as donor
            dist = np.linalg.norm(other_acc.positions - pos[d], axis=1)
            for k in np.where(dist < HB_DIST)[0]:
                acc = other_acc[k]
                if angle_deg(pos[d], pos[h], acc.position) >= HB_ANGLE:
                    seen.add(f"{u.atoms[h].name}->{acc.resname}{acc.resid}:{acc.name}")
        for acc in site_acc:                          # site O as acceptor
            for h, d in other_donors:
                if np.linalg.norm(pos[d] - acc.position) < HB_DIST and \
                        angle_deg(pos[d], pos[h], acc.position) >= HB_ANGLE:
                    da = u.atoms[d]
                    seen.add(f"{da.resname}{da.resid}:{da.name}->{acc.name}")
        for key in seen:
            hb_counts[key] += 1

    if n_frames == 0:
        raise SystemExit(f"{rep_dir}: no frames after {args.skip_ns} ns")

    ca_xyz = np.array(ca_xyz)
    mean = ca_xyz.mean(0)
    for _ in range(2):                                # refine the average structure
        mean_c = mean - mean.mean(0)
        for i in range(len(ca_xyz)):
            x = ca_xyz[i] - ca_xyz[i].mean(0)
            R, _ = rotation_matrix(x, mean_c)
            ca_xyz[i] = x @ R.T
        mean = ca_xyz.mean(0)
    rmsf = np.sqrt(((ca_xyz - mean) ** 2).sum(-1).mean(0))

    helicity = None
    try:
        from MDAnalysis.analysis.dssp import DSSP
        dssp_res = u.select_atoms(f"resid {args.first}:{args.last}")
        dssp = DSSP(dssp_res).run(frames=frame_ids)
        helicity = (np.asarray(dssp.results.dssp) == "H").mean(axis=0)
    except Exception as exc:                          # older MDAnalysis or missing atoms
        print(f"   DSSP skipped for {rep_dir}: {exc}")

    chi1, chi2 = np.array(chi1), np.array(chi2)
    dchi2 = (chi2 - chi2_start + 180.0) % 360.0 - 180.0
    sasa_file = rep_dir / f"sasa_res{args.mut}.xvg"
    sasa = None
    if sasa_file.exists():
        s = read_xvg(sasa_file)
        sasa = s[s[:, 0] - s[0, 0] >= t0, 1]

    return {
        "resname": resname,
        "time": np.array(times),
        "rmsd": np.array(rmsd),
        "rmsf_resids": ca.resids.copy(),
        "rmsf": rmsf,
        "helicity": helicity,
        "chi1": chi1,
        "chi2": chi2,
        "scalars": {
            "core_rmsd_A": float(np.mean(rmsd)),
            "chi1_frac_m(-120..0)": float(np.mean((chi1 >= -120) & (chi1 < 0))),
            "chi1_frac_p(0..120)": float(np.mean((chi1 >= 0) & (chi1 < 120))),
            "chi1_frac_t": float(np.mean(np.abs(chi1) >= 120)),
            "ring_flipped_frac": float(np.mean(np.abs(dchi2) > 90.0)),
            **{f"ring_dist_{args.mut}-{p}_A": float(np.mean(v)) for p, v in cent.items()},
            **({"sasa_sidechain_nm2": float(np.mean(sasa))} if sasa is not None and len(sasa) else {}),
            **({f"helicity_{args.helix}": float(np.mean(helicity[
                (ca.resids >= int(args.helix.split(':')[0])) & (ca.resids <= int(args.helix.split(':')[1]))]))}
               if helicity is not None else {}),
        },
        "contacts": {k: v / n_frames for k, v in contact_counts.items()},
        "hbonds": {k: v / n_frames for k, v in hb_counts.items()},
        "n_frames": n_frames,
    }


def mean_sem(values):
    v = np.asarray(values, dtype=float)
    if len(v) < 2:
        return float(v.mean()), float("nan")
    return float(v.mean()), float(v.std(ddof=1) / np.sqrt(len(v)))


def welch_p(a, b):
    try:
        from scipy.stats import ttest_ind
    except ImportError:
        return float("nan")
    if len(a) < 2 or len(b) < 2:
        return float("nan")
    return float(ttest_ind(a, b, equal_var=False).pvalue)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", required=True)
    ap.add_argument("--systems", nargs="+", required=True)
    ap.add_argument("--nrep", type=int, required=True)
    ap.add_argument("--first", type=int, required=True)
    ap.add_argument("--last", type=int, required=True)
    ap.add_argument("--mut", type=int, required=True)
    ap.add_argument("--skip-ns", type=float, default=50.0)
    ap.add_argument("--stride", type=int, default=5, help="use every n-th frame (10 ps each)")
    ap.add_argument("--core", default="37:140", help="residues for fitting / core RMSD")
    ap.add_argument("--helix", default="64:84", help="helix containing the site (helicity metric)")
    ap.add_argument("--partners", type=int, nargs="*", default=[112, 124, 127])
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    data = {}
    for system in args.systems:
        reps = []
        for r in range(1, args.nrep + 1):
            rep_dir = Path(args.runs) / system / f"rep{r}"
            if not (rep_dir / "md_center.xtc").exists():
                print(f"skip {system}/rep{r}: no md_center.xtc")
                continue
            print(f"analyzing {system}/rep{r}")
            reps.append(analyze_replica(rep_dir, args))
        if reps:
            data[system] = reps
    if not data:
        raise SystemExit("no replicas found")
    ref_sys = args.systems[0] if args.systems[0] in data else next(iter(data))

    # scalar summary
    metrics = sorted({k for reps in data.values() for rep in reps for k in rep["scalars"]})
    with open(out / "summary.tsv", "w") as fh:
        fh.write(f"metric\tsystem\tn_rep\tmean\tsem\tp_vs_{ref_sys}\tper_replica\n")
        for m in metrics:
            ref_vals = [rep["scalars"][m] for rep in data[ref_sys] if m in rep["scalars"]]
            for system, reps in data.items():
                vals = [rep["scalars"][m] for rep in reps if m in rep["scalars"]]
                if not vals:
                    continue
                mu, se = mean_sem(vals)
                p = welch_p(vals, ref_vals) if system != ref_sys else float("nan")
                fh.write(f"{m}\t{system}\t{len(vals)}\t{mu:.4f}\t{se:.4f}\t{p:.3g}\t"
                         + ",".join(f"{v:.4f}" for v in vals) + "\n")

    # per-residue RMSF and helicity
    resids = data[ref_sys][0]["rmsf_resids"]
    with open(out / "per_residue.tsv", "w") as fh:
        fh.write("system\tresid\trmsf_A_mean\trmsf_A_sem\thelicity_mean\thelicity_sem\n")
        for system, reps in data.items():
            rmsf = np.array([rep["rmsf"] for rep in reps])
            hel = [rep["helicity"] for rep in reps if rep["helicity"] is not None]
            hel = np.array(hel) if hel else None
            for i, rid in enumerate(resids):
                rm, rs = mean_sem(rmsf[:, i])
                hm, hs = mean_sem(hel[:, i]) if hel is not None else (float("nan"), float("nan"))
                fh.write(f"{system}\t{rid}\t{rm:.3f}\t{rs:.3f}\t{hm:.3f}\t{hs:.3f}\n")

    # contacts and hydrogen bonds of the site (occupancy, mean over replicas)
    for kind, fname in (("contacts", "contacts.tsv"), ("hbonds", "hbonds.tsv")):
        keys = sorted({k for reps in data.values() for rep in reps for k in rep[kind]})
        with open(out / fname, "w") as fh:
            fh.write("partner\t" + "\t".join(f"{s}_mean\t{s}_sem" for s in data) + "\n")
            rows = []
            for k in keys:
                cols, best = [], 0.0
                for system, reps in data.items():
                    mu, se = mean_sem([rep[kind].get(k, 0.0) for rep in reps])
                    cols.append(f"{mu:.3f}\t{se:.3f}")
                    best = max(best, mu)
                if best >= 0.05:
                    rows.append((best, k, cols))
            for _, k, cols in sorted(rows, key=lambda r: -r[0]):
                fh.write(f"{k}\t" + "\t".join(cols) + "\n")

    make_plots(data, args, out, resids)
    print(f"\nresults in {out}")
    print(open(out / "summary.tsv").read())


def make_plots(data, args, out, resids):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available, no plots")
        return
    colors = dict(zip(data, ["#4a4a4a", "#d1495b", "#00798c", "#edae49"]))

    fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    for system, reps in data.items():
        rmsf = np.array([rep["rmsf"] for rep in reps])
        m, s = rmsf.mean(0), rmsf.std(0, ddof=1) / np.sqrt(len(reps)) if len(reps) > 1 else 0
        axes[0].plot(resids, m, color=colors[system], label=system)
        axes[0].fill_between(resids, m - s, m + s, color=colors[system], alpha=0.25, lw=0)
        hel = [rep["helicity"] for rep in reps if rep["helicity"] is not None]
        if hel:
            h = np.array(hel).mean(0)
            axes[1].plot(resids, h, color=colors[system], label=system)
    for ax in axes:
        ax.axvline(args.mut, color="k", ls=":", lw=1)
    axes[0].set_ylabel("C$\\alpha$ RMSF (Å)")
    axes[1].set_ylabel("helicity (DSSP H)")
    axes[1].set_xlabel("residue")
    axes[0].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out / "rmsf_helicity.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 4))
    for system, reps in data.items():
        for i, rep in enumerate(reps):
            ax.plot(rep["time"], rep["rmsd"], color=colors[system], lw=0.6,
                    label=system if i == 0 else None)
    ax.set_xlabel("time (ns)")
    ax.set_ylabel(f"core C$\\alpha$ RMSD (Å), res {args.core}")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out / "rmsd_time.png", dpi=150)
    plt.close(fig)

    fig, axes = plt.subplots(1, len(data), figsize=(4 * len(data), 4), sharex=True, sharey=True)
    axes = np.atleast_1d(axes)
    for ax, (system, reps) in zip(axes, data.items()):
        c1 = np.concatenate([rep["chi1"] for rep in reps])
        c2 = np.concatenate([rep["chi2"] for rep in reps])
        ax.hist2d(c1, c2, bins=72, range=[[-180, 180], [-180, 180]], cmap="Greys")
        ax.set_title(f"{system} ({reps[0]['resname']}{args.mut})")
        ax.set_xlabel("χ1 (°)")
    axes[0].set_ylabel("χ2 (°)")
    fig.tight_layout()
    fig.savefig(out / "chi1_chi2.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
