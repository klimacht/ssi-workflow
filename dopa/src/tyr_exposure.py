#!/usr/bin/env python3
"""Time-resolved solvent exposure of the reactive tyrosine atoms and Tyr-Tyr contacts.

Bothen et al. (Free Radic. Res. 2026) relate the site selectivity of nitration and
hydroxylation in Phl p 5a to the atom-level ASA of one relaxed AF3 structure (Table 1):
the phenolic O and the ortho carbons CE1/CE2 (tyrosyl radical, nitration, cross-links)
and the meta carbons CD1/CD2 (OH addition). Free MD shows how much of that exposure
is transient. Per frame and tyrosine:
  ASA of OH, CE1, CE2, CD1, CD2 and of the residue (RSA = ASA / 217.67 A^2, as in the paper)
and per tyrosine pair:
  ring-centroid distance, shortest ortho-ortho C...C distance (3,3'-dityrosine, C-C
  coupling) and shortest ortho C...O distance (isodityrosine, C-O coupling).

ASA: Shrake-Rupley with 960 points per atom, heavy atoms only with united-atom radii
(Chothia 1976, as in NACCESS) and a 1.4 A probe, like SurfaceRacer (used in the paper);
`static` checks the method against Table 1 on a structure file. On the AF3 model this
gives the RSA within 0.05 for 9 of 12 tyrosines and the atom ASA with an RMS deviation of
5.6 A^2 (the paper used the model after a 2 ns restrained relaxation). Explicit hydrogens
make the ring carbons far smaller than in Table 1, so the paper's values are heavy-atom ASA.

  tyr_exposure.py static --pdb MODEL.pdb                          one structure vs. Table 1
  tyr_exposure.py md --runs RUNS --systems wt --nrep 4 --out DIR   replicas (14_analysis.slurm output)
  tyr_exposure.py md ... --reuse                                  re-summarize frames_rep*.tsv

Outputs of `md` in DIR/<system>/:
  frames_rep<i>.tsv   per-frame values
  exposure.tsv        per Tyr and atom: paper, MD start, MD mean +- SEM, fraction of time accessible
  tyr_vs_paper.tsv    one row per Tyr: modification degrees next to static and MD exposure
  pairs.tsv           Tyr-Tyr distances and contact fractions, cross-links / NMR contacts from the paper
  exposure_bfac.pdb   start structure, Tyr atoms with MD mean ASA as B-factor and the
                      fraction of time accessible as occupancy (for PyMOL renderings)
  *.png               time series, distributions, exposure vs. modification degree, pair distances
"""
import argparse
import math
from pathlib import Path

import numpy as np

DATA = Path(__file__).resolve().parents[1] / "data"
PROBE = 1.4
RSA_MAX = 217.67              # A^2, Gly-Tyr-Gly in the paper (SurfaceRacer)
ACCESSIBLE = 1.0              # A^2, "buried" below this in the paper (Singh et al. 2006)
EXPOSED = 5.0                 # A^2, "partially exposed" and above
CONTACT = 4.5                 # A, ortho-ortho / ortho-O distance counted as coupling geometry
REACTIVE = ("OH", "CE1", "CE2", "CD1", "CD2")
PAPER_NAME = {"OH": "O", "CE1": "CE1", "CE2": "CE2", "CD1": "CD1", "CD2": "CD2"}
RING = ("CG", "CD1", "CD2", "CE1", "CE2", "CZ")
ORTHO = ("CE1", "CE2")

# united-atom radii (A) by atom class: NACCESS defaults (Chothia 1976) and ProtOr (Tsai et al. 1999)
RADII = {
    "naccess": {"C3H0": 1.76, "C3H1": 1.76, "C4": 1.87, "N": 1.65, "O1": 1.40, "O2H1": 1.40, "S": 1.85},
    "protor": {"C3H0": 1.61, "C3H1": 1.76, "C4": 1.88, "N": 1.64, "O1": 1.42, "O2H1": 1.46, "S": 1.77},
}
AROMATIC = {"PHE", "TYR", "TRP", "HIS", "HID", "HIE", "HIP", "DAH", "TOP", "DTY"}
SP2_NOH = {                    # trigonal carbons without H, besides the backbone C
    "CG": AROMATIC | {"ASP", "ASH", "ASN"}, "CD": {"GLU", "GLH", "GLN"}, "CZ": {"ARG", "TYR", "DAH", "TOP", "DTY"},
    "CD2": {"TRP"}, "CE2": {"TRP", "DAH", "TOP", "DTY"}, "CD1": {"TOP"},
}
SP2_CH = {                     # aromatic CH
    "CD1": {"PHE", "TYR", "TRP", "DAH", "DTY"},
    "CD2": {"PHE", "TYR", "HIS", "HID", "HIE", "HIP", "DAH", "TOP", "DTY"},
    "CE1": {"PHE", "TYR", "HIS", "HID", "HIE", "HIP", "DAH", "TOP", "DTY"},
    "CE2": {"PHE", "TYR"}, "CZ": {"PHE"}, "CE3": {"TRP"}, "CZ2": {"TRP"}, "CZ3": {"TRP"}, "CH2": {"TRP"},
}


def atom_class(resname, name):
    el = name.lstrip("0123456789")[0]
    if el == "N":
        return "N"
    if el == "S":
        return "S"
    if el == "O":
        hydroxyl = name in ("OG", "OG1", "OH") or (resname in ("DAH", "TOP") and name in ("OE2", "OD1"))
        return "O2H1" if hydroxyl else "O1"
    if el == "C":
        if name == "C" or resname in SP2_NOH.get(name, ()):
            return "C3H0"
        if resname in SP2_CH.get(name, ()):
            return "C3H1"
        return "C4"
    raise ValueError(f"no radius for {resname} {name}")


def radii(resnames, names, scheme="naccess"):
    table = RADII[scheme]
    return np.array([table[atom_class(r, n)] for r, n in zip(resnames, names)])


def sphere_points(n=960):
    """Golden-spiral points on the unit sphere."""
    k = np.arange(n) + 0.5
    z = 1.0 - 2.0 * k / n
    r = np.sqrt(1.0 - z * z)
    phi = np.pi * (3.0 - np.sqrt(5.0)) * k
    return np.c_[r * np.cos(phi), r * np.sin(phi), z]


SPHERE = sphere_points()


def atom_asa(xyz, rad, targets, probe=PROBE, dots=SPHERE):
    """Shrake-Rupley ASA (A^2) of the atoms `targets` (indices into xyz, in A)."""
    R = rad + probe
    out = np.empty(len(targets))
    for k, i in enumerate(targets):
        d2 = ((xyz - xyz[i]) ** 2).sum(1)
        nb = np.flatnonzero(d2 < (R[i] + R) ** 2)
        nb = nb[nb != i]
        pts = xyz[i] + R[i] * dots
        buried = (((pts[:, None, :] - xyz[nb][None]) ** 2).sum(-1) < R[nb] ** 2).any(1)
        out[k] = 4.0 * np.pi * R[i] ** 2 * (1.0 - buried.mean())
    return out


class TyrFrames:
    """Per-frame exposure and pair geometry of the tyrosines in one heavy-atom set."""

    def __init__(self, resnames, resids, names, tyrs=None):
        resnames, resids, names = list(resnames), [int(r) for r in resids], list(names)
        self.tyrs = sorted(tyrs) if tyrs else sorted({r for r, n in zip(resids, resnames) if n == "TYR"})
        idx = {(r, n): i for i, (r, n) in enumerate(zip(resids, names))}
        self.res_atoms = {t: [i for i, r in enumerate(resids) if r == t] for t in self.tyrs}
        for t in self.tyrs:
            missing = [n for n in REACTIVE + RING if (t, n) not in idx]
            if missing:
                raise SystemExit(f"Tyr {t}: atoms {missing} not found (residue {t} is not a Tyr?)")
        self.atom = {(t, n): idx[(t, n)] for t in self.tyrs for n in REACTIVE + RING}
        self.targets = np.array([i for t in self.tyrs for i in self.res_atoms[t]])
        self.pairs = [(a, b) for k, a in enumerate(self.tyrs) for b in self.tyrs[k + 1:]]
        self.columns = [f"Y{t}:{c}" for t in self.tyrs for c in ("RSA",) + tuple(PAPER_NAME[n] for n in REACTIVE)]
        self.columns += [f"{a}-{b}:{c}" for a, b in self.pairs for c in ("ring", "oo", "oO")]

    def values(self, xyz, rad):
        asa = dict(zip(self.targets, atom_asa(xyz, rad, self.targets)))
        row = []
        for t in self.tyrs:
            row.append(sum(asa[i] for i in self.res_atoms[t]) / RSA_MAX)
            row += [asa[self.atom[(t, n)]] for n in REACTIVE]
        x = lambda t, n: xyz[self.atom[(t, n)]]
        for a, b in self.pairs:
            ca = np.mean([x(a, n) for n in RING], 0)
            cb = np.mean([x(b, n) for n in RING], 0)
            oo = min(np.linalg.norm(x(a, p) - x(b, q)) for p in ORTHO for q in ORTHO)
            o_o = min([np.linalg.norm(x(a, p) - x(b, "OH")) for p in ORTHO]
                      + [np.linalg.norm(x(a, "OH") - x(b, q)) for q in ORTHO])
            row += [np.linalg.norm(ca - cb), oo, o_o]
        return np.array(row)


# ---------------------------------------------------------------------------
# reference data from the paper
# ---------------------------------------------------------------------------

def read_tsv(path):
    lines = [l.rstrip("\n") for l in open(path) if l.strip() and not l.startswith("#")]
    head = lines[0].split("\t")
    return [dict(zip(head, l.split("\t"))) for l in lines[1:]]


def num(s):
    return float("nan") if s in ("", "NA", None) else float(s)


def read_ref(path):
    return {int(r["tyr"]): r for r in read_tsv(path)}


def degree(row, kind):
    """(value %, note) of a modification degree at 5:1: total, else single + multi, else single."""
    if row is None:
        return float("nan"), ""
    tot, single, multi = (num(row[f"{kind}_{k}"]) for k in ("total", "single", "multi"))
    if not math.isnan(tot):
        return tot, ""
    if not math.isnan(single) and not math.isnan(multi):
        return single + multi, "single+multi"
    if not math.isnan(single):
        return single, "single only"
    return float("nan"), ""


def paper_value(row, quantity):
    if row is None:
        return float("nan")
    if quantity == "RSA":
        return num(row["RSA"])
    if quantity == "ortho":
        return num(row["ASA_CE1"]) + num(row["ASA_CE2"])
    if quantity == "meta":
        return num(row["ASA_CD1"]) + num(row["ASA_CD2"])
    return num(row[f"ASA_{quantity}"])


# ---------------------------------------------------------------------------
# static: one structure vs. Table 1
# ---------------------------------------------------------------------------

def read_pdb(path):
    resnames, resids, names, xyz = [], [], [], []
    for l in open(path):
        if not l.startswith(("ATOM", "HETATM")) or l[16] not in " A":
            continue
        name, resname = l[12:16].strip(), l[17:20].strip()
        element = l[76:78].strip() or name.lstrip("0123456789")[0]
        if element == "H" or resname in ("HOH", "SOL", "WAT"):
            continue
        resnames.append(resname)
        resids.append(int(l[22:26]))
        names.append(name)
        xyz.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return resnames, resids, names, np.array(xyz)


def cmd_static(args):
    resnames, resids, names, xyz = read_pdb(args.pdb)
    ref = read_ref(args.ref)
    tf = TyrFrames(resnames, resids, names, args.tyr)
    v = dict(zip(tf.columns, tf.values(xyz, radii(resnames, names, args.radii))))
    cols = ["RSA", "O", "CE1", "CE2", "CD1", "CD2"]
    print(f"{args.pdb}: {len(names)} heavy atoms, radii {args.radii}, probe {PROBE} A")
    print("tyr  " + "  ".join(f"{c:>13s}" for c in cols) + "      (paper / here)")
    diffs = []
    for t in tf.tyrs:
        row = ref.get(t)
        cells = []
        for c in cols:
            p, h = paper_value(row, c), v[f"Y{t}:{c}"]
            fmt = "{:.2f}" if c == "RSA" else "{:.1f}"
            cells.append(f"{fmt.format(p) if not math.isnan(p) else 'NA':>6s}/{fmt.format(h):<6s}")
            if c != "RSA" and not math.isnan(p):
                diffs.append(h - p)
        print(f"Y{t:<4d}" + "  ".join(f"{c:>13s}" for c in cells))
    if diffs:
        d = np.array(diffs)
        print(f"atoms: mean diff {d.mean():+.2f}, RMS {np.sqrt((d ** 2).mean()):.2f}, "
              f"max |diff| {np.abs(d).max():.2f} A^2 (n = {len(d)})")


# ---------------------------------------------------------------------------
# md: replicas
# ---------------------------------------------------------------------------

def read_replica(rep_dir, stride):
    """(resnames, resids, names, iterator of (time_ns, xyz)) for the heavy atoms of one replica."""
    import MDAnalysis as mda
    u = mda.Universe(str(rep_dir / "protein_ref.gro"), str(rep_dir / "md_center.xtc"))
    heavy = u.select_atoms("not name H*")
    t_first = u.trajectory[0].time

    def frames():
        for ts in u.trajectory[::stride]:
            yield (ts.time - t_first) / 1000.0, heavy.positions.copy()
    return list(heavy.resnames), list(heavy.resids), list(heavy.names), frames()


def mean_sem(values):
    v = np.asarray(values, dtype=float)
    v = v[~np.isnan(v)]
    if len(v) == 0:
        return float("nan"), float("nan")
    if len(v) == 1:
        return float(v[0]), float("nan")
    return float(v.mean()), float(v.std(ddof=1) / np.sqrt(len(v)))


def quantity_series(rep, t, q):
    """Per-frame values of one quantity; ortho/meta = sum of the two carbons (sum and max)."""
    col = {c: i for i, c in enumerate(rep["columns"])}
    get = lambda n: rep["data"][:, col[f"Y{t}:{n}"]]
    if q == "ortho":
        return get("CE1") + get("CE2"), np.maximum(get("CE1"), get("CE2"))
    if q == "meta":
        return get("CD1") + get("CD2"), np.maximum(get("CD1"), get("CD2"))
    return get(q), get(q)


QUANTITIES = ("RSA", "O", "CE1", "CE2", "CD1", "CD2", "ortho", "meta")


def summarize(reps, tyrs, pairs, ref, pair_ref, skip_ns, out, start_heavy):
    out.mkdir(parents=True, exist_ok=True)
    stats = {}
    for t in tyrs:
        for q in QUANTITIES:
            per = []
            for rep in reps:
                val, gate = quantity_series(rep, t, q)
                m = rep["time"] >= skip_ns
                thr = 0.2 if q == "RSA" else ACCESSIBLE
                per.append((val[0], val[m].mean(), (gate[m] > thr).mean(),
                            (gate[m] > EXPOSED).mean() if q != "RSA" else float("nan")))
            per = np.array(per)
            stats[(t, q)] = {
                "paper": paper_value(ref.get(t), q),
                "start": float(per[:, 0].mean()),
                "mean": mean_sem(per[:, 1]), "acc": mean_sem(per[:, 2]), "exp5": mean_sem(per[:, 3]),
                "per_rep": per[:, 1],
            }

    f2 = lambda x: "NA" if math.isnan(x) else f"{x:.3f}"
    with open(out / "exposure.tsv", "w") as fh:
        fh.write("# ASA in A^2 (RSA dimensionless); ortho = CE1+CE2, meta = CD1+CD2; md_start = first production\n"
                 f"# frame; md_mean over t >= {skip_ns} ns, mean +- SEM over replicas; frac_accessible = fraction of\n"
                 f"# frames with ASA > {ACCESSIBLE} A^2 (ortho/meta: either carbon; RSA: > 0.2), frac_gt5: > {EXPOSED} A^2\n")
        fh.write("tyr\tquantity\tpaper\tmd_start\tmd_mean\tmd_sem\tfrac_accessible\tfrac_accessible_sem"
                 "\tfrac_gt5\tfrac_gt5_sem\tper_replica\n")
        for (t, q), s in stats.items():
            fh.write(f"{t}\t{q}\t{f2(s['paper'])}\t{f2(s['start'])}\t{f2(s['mean'][0])}\t{f2(s['mean'][1])}"
                     f"\t{f2(s['acc'][0])}\t{f2(s['acc'][1])}\t{f2(s['exp5'][0])}\t{f2(s['exp5'][1])}\t"
                     + ",".join(f"{x:.3f}" for x in s["per_rep"]) + "\n")

    with open(out / "tyr_vs_paper.tsv", "w") as fh:
        fh.write(f"# modification degrees (%) at ONOO-:Y = 5:1 from the paper; exposure: paper (static) / MD mean;"
                 f" acc = fraction of time accessible (> {ACCESSIBLE} A^2)\n")
        fh.write("tyr\tsecondary\tND\tND_note\tHD\tHD_note\tRSA_paper\tRSA_md\tO_paper\tO_md\tO_acc"
                 "\tortho_paper\tortho_md\tortho_acc\tmeta_paper\tmeta_md\tmeta_acc\n")
        for t in tyrs:
            row = ref.get(t)
            nd, ndn = degree(row, "ND")
            hd, hdn = degree(row, "HD")
            cells = [str(t), row["secondary"] if row else "NA", f2(nd), ndn or "-", f2(hd), hdn or "-"]
            for q in ("RSA", "O", "ortho", "meta"):
                s = stats[(t, q)]
                cells += [f2(s["paper"]), f2(s["mean"][0])] + ([] if q == "RSA" else [f2(s["acc"][0])])
            fh.write("\t".join(cells) + "\n")

    pstats = {}
    with open(out / "pairs.tsv", "w") as fh:
        fh.write(f"# distances in A; oo = shortest ortho-ortho C...C (3,3'-dityrosine), oO = shortest ortho C...OH\n"
                 f"# (isodityrosine); contact = fraction of frames <= {CONTACT} A; mean +- SEM over replicas,\n"
                 f"# t >= {skip_ns} ns; crosslink / nmr_contact from the paper (NA = not named)\n")
        fh.write("pair\tcrosslink\tpersists\tnmr_contact\tring_start\tring_mean\tring_sem\too_start\too_mean\too_sem"
                 "\too_min\too_contact\too_contact_sem\toO_mean\toO_sem\toO_contact\toO_contact_sem\n")
        for a, b in pairs:
            key = f"{a}-{b}"
            per = []
            for rep in reps:
                col = {c: i for i, c in enumerate(rep["columns"])}
                d = {c: rep["data"][:, col[f"{key}:{c}"]] for c in ("ring", "oo", "oO")}
                m = rep["time"] >= skip_ns
                per.append((d["ring"][0], d["ring"][m].mean(), d["oo"][0], d["oo"][m].mean(), d["oo"][m].min(),
                            (d["oo"][m] <= CONTACT).mean(), d["oO"][m].mean(), (d["oO"][m] <= CONTACT).mean()))
            per = np.array(per)
            pstats[key] = per
            pr = pair_ref.get(key, {})
            ring, oo, ooc, oO, oOc = (mean_sem(per[:, k]) for k in (1, 3, 5, 6, 7))
            fh.write("\t".join([key, pr.get("crosslink", "NA"), pr.get("persists", "NA"), pr.get("nmr_contact", "NA"),
                                f2(per[:, 0].mean()), f2(ring[0]), f2(ring[1]), f2(per[:, 2].mean()), f2(oo[0]),
                                f2(oo[1]), f2(per[:, 4].min()), f2(ooc[0]), f2(ooc[1]), f2(oO[0]), f2(oO[1]),
                                f2(oOc[0]), f2(oOc[1])]) + "\n")

    resnames, resids, names, xyz = start_heavy
    write_bfac_pdb(out / "exposure_bfac.pdb", resnames, resids, names, xyz, stats, tyrs, skip_ns)
    return stats, pstats


def write_bfac_pdb(path, resnames, resids, names, xyz, stats, tyrs, skip_ns):
    """Start structure; Tyr reactive atoms: B = MD mean ASA, occupancy = fraction of time accessible."""
    bq = {}
    for t in tyrs:
        for n in REACTIVE:
            s = stats[(t, PAPER_NAME[n])]
            bq[(t, n)] = (s["mean"][0], s["acc"][0])
    lines = [f"REMARK   1 tyr_exposure.py: B = mean ASA (A^2) over t >= {skip_ns} ns and replicas,",
             f"REMARK   1 occupancy = fraction of frames with ASA > {ACCESSIBLE} A^2; only Tyr OH/CE1/CE2/CD1/CD2"]
    for k, (rn, ri, nm, x) in enumerate(zip(resnames, resids, names, xyz), start=1):
        b, q = bq.get((int(ri), nm), (0.0, 0.0))
        nm4 = nm if len(nm) == 4 else f" {nm:<3s}"
        lines.append(f"ATOM  {k % 100000:5d} {nm4} {rn:>3s} A{int(ri) % 10000:4d}    "
                     f"{x[0]:8.3f}{x[1]:8.3f}{x[2]:8.3f}{q:6.2f}{b:6.2f}          {nm.lstrip('0123456789')[0]:>2s}")
    lines.append("END")
    Path(path).write_text("\n".join(lines) + "\n")


def read_frames(path):
    lines = open(path).read().split("\n")
    head = lines[0].split("\t")
    data = np.array([[float(x) for x in l.split("\t")] for l in lines[1:] if l.strip()])
    return data[:, 0], head[1:], data[:, 1:]


def cmd_md(args, reader=read_replica):
    ref = read_ref(args.ref)
    pair_ref = {r["pair"]: r for r in read_tsv(args.pairs)}
    for system in args.systems:
        out = Path(args.out) / system
        out.mkdir(parents=True, exist_ok=True)
        if args.reuse:                  # new reference values or skip time: no SASA recomputation
            start = read_pdb(out / "exposure_bfac.pdb")
            tf = TyrFrames(*start[:3], args.tyr)
            reps = []
            for r in range(1, args.nrep + 1):
                f = out / f"frames_rep{r}.tsv"
                if f.exists():
                    times, columns, data = read_frames(f)
                    if columns != tf.columns:
                        raise SystemExit(f"{f}: columns do not match the tyrosines {tf.tyrs}")
                    reps.append({"time": times, "data": data, "columns": columns, "rep": r})
            print(f"{system}: reusing {len(reps)} frame files in {out}")
        else:
            reps, tf, start = compute_replicas(args, system, out, reader)
        if not reps:
            print(f"{system}: no replicas")
            continue
        stats, pstats = summarize(reps, tf.tyrs, tf.pairs, ref, pair_ref, args.skip_ns, out, start)
        make_plots(reps, tf, ref, stats, pair_ref, args.skip_ns, out)
        print(f"results in {out}")
        print(open(out / "tyr_vs_paper.tsv").read())
        print(open(out / "pairs.tsv").read())


def compute_replicas(args, system, out, reader):
    reps, tf, start = [], None, None
    for r in range(1, args.nrep + 1):
        rep_dir = Path(args.runs) / system / f"rep{r}"
        if not (rep_dir / "md_center.xtc").exists() and reader is read_replica:
            print(f"skip {system}/rep{r}: no md_center.xtc")
            continue
        resnames, resids, names, frames = reader(rep_dir, args.stride)
        if tf is None:
            tf = TyrFrames(resnames, resids, names, args.tyr)
            rad = radii(resnames, names, args.radii)
            print(f"{system}: tyrosines {tf.tyrs}, {len(names)} heavy atoms")
        times, rows = [], []
        for t_ns, xyz in frames:
            if start is None:
                start = (resnames, resids, names, xyz)
            times.append(t_ns)
            rows.append(tf.values(xyz, rad))
        data = np.array(rows)
        with open(out / f"frames_rep{r}.tsv", "w") as fh:
            fh.write("time_ns\t" + "\t".join(tf.columns) + "\n")
            for t_ns, row in zip(times, data):
                fh.write(f"{t_ns:.3f}\t" + "\t".join(f"{v:.3f}" for v in row) + "\n")
        print(f"   rep{r}: {len(times)} frames, {times[-1]:.1f} ns")
        reps.append({"time": np.array(times), "data": data, "columns": tf.columns, "rep": r})
    return reps, tf, start


# ---------------------------------------------------------------------------
# plots
# ---------------------------------------------------------------------------

def running_mean(y, n):
    if n <= 1 or len(y) < n:
        return y
    c = np.cumsum(np.insert(y, 0, 0.0))
    return np.r_[np.full(n - 1, np.nan), (c[n:] - c[:-n]) / n]


def make_plots(reps, tf, ref, stats, pair_ref, skip_ns, out):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available, no plots")
        return
    colors = ["#4a4a4a", "#d1495b", "#00798c", "#edae49", "#66a182", "#2e4057"]
    tyrs = tf.tyrs
    dt = np.median(np.diff(reps[0]["time"])) if len(reps[0]["time"]) > 1 else 1.0
    win = max(1, int(round(1.0 / dt)))                     # 1 ns running mean

    # time series
    qs = [("O", "phenol O"), ("ortho", "ortho CE1+CE2"), ("meta", "meta CD1+CD2")]
    fig, axes = plt.subplots(len(tyrs), 3, figsize=(12, 1.9 * len(tyrs) + 0.6), sharex=True, squeeze=False)
    for i, t in enumerate(tyrs):
        for j, (q, title) in enumerate(qs):
            ax = axes[i, j]
            for k, rep in enumerate(reps):
                val, _ = quantity_series(rep, t, q)
                ax.plot(rep["time"], running_mean(val, win), color=colors[k % len(colors)], lw=0.7,
                        label=f"rep{rep['rep']}" if i == 0 and j == 0 else None)
            p = stats[(t, q)]["paper"]
            if not math.isnan(p):
                ax.axhline(p, color="k", ls="--", lw=1, label="paper (static)" if i == 0 and j == 0 else None)
            ax.axvline(skip_ns, color="k", ls=":", lw=0.8)
            if i == 0:
                ax.set_title(title)
            if j == 0:
                ax.set_ylabel(f"Y{t}\nASA (Å²)")
    for ax in axes[-1]:
        ax.set_xlabel("time (ns)")
    axes[0, 0].legend(frameon=False, fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(out / "tyr_asa_time.png", dpi=150)
    plt.close(fig)

    # distributions per atom vs. paper
    atoms = ["O", "CE1", "CE2", "CD1", "CD2"]
    fig, axes = plt.subplots(1, len(tyrs), figsize=(2.6 * len(tyrs), 3.4), sharey=True, squeeze=False)
    for ax, t in zip(axes[0], tyrs):
        vals = []
        for q in atoms:
            v = np.concatenate([quantity_series(rep, t, q)[0][rep["time"] >= skip_ns] for rep in reps])
            vals.append(v)
        ax.boxplot(vals, whis=(5, 95), showfliers=False)
        ax.plot(range(1, len(atoms) + 1), [stats[(t, q)]["paper"] for q in atoms], "D", color="#d1495b",
                ms=5, label="paper (static)")
        ax.plot(range(1, len(atoms) + 1), [stats[(t, q)]["start"] for q in atoms], "o", mfc="none",
                color="#00798c", ms=5, label="MD start")
        ax.set_xticks(range(1, len(atoms) + 1))
        ax.set_xticklabels(atoms, fontsize=8)
        ax.axhline(ACCESSIBLE, color="k", ls=":", lw=0.8)
        ax.set_title(f"Y{t}")
    axes[0, 0].set_ylabel("ASA (Å²), MD 5–95 %")
    axes[0, 0].legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(out / "tyr_asa_dist.png", dpi=150)
    plt.close(fig)

    # exposure vs. modification degree
    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    for ax, q, kind, label in ((axes[0], "ortho", "ND", "nitration"), (axes[1], "meta", "HD", "hydroxylation")):
        for t in tyrs:
            y, note = degree(ref.get(t), kind)
            if math.isnan(y):
                continue
            s = stats[(t, q)]
            ax.plot([s["paper"], s["mean"][0]], [y, y], color="0.6", lw=0.8)
            ax.plot(s["paper"], y, "o", mfc="none", color="k")
            ax.errorbar(s["mean"][0], y, xerr=None if math.isnan(s["mean"][1]) else s["mean"][1],
                        fmt="o", color="#d1495b")
            ax.annotate(f"Y{t}" + ("*" if note else ""), (s["mean"][0], y), textcoords="offset points",
                        xytext=(4, 3), fontsize=8)
        ax.set_xlabel(f"{q} ASA (Å²): open = paper, filled = MD mean")
        ax.set_ylabel(f"{kind} at 5:1 (%)")
        ax.set_title(label)
    fig.text(0.01, 0.01, "* single+multi or single only (no total in the main text)", fontsize=7)
    fig.tight_layout()
    fig.savefig(out / "exposure_vs_modification.png", dpi=150)
    plt.close(fig)

    # pair distances
    show = [(a, b) for a, b in tf.pairs
            if min(np.nanmin(rep["data"][:, rep["columns"].index(f"{a}-{b}:oo")]) for rep in reps) < 8.0]
    if show:
        fig, axes = plt.subplots(len(show), 1, figsize=(9, 1.8 * len(show) + 0.6), sharex=True, squeeze=False)
        for ax, (a, b) in zip(axes[:, 0], show):
            for k, rep in enumerate(reps):
                d = rep["data"][:, rep["columns"].index(f"{a}-{b}:oo")]
                ax.plot(rep["time"], running_mean(d, win), color=colors[k % len(colors)], lw=0.7)
            ax.axhline(CONTACT, color="k", ls="--", lw=0.8)
            pr = pair_ref.get(f"{a}-{b}", {})
            tag = ", cross-link in paper" if pr.get("crosslink") == "yes" else ""
            tag += ", NMR contact" if pr.get("nmr_contact") == "yes" else ""
            ax.set_ylabel("Å")
            ax.set_title(f"Y{a}–Y{b} shortest ortho C···C{tag}", fontsize=9)
        axes[-1, 0].set_xlabel("time (ns)")
        fig.tight_layout()
        fig.savefig(out / "tyr_pairs_time.png", dpi=150)
        plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("static", "md"):
        p = sub.add_parser(name)
        p.add_argument("--ref", default=str(DATA / "bothen2026_tyr.tsv"))
        p.add_argument("--tyr", type=int, nargs="*", help="tyrosine residue numbers (default: all TYR)")
        p.add_argument("--radii", choices=sorted(RADII), default="naccess")
        if name == "static":
            p.add_argument("--pdb", required=True)
        else:
            p.add_argument("--pairs", default=str(DATA / "bothen2026_pairs.tsv"))
            p.add_argument("--runs", required=True)
            p.add_argument("--systems", nargs="+", default=["wt"])
            p.add_argument("--nrep", type=int, required=True)
            p.add_argument("--skip-ns", type=float, default=50.0)
            p.add_argument("--stride", type=int, default=5, help="use every n-th frame (10 ps each)")
            p.add_argument("--out", required=True)
            p.add_argument("--reuse", action="store_true",
                           help="summarize existing frames_rep*.tsv (e.g. after editing the reference table)")
    args = ap.parse_args()
    {"static": cmd_static, "md": cmd_md}[args.cmd](args)


if __name__ == "__main__":
    main()
