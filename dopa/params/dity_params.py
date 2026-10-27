#!/usr/bin/env python3
"""Parameters for the 3,3'-dityrosine cross-link (residue DTY, bridge carbon type CP),
compatible with Amber ff99SB(-ILDN) and the DOPA/TOPA parameters of this study.

Model compounds
  dimer   2,2'-dihydroxy-5,5'-dimethylbiphenyl = dityrosine core (methyl for CB),
          halves A and B, bridge CE2(A)-CE2(B)
  cresol  4-methylphenol = Tyr side chain (reference for Tyr-referenced charges)

QM (PySCF, HF/6-31G* with Cartesian d, geometries optimized with geomeTRIC)
  scan_030 ... scan_180   bridge torsion CD2-CE2-CE2'-CD2' frozen in 15 deg steps,
                          both OH hydrogens frozen pointing away from the bridge
                          (no gas-phase H-bond across the bridge), MP2/6-31G* single
                          points on the HF geometries
  resp_*                  dimer at 60/120 deg with two OH patterns, and p-cresol:
                          ESP on Merz-Kollman points for RESP

Fit
  charges   RESP (same protocol as resp_dopa.py) for dimer and cresol; DTY charges
            = ff99SB Tyr + (dimer half - cresol) for the ring/OH atoms, so the
            WT-to-DTY difference only reflects the cross-link (Tyr-referenced)
  bonded    CP-CP bond length and CA-CP-CP angle from the QM geometries; other CP
            terms copied from the aromatic CA terms of parm99; X-CP-CP-X torsion
            (n = 2, 4) fitted to E(MP2) - E(MM without this torsion)
Output: ff/dty.rtp, ff/dty.hdb, ff/dity_bonded.itp, ff/dity_nonbonded.itp,
        ff/dity_atomtypes.atp, params/results/dity_report.txt (+ json/tsv)

  dity_params.py list
  dity_params.py qm  --task NAME --outdir DIR
  dity_params.py fit --workdir DIR --ffbase GMX_TOP/amber99sb-ildn.ff --outdir RESULTS --ffdir FF
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "src"))
import resp_dopa as R                                    # noqa: E402
from amberff import AmberFF, Molecule, angle            # noqa: E402

HARTREE = 2625.499639                                     # kJ/mol
B_CC = 1.49                                               # A, start value of the bridge

TYPES = {"CB": "CT", "HB1": "HC", "HB2": "HC", "HB3": "HC", "CG": "CA", "CD1": "CA", "HD1": "HA",
         "CD2": "CA", "HD2": "HA", "CE1": "CA", "HE1": "HA", "CE2": "CA", "HE2": "HA", "CZ": "CA",
         "OH": "OH", "HH": "HO"}
HALF_BONDS = [("CB", "HB1"), ("CB", "HB2"), ("CB", "HB3"), ("CB", "CG"), ("CG", "CD1"), ("CG", "CD2"),
              ("CD1", "HD1"), ("CD1", "CE1"), ("CE1", "HE1"), ("CE1", "CZ"), ("CZ", "OH"), ("OH", "HH"),
              ("CZ", "CE2"), ("CE2", "CD2"), ("CD2", "HD2")]
HALF_IMPROPERS = [("CD1", "CD2", "CG", "CB"), ("CG", "CE1", "CD1", "HD1"), ("CD1", "CZ", "CE1", "HE1"),
                  ("CE1", "CE2", "CZ", "OH"), ("CG", "CE2", "CD2", "HD2")]
RING_ANGLE = {"CZ": 0.0, "CE2": -60.0, "CD2": -120.0, "CG": 180.0, "CD1": 120.0, "CE1": 60.0}
RING_OH = ["CG", "CD1", "HD1", "CD2", "HD2", "CE1", "HE1", "CE2", "CZ", "OH", "HH"]

# 0-15 deg (cis) would put the two ortho oxygens on top of each other; the even
# Fourier terms of the bridge torsion are fixed by 30-180 deg
TASKS = {f"scan_{p:03d}": ("dimer", p, "away", "away", True) for p in range(30, 181, 15)}
FIT_EMAX = 50.0                                           # kJ/mol: scan points above are not fitted
TASKS.update({
    "resp_060_aa": ("dimer", 60, "away", "away", False),
    "resp_060_ta": ("dimer", 60, "toward", "away", False),
    "resp_120_aa": ("dimer", 120, "away", "away", False),
    "resp_120_ta": ("dimer", 120, "toward", "away", False),
    "resp_cresol": ("cresol", None, "away", None, False),
})

# DTY residue (Tyr with the cross-linked CE2, no HE2); CZ typed CA as in DAH/TOP
DTY_TYPES = [("N", "N"), ("H", "H"), ("CA", "CT"), ("HA", "H1"), ("CB", "CT"), ("HB1", "HC"), ("HB2", "HC"),
             ("CG", "CA"), ("CD1", "CA"), ("HD1", "HA"), ("CE1", "CA"), ("HE1", "HA"), ("CZ", "CA"),
             ("OH", "OH"), ("HH", "HO"), ("CE2", "CP"), ("CD2", "CA"), ("HD2", "HA"), ("C", "C"), ("O", "O")]
DTY_BONDS = [("N", "H"), ("N", "CA"), ("CA", "HA"), ("CA", "CB"), ("CA", "C"), ("CB", "HB1"), ("CB", "HB2"),
             ("CB", "CG"), ("CG", "CD1"), ("CG", "CD2"), ("CD1", "HD1"), ("CD1", "CE1"), ("CE1", "HE1"),
             ("CE1", "CZ"), ("CZ", "OH"), ("CZ", "CE2"), ("OH", "HH"), ("CE2", "CD2"), ("CD2", "HD2"),
             ("C", "O"), ("-C", "N")]
DTY_IMPROPERS = [("-C", "CA", "N", "H"), ("CA", "+N", "C", "O"), ("CG", "CE2", "CD2", "HD2"),
                 ("CD1", "CZ", "CE1", "HE1"), ("CG", "CE1", "CD1", "HD1"), ("CD1", "CD2", "CG", "CB"),
                 ("CE1", "CE2", "CZ", "OH")]
DTY_HDB = ["1\t1\tH\tN\t-C\tCA", "1\t5\tHA\tCA\tN\tCB\tC", "2\t6\tHB\tCB\tCA\tCG",
           "1\t1\tHD1\tCD1\tCG\tCE1", "1\t1\tHE1\tCE1\tCD1\tCZ", "1\t2\tHH\tOH\tCZ\tCE1",
           "1\t1\tHD2\tCD2\tCG\tCE2"]


# ---------------------------------------------------------------------------
# model geometries (A)
# ---------------------------------------------------------------------------
def unit(a):
    return np.array([np.cos(np.radians(a)), np.sin(np.radians(a)), 0.0])


def build_half(oh):
    """p-cresol in the xy plane, ring centre at the origin; HH towards CE1 ('away') or CE2 ('toward')."""
    X = {n: 1.39 * unit(a) for n, a in RING_ANGLE.items()}
    X["CB"] = X["CG"] + 1.51 * unit(180)
    for h, c, a in (("HD1", "CD1", 120), ("HE1", "CE1", 60), ("HD2", "CD2", -120), ("HE2", "CE2", -60)):
        X[h] = X[c] + 1.08 * unit(a)
    X["OH"] = X["CZ"] + 1.36 * unit(0)
    X["HH"] = R.place(X["CE1"] if oh == "away" else X["CE2"], X["CZ"], X["OH"], 0.96, 113.0, 0.0)
    for k, tor in enumerate((90.0, 210.0, 330.0), start=1):
        X[f"HB{k}"] = R.place(X["CD1"], X["CG"], X["CB"], 1.09, 109.5, tor)
    return X


def rotate(points, p0, axis, ang):
    k = axis / np.linalg.norm(axis)
    v = points - p0
    c, s = np.cos(np.radians(ang)), np.sin(np.radians(ang))
    return p0 + v * c + np.cross(k, v) * s + np.outer(v @ k, k) * (1 - c)


def build_dimer(phi, oh_a, oh_b):
    A, B = build_half(oh_a), build_half(oh_b)
    del A["HE2"], B["HE2"]
    names = list(B)
    pb = np.array([B[n] for n in names]) * np.array([-1.0, -1.0, 1.0])     # 180 deg about z
    u = A["CE2"] / np.linalg.norm(A["CE2"])
    pb += A["CE2"] + B_CC * u - pb[names.index("CE2")]
    i_ce, i_cd = names.index("CE2"), names.index("CD2")
    cur = R.dihedral(A["CD2"], A["CE2"], pb[i_ce], pb[i_cd])
    for sign in (1.0, -1.0):
        trial = rotate(pb, A["CE2"], pb[i_ce] - A["CE2"], sign * (phi - cur))
        if abs((R.dihedral(A["CD2"], A["CE2"], trial[i_ce], trial[i_cd]) - phi + 180) % 360 - 180) < 1e-6:
            break
    else:
        raise RuntimeError("could not set the bridge torsion")
    atoms = [("A", n, n[0] if not n.startswith("O") else "O", A[n]) for n in A]
    atoms += [("B", n, n[0] if not n.startswith("O") else "O", trial[i]) for i, n in enumerate(names)]
    return atoms


def build_task(name):
    kind, phi, oh_a, oh_b, _ = TASKS[name]
    if kind == "cresol":
        X = build_half(oh_a)
        return [("C", n, "O" if n.startswith("O") else n[0], X[n]) for n in X]
    return build_dimer(phi, oh_a, oh_b)


def idx(atoms, res, name):
    return next(i for i, a in enumerate(atoms) if a[0] == res and a[1] == name)


def constraints(name, atoms):
    """1-based dihedrals to freeze: bridge torsion (dimer) and the OH hydrogen orientations."""
    kind, _, oh_a, oh_b, _ = TASKS[name]
    out = []
    if kind == "dimer":
        out.append([idx(atoms, r, n) for r, n in (("A", "CD2"), ("A", "CE2"), ("B", "CE2"), ("B", "CD2"))])
        halves = (("A", oh_a), ("B", oh_b))
    else:
        halves = (("C", oh_a),)
    for res, oh in halves:
        ref = "CE1" if oh == "away" else "CE2"
        out.append([idx(atoms, res, n) for n in (ref, "CZ", "OH", "HH")])
    return [[i + 1 for i in d] for d in out]


# ---------------------------------------------------------------------------
# QM
# ---------------------------------------------------------------------------
def cmd_qm(args):
    from pyscf import gto, scf, mp
    from pyscf.geomopt.geometric_solver import optimize

    name, out = args.task, Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    atoms = build_task(name)
    elems = [a[2] for a in atoms]
    xyz = np.array([a[3] for a in atoms])
    R.write_xyz(out / f"{name}_start.xyz", elems, xyz, f"{name} start")
    cons = out / f"{name}_constraints.txt"
    cons.write_text("$freeze\n" + "".join(f"dihedral {' '.join(map(str, d))}\n" for d in constraints(name, atoms)))

    def make_mol(coords, log):
        return gto.M(atom=[(e, tuple(x)) for e, x in zip(elems, coords)], unit="Angstrom", basis="6-31g*",
                     cart=True, charge=0, spin=0, verbose=4, output=str(out / log), max_memory=args.memory_mb)

    mf = scf.RHF(make_mol(xyz, f"{name}_pyscf.log"))
    mf.conv_tol = 1e-9
    mol_eq = optimize(mf, constraints=str(cons), maxsteps=args.maxsteps, prefix=str(out / f"{name}_geomeTRIC"))
    mf_eq = scf.RHF(mol_eq)
    mf_eq.conv_tol = 1e-10
    e_hf = mf_eq.kernel()
    if not mf_eq.converged:
        raise SystemExit(f"{name}: final SCF not converged")
    meta = {"name": name, "task": TASKS[name], "atoms": [[a[0], a[1], a[2]] for a in atoms], "e_hf": e_hf}
    x_eq = mol_eq.atom_coords(unit="Angstrom")
    if TASKS[name][4]:                                      # scan point: MP2 single point
        meta["e_mp2"] = mp.MP2(mf_eq).run().e_tot
    else:                                                   # RESP conformer: ESP on MK points
        pts = R.mk_points(elems, x_eq) * R.ANG2BOHR
        V = R.esp_at_points(mol_eq, mf_eq.make_rdm1(), pts)
        with open(out / f"{name}_vpot.out", "w") as fh:
            fh.write(f"{len(pts)}\n")
            for p, v in zip(pts, V):
                fh.write(f"{p[0]:16.8f} {p[1]:16.8f} {p[2]:16.8f} {v:18.10e}\n")
    if TASKS[name][0] == "dimer":
        d = [idx(atoms, r, n) for r, n in (("A", "CD2"), ("A", "CE2"), ("B", "CE2"), ("B", "CD2"))]
        meta["phi"] = R.dihedral(*x_eq[d])
    R.write_xyz(out / f"{name}.xyz", elems, x_eq, f"{name} HF/6-31G* optimized E={e_hf:.8f}")
    (out / f"{name}.json").write_text(json.dumps(meta, indent=1))
    print(f"{name}: E(HF) = {e_hf:.8f}" + (f", E(MP2) = {meta['e_mp2']:.8f}" if "e_mp2" in meta else "")
          + (f", phi = {meta['phi']:.1f}" if "phi" in meta else ""))


# ---------------------------------------------------------------------------
# fit
# ---------------------------------------------------------------------------
def resp_two_stage(atoms, confs, equiv1, methyl_groups):
    """Two-stage RESP, molecule neutral. equiv1: classes for stage 1 (symmetry);
    methyl_groups: [(carbon indices, hydrogen indices)] refitted in stage 2."""
    n = len(atoms)
    groups = [(list(range(n)), 0.0)]
    heavy = [i for i, a in enumerate(atoms) if a[2] != "H"]
    q1 = R.resp_fit(confs, n, {}, groups, equiv1, heavy, a=0.0005)
    refit = {i for c, h in methyl_groups for i in list(c) + list(h)}
    fixed = {i: q1[i] for i in range(n) if i not in refit}
    equiv2 = [list(c) for c, h in methyl_groups] + [list(h) for c, h in methyl_groups]
    restr = [i for c, h in methyl_groups for i in c]
    q2 = R.resp_fit(confs, n, fixed, groups, equiv2, restr, a=0.001)
    return q2, R.rrms(confs, q2)


def load(workdir, names):
    meta, confs = {}, []
    for nm in names:
        d = Path(workdir) / nm
        m = json.loads((d / f"{nm}.json").read_text())
        elems, xyz = R.read_xyz(d / f"{nm}.xyz")
        meta[nm] = (m, xyz)
        if (d / f"{nm}_vpot.out").exists():
            P, V = R.read_vpot(d / f"{nm}_vpot.out")
            confs.append((xyz * R.ANG2BOHR, P, V))
    return meta, confs


def dimer_molecule(atoms, charges):
    types = [("CP" if a[1] == "CE2" else TYPES[a[1]]) for a in atoms]
    bonds = [(idx(atoms, r, i), idx(atoms, r, j)) for r in ("A", "B") for i, j in HALF_BONDS]
    bonds.append((idx(atoms, "A", "CE2"), idx(atoms, "B", "CE2")))
    imps = [tuple(idx(atoms, r, x) for x in imp) for r in ("A", "B") for imp in HALF_IMPROPERS]
    for r, o in (("A", "B"), ("B", "A")):
        imps.append((idx(atoms, r, "CZ"), idx(atoms, r, "CD2"), idx(atoms, r, "CE2"), idx(atoms, o, "CE2")))
    return Molecule(types, charges, bonds, imps)


def bonded_itp(b0, th0, torsion_terms):
    lines = ["; 3,3'-dityrosine bridge (type CP), see dopa/params/dity_params.py",
             "[ bondtypes ]",
             "  CA CP         1    0.14000   392459.2 ; as CA-CA (parm99)",
             f"  CP CP         1    {b0:.5f}   292880.0 ; b0 from HF/6-31G*, k for a conjugated C-C single bond",
             "[ angletypes ]",
             "CA  CP  CA           1   120.000    527.184 ; as CA-CA-CA",
             "CA  CA  CP           1   120.000    527.184 ; as CA-CA-CA",
             f"CA  CP  CP           1   {th0:7.3f}    527.184 ; theta0 from HF/6-31G*",
             "CP  CA  HA           1   120.000    418.400 ; as CA-CA-HA",
             "CP  CA  OH           1   120.000    585.760 ; as CA-CA-OH",
             "[ dihedraltypes ]",
             "CA  CP  CA  OH       4      180.00     4.60240     2 ; as CA-CA-CA-OH",
             "CA  CA  CP  CP       4      180.00     4.60240     2 ; bridge substituent, as CA-CA-CA-CT",
             " X   CA  CP  X     9     180.0     15.16700     2 ; as X-CA-CA-X"]
    for ph, k, n in torsion_terms:
        lines.append(f" X   CP  CP  X     9     {ph:5.1f}    {k:9.5f}     {n:d} ; fitted to MP2/6-31G*//HF/6-31G* scan")
    return "\n".join(lines) + "\n"


def cmd_fit(args):
    work = Path(args.workdir)
    missing = [t for t in TASKS if not (work / t / f"{t}.json").exists()]
    if missing:
        raise SystemExit(f"ERROR: QM results missing: {missing}")
    out, ffdir = Path(args.outdir), Path(args.ffdir)
    out.mkdir(parents=True, exist_ok=True)
    report = ["3,3'-dityrosine (DTY) parameters: HF/6-31G* (Cartesian d) geometries, MP2/6-31G* scan "
              f"energies, MK {R.DENSITY} points/A^2, two-stage RESP", ""]

    # --- charges
    resp_d = [t for t in TASKS if t.startswith("resp_") and TASKS[t][0] == "dimer"]
    meta_d, confs_d = load(work, resp_d)
    atoms_d = [tuple(a) for a in meta_d[resp_d[0]][0]["atoms"]]
    equiv_d = [[idx(atoms_d, "A", n), idx(atoms_d, "B", n)] for n in TYPES if n != "HE2"]
    methyl_d = [([idx(atoms_d, "A", "CB"), idx(atoms_d, "B", "CB")],
                 [idx(atoms_d, r, f"HB{k}") for r in ("A", "B") for k in (1, 2, 3)])]
    q_d, err_d = resp_two_stage(atoms_d, confs_d, equiv_d, methyl_d)

    meta_c, confs_c = load(work, ["resp_cresol"])
    atoms_c = [tuple(a) for a in meta_c["resp_cresol"][0]["atoms"]]
    equiv_c = [[idx(atoms_c, "C", a), idx(atoms_c, "C", b)]
               for a, b in (("CD1", "CD2"), ("CE1", "CE2"), ("HD1", "HD2"), ("HE1", "HE2"))]
    methyl_c = [([idx(atoms_c, "C", "CB")], [idx(atoms_c, "C", f"HB{k}") for k in (1, 2, 3)])]
    q_c, err_c = resp_two_stage(atoms_c, confs_c, equiv_c, methyl_c)
    report += [f"RRMS dimer {nm}: {e:.4f}" for nm, e in zip(resp_d, err_d)] + [f"RRMS cresol: {err_c[0]:.4f}", ""]

    half = {n: q_d[idx(atoms_d, "A", n)] for n in TYPES if n != "HE2"}
    cres = {n: q_c[idx(atoms_c, "C", n)] for n in TYPES}
    dty = {k: v for k, v in R.TYR_FF99SB.items() if k != "HE2"}
    for n in RING_OH:
        dty[n] += half[n] - cres[n]
    leftover = sum(dty.values())
    dty["CE2"] -= leftover
    names = [n for n, _ in DTY_TYPES]
    dty = R.round_neutral(dty, names, sink="CE2")
    report.append(f"== DTY charges (Tyr-referenced; leftover {leftover:+.4f} e put on CE2)")
    report.append(f"   {'atom':5s} {'type':4s} {'DTY':>8s} {'TYR ff':>8s} {'dimer':>8s} {'cresol':>8s}")
    for n, t in DTY_TYPES:
        ref = R.TYR_FF99SB.get(n)
        dd = f"{half[n]:8.4f} {cres[n]:8.4f}" if n in half and n in cres else ""
        report.append(f"   {n:5s} {t:4s} {dty[n]:8.4f} {ref:8.4f} {dd}")

    # --- bridge geometry from all optimized dimers
    dims = [t for t in TASKS if TASKS[t][0] == "dimer"]
    meta_all, _ = load(work, dims)
    b, th = [], []
    for nm in dims:
        x = meta_all[nm][1]
        a, c = idx(atoms_d, "A", "CE2"), idx(atoms_d, "B", "CE2")
        b.append(np.linalg.norm(x[a] - x[c]))
        for r, o in (("A", "B"), ("B", "A")):
            for nb in ("CD2", "CZ"):
                th.append(angle(x[idx(atoms_d, r, nb)], x[idx(atoms_d, r, "CE2")], x[idx(atoms_d, o, "CE2")]))
    b0, th0 = float(np.mean(b)) / 10.0, float(np.mean(th))
    report += ["", f"bridge CE2-CE2': {np.mean(b):.4f} +- {np.std(b):.4f} A; CA-CP-CP angle {th0:.1f} +- {np.std(th):.1f} deg"]

    # --- torsion fit: E(MP2) - E(MM without X-CP-CP-X) = const + sum_n c_n sum_d cos(n phi_d)
    tmp = out / "dity_bonded_notorsion.itp"
    tmp.write_text(bonded_itp(b0, th0, []))
    nb = out / "dity_nonbonded_tmp.itp"
    nb.write_text("[ atomtypes ]\nCP 6 12.01 0.0000 A 3.39967e-01 3.59824e-01\n")
    ff = AmberFF(args.ffbase, extra=[nb, tmp])
    scans = sorted(t for t in TASKS if TASKS[t][4])
    mol = dimer_molecule(atoms_d, [q_d[i] for i in range(len(atoms_d))])
    bridge = [d for d in mol.propers if {d[1], d[2]} == {idx(atoms_d, "A", "CE2"), idx(atoms_d, "B", "CE2")}]
    rows, A, y = [], [], []
    for nm in scans:
        m, x = meta_all[nm]
        xn = x / 10.0
        e_mm = mol.energy(ff, xn, skip_proper=lambda t: t[1] == "CP" and t[2] == "CP")
        phis = [R.dihedral(*x[list(d)]) for d in bridge]
        A.append([1.0] + [sum(np.cos(np.radians(n * p)) for p in phis) for n in (2, 4)])
        y.append(m["e_mp2"] * HARTREE - e_mm)
        rows.append((m["phi"], m["e_mp2"] * HARTREE, e_mm))
    A, y = np.array(A), np.array(y)
    e_qm_all = np.array([r[1] for r in rows])
    use = e_qm_all - e_qm_all.min() <= FIT_EMAX
    coef, *_ = np.linalg.lstsq(A[use], y[use], rcond=None)
    terms = []
    for n, c in zip((2, 4), coef[1:]):
        terms.append((0.0, c, n) if c >= 0 else (180.0, -c, n))   # k(1+cos(n phi - ph))
    fit = A @ coef
    e_qm = np.array([r[1] for r in rows])
    e_mm0 = np.array([r[2] for r in rows])
    e_mm = e_mm0 + fit - coef[0]
    rel = lambda v: v - v.min()
    rms = np.sqrt(np.mean(((rel(e_qm) - rel(e_mm)) ** 2)[use]))
    rms0 = np.sqrt(np.mean(((rel(e_qm) - rel(e_mm0)) ** 2)[use]))
    report += ["", "== bridge torsion X-CP-CP-X (kJ/mol): " +
               ", ".join(f"n={n}: k={k:.3f}, phase={ph:.0f}" for ph, k, n in terms),
               f"   RMS(QM - MM) over the fitted points (dE_MP2 <= {FIT_EMAX:g} kJ/mol, {use.sum()} of {len(use)}): "
               f"{rms:.2f} kJ/mol with the fitted torsion, {rms0:.2f} kJ/mol without",
               f"   {'phi':>6s} {'dE_MP2':>8s} {'dE_MM':>8s} {'dE_MM0':>8s}"]
    tsv = ["phi\tdE_MP2\tdE_MM_fit\tdE_MM_no_torsion"]
    for (phi, _, _), q, mm, m0 in zip(rows, rel(e_qm), rel(e_mm), rel(e_mm0)):
        report.append(f"   {phi:6.1f} {q:8.2f} {mm:8.2f} {m0:8.2f}")
        tsv.append(f"{phi:.1f}\t{q:.3f}\t{mm:.3f}\t{m0:.3f}")
    tmp.unlink()
    nb.unlink()

    # --- write force-field pieces
    lines = ["", "[ DTY ]", " ; 3,3'-dityrosine half (cross-linked Tyr, bridge carbon CE2 of type CP), "
             "RESP HF/6-31G* Tyr-referenced, see dopa/README.md", " [ atoms ]"]
    for k, (n, t) in enumerate(DTY_TYPES, start=1):
        lines.append(f"  {n:>5s}  {t:>4s}  {dty[n]:10.5f}  {k:3d}")
    lines.append(" [ bonds ]")
    lines += [f"  {i:>5s}  {j:>5s}" for i, j in DTY_BONDS]
    lines.append(" [ impropers ]")
    lines += ["  " + "  ".join(f"{x:>5s}" for x in imp) for imp in DTY_IMPROPERS]
    (ffdir / "dty.rtp").write_text("\n".join(lines) + "\n")
    (ffdir / "dty.hdb").write_text("\n".join([f"DTY\t{len(DTY_HDB)}"] + DTY_HDB) + "\n")
    (ffdir / "dity_bonded.itp").write_text(bonded_itp(b0, th0, terms))
    (ffdir / "dity_nonbonded.itp").write_text(
        "; biaryl bridge carbon of 3,3'-dityrosine, LJ as CA\n[ atomtypes ]\n"
        "CP           6      12.01    0.0000  A   3.39967e-01  3.59824e-01\n")
    (ffdir / "dity_atomtypes.atp").write_text("CP                12.01000\t; biaryl bridge C (3,3'-dityrosine)\n")
    (out / "dity_torsion.tsv").write_text("\n".join(tsv) + "\n")
    (out / "dty_charges.json").write_text(json.dumps(dty, indent=1))
    text = "\n".join(report)
    (out / "dity_report.txt").write_text(text + "\n")
    print(text)
    print(f"\nwrote {ffdir}/dty.rtp, dty.hdb, dity_bonded.itp, dity_nonbonded.itp, dity_atomtypes.atp")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("list")
    p = sub.add_parser("build")
    p.add_argument("--task", choices=list(TASKS), required=True)
    p.add_argument("--out", required=True)
    p = sub.add_parser("qm")
    p.add_argument("--task", choices=list(TASKS), required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--memory-mb", type=int, default=20000)
    p.add_argument("--maxsteps", type=int, default=300)
    p = sub.add_parser("fit")
    p.add_argument("--workdir", required=True)
    p.add_argument("--ffbase", required=True, help="GROMACS amber99sb-ildn.ff directory")
    p.add_argument("--outdir", required=True)
    p.add_argument("--ffdir", required=True)
    args = ap.parse_args()
    if args.cmd == "list":
        print("\n".join(TASKS))
    elif args.cmd == "build":
        atoms = build_task(args.task)
        R.write_xyz(args.out, [a[2] for a in atoms], np.array([a[3] for a in atoms]), args.task)
    else:
        {"qm": cmd_qm, "fit": cmd_fit}[args.cmd](args)


if __name__ == "__main__":
    main()
