#!/usr/bin/env python3
"""RESP charges for hydroxylated tyrosines compatible with Amber ff94/ff99SB(-ILDN):
L-DOPA (3,4-dihydroxyphenylalanine, residue DAH) and TOPA (2,4,5-trihydroxy-
phenylalanine, residue TOP).

Recipe (Cornell et al. 1995; Cieplak et al. 1995; Bayly et al. 1993):
  * capped dipeptides Ace-X-NMe, backbone in alpha (phi,psi = -60,-40) and
    beta (-120,140) conformation, phi/psi constrained, everything else relaxed
  * HF/6-31G* (Cartesian d, as Gaussian's default 6D) geometry and electrostatic
    potential with PySCF; constrained optimization with geomeTRIC
  * Merz-Kollman sampling shells at 1.4/1.6/1.8/2.0 x vdW radius, 1 point/A^2
  * two-stage multi-conformation RESP (a = 0.0005 / 0.001, b = 0.1),
    backbone N, H, C, O fixed to the ff99SB values, residue and caps neutral
  * DOPA: both intramolecular catechol H-bond patterns (a: 3-OH donates,
    b: 4-OH donates) in both backbone conformations -> 4 conformers;
    TOPA: additionally both in-plane orientations of the 2-OH H (e: towards
    CE1, g: towards CG) -> 8 conformers
  * Tyr is fitted with the same protocol as a check against ff99SB
    (symmetric ring atoms equivalenced, as in ff94)
  * charges written to the force field (default): Tyr-referenced,
    ff99SB Tyr + (fit - Tyr fit), so that the residual difference of our
    protocol to the original ff94 fit cancels between WT and the modified residue

Subcommands:
  list    print conformer names (optionally of one residue)
  build   write start geometries (xyz) + atom lists (json), for inspection
  qm      build one conformer, optimize it and write <name>.xyz + <name>_vpot.out
          (needs pyscf and geometric, see make_qm_venv.sh)
  fit     fit charges over all conformers, write report + GROMACS rtp/hdb
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

ANG2BOHR = 1.0 / 0.529177210903

# ff99SB backbone charges, identical for all neutral residues
BACKBONE = {"N": -0.4157, "H": 0.2719, "C": 0.5973, "O": -0.5679}

TYR_FF99SB = {
    "N": -0.4157, "H": 0.2719, "CA": -0.0014, "HA": 0.0876, "CB": -0.0152,
    "HB1": 0.0295, "HB2": 0.0295, "CG": -0.0011, "CD1": -0.1906, "HD1": 0.1699,
    "CE1": -0.2341, "HE1": 0.1656, "CZ": 0.3226, "OH": -0.5579, "HH": 0.3992,
    "CE2": -0.2341, "HE2": 0.1656, "CD2": -0.1906, "HD2": 0.1699,
    "C": 0.5973, "O": -0.5679,
}

# Amber atom types. All hydroxylated ring carbons are typed CA so that every
# bonded term exists in parm99 (CA-OH, CA-CA-OH, CA-OH-HO, X-CA-OH-X,
# improper CA-CA-CA-OH). Tyr uses type C for CZ, but C-CA-OH is not defined.
DAH_TYPES = [
    ("N", "N"), ("H", "H"), ("CA", "CT"), ("HA", "H1"), ("CB", "CT"),
    ("HB1", "HC"), ("HB2", "HC"), ("CG", "CA"), ("CD1", "CA"), ("HD1", "HA"),
    ("CE1", "CA"), ("HE1", "HA"), ("CZ", "CA"), ("OH", "OH"), ("HH", "HO"),
    ("CE2", "CA"), ("OE2", "OH"), ("HOE2", "HO"), ("CD2", "CA"), ("HD2", "HA"),
    ("C", "C"), ("O", "O"),
]
TOP_TYPES = [
    ("N", "N"), ("H", "H"), ("CA", "CT"), ("HA", "H1"), ("CB", "CT"),
    ("HB1", "HC"), ("HB2", "HC"), ("CG", "CA"), ("CD1", "CA"), ("OD1", "OH"), ("HOD1", "HO"),
    ("CE1", "CA"), ("HE1", "HA"), ("CZ", "CA"), ("OH", "OH"), ("HH", "HO"),
    ("CE2", "CA"), ("OE2", "OH"), ("HOE2", "HO"), ("CD2", "CA"), ("HD2", "HA"),
    ("C", "C"), ("O", "O"),
]
DAH_BONDS = [("N", "H"), ("N", "CA"), ("CA", "HA"), ("CA", "CB"), ("CA", "C"), ("CB", "HB1"),
             ("CB", "HB2"), ("CB", "CG"), ("CG", "CD1"), ("CG", "CD2"), ("CD1", "HD1"),
             ("CD1", "CE1"), ("CE1", "HE1"), ("CE1", "CZ"), ("CZ", "OH"), ("CZ", "CE2"),
             ("OH", "HH"), ("CE2", "OE2"), ("OE2", "HOE2"), ("CE2", "CD2"), ("CD2", "HD2"),
             ("C", "O"), ("-C", "N")]
DAH_IMPROPERS = [("-C", "CA", "N", "H"), ("CA", "+N", "C", "O"), ("CG", "CE2", "CD2", "HD2"),
                 ("CZ", "CD2", "CE2", "OE2"), ("CD1", "CZ", "CE1", "HE1"),
                 ("CG", "CE1", "CD1", "HD1"), ("CD1", "CD2", "CG", "CB"), ("CE1", "CE2", "CZ", "OH")]
DAH_HDB = ["1\t1\tH\tN\t-C\tCA", "1\t5\tHA\tCA\tN\tCB\tC", "2\t6\tHB\tCB\tCA\tCG",
           "1\t1\tHD1\tCD1\tCG\tCE1", "1\t1\tHE1\tCE1\tCD1\tCZ", "1\t2\tHH\tOH\tCZ\tCE1",
           "1\t2\tHOE2\tOE2\tCE2\tCZ", "1\t1\tHD2\tCD2\tCG\tCE2"]

RESIDUES = {
    "DAH": {
        "title": "L-3,4-dihydroxyphenylalanine (L-DOPA)",
        "types": DAH_TYPES, "bonds": DAH_BONDS, "impropers": DAH_IMPROPERS, "hdb": DAH_HDB,
        "replaced": {"HE2": "CE2"},           # Tyr H replaced by an OH -> its carbon
    },
    "TOP": {
        "title": "2,4,5-trihydroxyphenylalanine (TOPA)",
        "types": TOP_TYPES,
        "bonds": [("CD1", "OD1") if b == ("CD1", "HD1") else b for b in DAH_BONDS] + [("OD1", "HOD1")],
        "impropers": [("CG", "CE1", "CD1", "OD1") if i == ("CG", "CE1", "CD1", "HD1") else i
                      for i in DAH_IMPROPERS],
        "hdb": ["1\t2\tHOD1\tOD1\tCD1\tCE1" if h.startswith("1\t1\tHD1") else h for h in DAH_HDB],
        "replaced": {"HE2": "CE2", "HD1": "CD1"},
    },
}

ALPHA = (-60.0, -40.0)
BETA = (-120.0, 140.0)
CHI1, CHI2 = 180.0, 90.0          # Y65 in the AF3 model: chi1 = 176, chi2 = -106

CONFORMERS = {
    # name: (residue, (phi, psi), H-bond pattern: catechol a/b [_ 2-OH e/g])
    "dah_alpha_a": ("DAH", ALPHA, "a"),
    "dah_alpha_b": ("DAH", ALPHA, "b"),
    "dah_beta_a": ("DAH", BETA, "a"),
    "dah_beta_b": ("DAH", BETA, "b"),
    "tyr_alpha": ("TYR", ALPHA, None),
    "tyr_beta": ("TYR", BETA, None),
}
for _bb, _phipsi in (("alpha", ALPHA), ("beta", BETA)):
    for _cat in "ab":
        for _oh in "eg":
            CONFORMERS[f"top_{_bb}_{_cat}_{_oh}"] = ("TOP", _phipsi, f"{_cat}_{_oh}")

# Symmetry-equivalent atoms (Tyr ring flips); the hydroxylated rings have none
SYMMETRY = {"TYR": [("CD1", "CD2"), ("CE1", "CE2"), ("HD1", "HD2"), ("HE1", "HE2")], "DAH": [], "TOP": []}

MK_RADII = {"H": 1.20, "C": 1.50, "N": 1.50, "O": 1.40}
SHELLS = (1.4, 1.6, 1.8, 2.0)
# Points per A^2: 1 = Merz-Kollman / Gaussian default, as used for ff94. The RESP
# restraint weights (a) are calibrated for it; denser grids weaken the restraints
# and let buried carbons drift (Tyr check: max |dq| 0.20 at 6/A^2 vs 0.06 at 1/A^2).
DENSITY = 1.0


# ----------------------------------------------------------------------------
# geometry
# ----------------------------------------------------------------------------
def place(a, b, c, bond, angle, torsion):
    """Atom d with |cd| = bond, angle(b,c,d) = angle, dihedral(a,b,c,d) = torsion."""
    ang, tor = np.radians(angle), np.radians(torsion)
    bc = (c - b) / np.linalg.norm(c - b)
    n = np.cross(b - a, bc)
    n /= np.linalg.norm(n)
    m = np.cross(n, bc)
    d = np.array([-bond * np.cos(ang),
                  bond * np.sin(ang) * np.cos(tor),
                  bond * np.sin(ang) * np.sin(tor)])
    return c + d[0] * bc + d[1] * m + d[2] * n


def dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w)))


def tetra_one(center, n1, n2, n3, bond):
    """Fourth tetrahedral substituent of a center with three neighbours."""
    u = sum((x - center) / np.linalg.norm(x - center) for x in (n1, n2, n3))
    return center - bond * u / np.linalg.norm(u)


def tetra_two(center, n1, n2, bond):
    """Two remaining tetrahedral substituents of a center with two neighbours."""
    u1 = (n1 - center) / np.linalg.norm(n1 - center)
    u2 = (n2 - center) / np.linalg.norm(n2 - center)
    bis = -(u1 + u2) / np.linalg.norm(u1 + u2)
    perp = np.cross(u1, u2) / np.linalg.norm(np.cross(u1, u2))
    half = np.radians(109.5 / 2)
    return (center + bond * (np.cos(half) * bis + np.sin(half) * perp),
            center + bond * (np.cos(half) * bis - np.sin(half) * perp))


def build_dipeptide(resname, phi, psi, pattern):
    """Ace-X-NMe with X = TYR or DAH. Returns list of (residue, name, element, xyz)."""
    X = {}
    X["ACE", "CH3"] = np.zeros(3)
    X["ACE", "C"] = np.array([1.522, 0.0, 0.0])
    X["ACE", "O"] = place(np.array([0.0, 1.0, 0.0]), X["ACE", "CH3"], X["ACE", "C"], 1.229, 121.0, 0.0)
    X["RES", "N"] = place(X["ACE", "O"], X["ACE", "CH3"], X["ACE", "C"], 1.335, 116.6, 180.0)
    X["RES", "H"] = place(X["ACE", "O"], X["ACE", "C"], X["RES", "N"], 1.010, 119.8, 180.0)
    X["RES", "CA"] = place(X["ACE", "O"], X["ACE", "C"], X["RES", "N"], 1.449, 121.9, 0.0)
    X["RES", "C"] = place(X["ACE", "C"], X["RES", "N"], X["RES", "CA"], 1.522, 111.1, phi)
    # L-configuration: dihedral C-N-CA-CB = -123 (measured in the PhLP5 model)
    X["RES", "CB"] = place(X["RES", "C"], X["RES", "N"], X["RES", "CA"], 1.526, 110.4, -123.0)
    X["RES", "HA"] = tetra_one(X["RES", "CA"], X["RES", "N"], X["RES", "C"], X["RES", "CB"], 1.090)
    X["RES", "O"] = place(X["RES", "N"], X["RES", "CA"], X["RES", "C"], 1.229, 120.4, psi + 180.0)
    X["NME", "N"] = place(X["RES", "N"], X["RES", "CA"], X["RES", "C"], 1.335, 116.6, psi)
    X["NME", "H"] = place(X["RES", "O"], X["RES", "C"], X["NME", "N"], 1.010, 119.8, 180.0)
    X["NME", "CH3"] = place(X["RES", "O"], X["RES", "C"], X["NME", "N"], 1.449, 121.9, 0.0)
    for k, tor in enumerate((60.0, 180.0, 300.0), start=1):
        X["NME", f"HH3{k}"] = place(X["RES", "C"], X["NME", "N"], X["NME", "CH3"], 1.090, 109.5, tor)
        X["ACE", f"HH3{k}"] = place(X["RES", "N"], X["ACE", "C"], X["ACE", "CH3"], 1.090, 109.5, tor - 60.0)

    # side chain
    g = lambda name: X["RES", name]
    X["RES", "CG"] = place(g("N"), g("CA"), g("CB"), 1.510, 114.0, CHI1)
    X["RES", "HB1"], X["RES", "HB2"] = tetra_two(g("CB"), g("CA"), g("CG"), 1.090)
    X["RES", "CD1"] = place(g("CA"), g("CB"), g("CG"), 1.390, 120.8, CHI2)
    X["RES", "CD2"] = place(g("CA"), g("CB"), g("CG"), 1.390, 120.8, CHI2 + 180.0)
    X["RES", "CE1"] = place(g("CB"), g("CG"), g("CD1"), 1.390, 120.0, 180.0)
    X["RES", "CE2"] = place(g("CB"), g("CG"), g("CD2"), 1.390, 120.0, 180.0)
    X["RES", "CZ"] = place(g("CG"), g("CD1"), g("CE1"), 1.390, 120.0, 0.0)
    X["RES", "OH"] = place(g("CD1"), g("CE1"), g("CZ"), 1.360, 120.0, 180.0)
    X["RES", "HD1"] = place(g("CB"), g("CG"), g("CD1"), 1.080, 120.0, 0.0)
    X["RES", "HD2"] = place(g("CB"), g("CG"), g("CD2"), 1.080, 120.0, 0.0)
    X["RES", "HE1"] = place(g("CG"), g("CD1"), g("CE1"), 1.080, 120.0, 180.0)
    if resname == "TYR":
        X["RES", "HE2"] = place(g("CG"), g("CD2"), g("CE2"), 1.080, 120.0, 180.0)
        X["RES", "HH"] = place(g("CE1"), g("CZ"), g("OH"), 0.960, 113.0, 0.0)
    else:
        X["RES", "OE2"] = place(g("CG"), g("CD2"), g("CE2"), 1.360, 120.0, 180.0)
        if pattern[0] == "a":   # 3-OH donates to the 4-O, 4-OH points away
            X["RES", "HOE2"] = place(g("CZ"), g("CE2"), g("OE2"), 0.960, 113.0, 0.0)
            X["RES", "HH"] = place(g("CE1"), g("CZ"), g("OH"), 0.960, 113.0, 0.0)
        else:                   # 4-OH donates to the 3-O, 3-OH points away
            X["RES", "HH"] = place(g("CE2"), g("CZ"), g("OH"), 0.960, 113.0, 0.0)
            X["RES", "HOE2"] = place(g("CD2"), g("CE2"), g("OE2"), 0.960, 113.0, 0.0)
    if resname == "TOP":        # 2-OH (on CD1, para to the 5-OH), H towards CE1 (e) or CG (g)
        X["RES", "OD1"] = place(g("CB"), g("CG"), g("CD1"), 1.360, 120.0, 0.0)
        ref = g("CE1") if pattern.endswith("e") else g("CG")
        X["RES", "HOD1"] = place(ref, g("CD1"), g("OD1"), 0.960, 113.0, 0.0)

    res_order = [n for n, _ in RESIDUES[resname]["types"]] if resname in RESIDUES else list(TYR_FF99SB)
    order = ([("ACE", n) for n in ("CH3", "HH31", "HH32", "HH33", "C", "O")]
             + [("RES", n) for n in res_order]
             + [("NME", n) for n in ("N", "H", "CH3", "HH31", "HH32", "HH33")])
    atoms = []
    for res, name in order:
        elem = "O" if name.startswith("O") else name[0]
        atoms.append((res, name, elem, X[res, name]))
    return atoms


def check_geometry(atoms):
    xyz = np.array([a[3] for a in atoms])
    d = np.linalg.norm(xyz[:, None] - xyz[None], axis=2)
    np.fill_diagonal(d, 9.0)
    if d.min() < 0.9:
        i, j = np.unravel_index(d.argmin(), d.shape)
        raise SystemExit(f"clash in start geometry: {atoms[i][:2]} - {atoms[j][:2]} {d.min():.2f} A")


def index_of(atoms, res, name):
    return next(i for i, a in enumerate(atoms) if a[0] == res and a[1] == name)


def write_xyz(path, elems, xyz, comment):
    with open(path, "w") as fh:
        fh.write(f"{len(elems)}\n{comment}\n")
        for e, x in zip(elems, xyz):
            fh.write(f"{e:2s} {x[0]:12.6f} {x[1]:12.6f} {x[2]:12.6f}\n")


def build_one(name, out):
    """Start geometry, atom list and phi/psi atom indices (0-based) of one conformer."""
    resname, (phi, psi), pattern = CONFORMERS[name]
    atoms = build_dipeptide(resname, phi, psi, pattern)
    check_geometry(atoms)
    idx_phi = [index_of(atoms, *k) for k in (("ACE", "C"), ("RES", "N"), ("RES", "CA"), ("RES", "C"))]
    idx_psi = [index_of(atoms, *k) for k in (("RES", "N"), ("RES", "CA"), ("RES", "C"), ("NME", "N"))]
    xyz = np.array([a[3] for a in atoms])
    got = (dihedral(*xyz[idx_phi]), dihedral(*xyz[idx_psi]))
    assert abs(got[0] - phi) < 0.01 and abs((got[1] - psi + 180) % 360 - 180) < 0.01, got
    elems = [a[2] for a in atoms]
    write_xyz(out / f"{name}_start.xyz", elems, xyz,
              f"{name} Ace-{resname}-NMe phi={phi} psi={psi} pattern={pattern}")
    meta = {"name": name, "resname": resname, "phi": phi, "psi": psi, "pattern": pattern,
            "atoms": [[a[0], a[1], a[2]] for a in atoms]}
    (out / f"{name}.json").write_text(json.dumps(meta, indent=1))
    return elems, xyz, idx_phi, idx_psi


def cmd_build(args):
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    names = [args.only] if args.only else list(CONFORMERS)
    for name in names:
        build_one(name, out)
    print(f"wrote {', '.join(names)} to {out}")


# ----------------------------------------------------------------------------
# ESP sampling points
# ----------------------------------------------------------------------------
def read_xyz(path):
    lines = Path(path).read_text().split("\n")
    n = int(lines[0])
    elems, xyz = [], []
    for line in lines[2:2 + n]:
        f = line.split()
        elems.append(f[0])
        xyz.append([float(v) for v in f[1:4]])
    return elems, np.array(xyz)


def fibonacci_sphere(n):
    i = np.arange(n) + 0.5
    polar = np.arccos(1.0 - 2.0 * i / n)
    azim = np.pi * (1.0 + 5.0 ** 0.5) * i
    return np.column_stack([np.cos(azim) * np.sin(polar), np.sin(azim) * np.sin(polar), np.cos(polar)])


def mk_points(elems, xyz):
    """Merz-Kollman style points (A) on scaled vdW shells, outside all inner spheres."""
    radii = np.array([MK_RADII[e] for e in elems])
    pts = []
    for scale in SHELLS:
        r = scale * radii
        for i in range(len(xyz)):
            n = max(int(DENSITY * 4.0 * np.pi * r[i] ** 2), 10)
            p = xyz[i] + r[i] * fibonacci_sphere(n)
            d = np.linalg.norm(p[:, None, :] - xyz[None, :, :], axis=2)
            pts.append(p[np.all(d >= r[None, :] - 1e-6, axis=1)])
    return np.vstack(pts)


def esp_at_points(mol, dm, pts, chunk=200):
    """Electrostatic potential (a.u.) of nuclei + electron density at points (bohr)."""
    R = mol.atom_coords()
    Z = mol.atom_charges()
    v_nuc = (Z[None, :] / np.linalg.norm(pts[:, None, :] - R[None, :, :], axis=2)).sum(axis=1)
    v_el = np.empty(len(pts))
    for s in range(0, len(pts), chunk):
        ints = mol.intor("int1e_grids", grids=pts[s:s + chunk])
        v_el[s:s + chunk] = np.einsum("gij,ij->g", ints, dm)
    return v_nuc - v_el


def cmd_qm(args):
    from pyscf import gto, scf
    from pyscf.geomopt.geometric_solver import optimize

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    name = args.name
    elems, xyz, idx_phi, idx_psi = build_one(name, out)

    def make_mol(coords, log):
        return gto.M(atom=[(e, tuple(x)) for e, x in zip(elems, coords)], unit="Angstrom",
                     basis="6-31g*", cart=True, charge=0, spin=0, verbose=4,
                     output=str(out / log), max_memory=args.memory_mb)

    opt_xyz = out / f"{name}.xyz"
    if opt_xyz.exists() and "optimized" in opt_xyz.read_text().split("\n")[1]:
        # reuse the optimized geometry (e.g. when only the ESP grid changed)
        elems_eq, xyz_opt = read_xyz(opt_xyz)
        assert elems_eq == elems, f"atom order changed in {opt_xyz}"
        print(f"{name}: reusing optimized geometry from {opt_xyz}")
        mol_eq = make_mol(xyz_opt, f"{name}_pyscf_esp.log")
    else:
        cons = out / f"{name}_constraints.txt"
        cons.write_text("$freeze\n"
                        f"dihedral {' '.join(str(i + 1) for i in idx_phi)}\n"
                        f"dihedral {' '.join(str(i + 1) for i in idx_psi)}\n")
        mf = scf.RHF(make_mol(xyz, f"{name}_pyscf.log"))
        mf.conv_tol = 1e-9
        mol_eq = optimize(mf, constraints=str(cons), maxsteps=args.maxsteps,
                          prefix=str(out / f"{name}_geomeTRIC"))

    mf_eq = scf.RHF(mol_eq)
    mf_eq.conv_tol = 1e-10
    e_tot = mf_eq.kernel()
    if not mf_eq.converged:
        raise SystemExit(f"{name}: final SCF not converged")
    xyz_eq = mol_eq.atom_coords(unit="Angstrom")
    phi, psi = dihedral(*xyz_eq[idx_phi]), dihedral(*xyz_eq[idx_psi])
    write_xyz(out / f"{name}.xyz", elems, xyz_eq,
              f"{name} HF/6-31G* optimized E={e_tot:.8f} phi={phi:.2f} psi={psi:.2f}")

    pts = mk_points(elems, xyz_eq) * ANG2BOHR
    V = esp_at_points(mol_eq, mf_eq.make_rdm1(), pts)
    with open(out / f"{name}_vpot.out", "w") as fh:
        fh.write(f"{len(pts)}\n")
        for p, v in zip(pts, V):
            fh.write(f"{p[0]:16.8f} {p[1]:16.8f} {p[2]:16.8f} {v:18.10e}\n")
    print(f"{name}: E(HF/6-31G*) = {e_tot:.8f} Eh, phi = {phi:.2f}, psi = {psi:.2f}, "
          f"{len(pts)} ESP points -> {out}/{name}_vpot.out")


# ----------------------------------------------------------------------------
# RESP fit
# ----------------------------------------------------------------------------
def read_vpot(path):
    """ESP file: first line = number of points, then x y z V (bohr, a.u.)."""
    rows = [l.split() for l in Path(path).read_text().split("\n") if l.strip()]
    n = int(rows[0][0])
    data = np.array([[float(v) for v in r[:4]] for r in rows[1:1 + n]])
    return data[:, :3], data[:, 3]


def resp_fit(confs, n, fixed, groups, equiv, restrained, a, b=0.1, tol=1e-7, maxiter=200):
    """Restrained ESP fit (Bayly et al. 1993) over several conformers.

    confs       list of (atom coordinates [bohr], points [bohr], ESP [a.u.])
    fixed       {atom index: charge} kept constant
    groups      [(atom indices, total charge)] exact linear constraints
    equiv       [[atom indices]] free atoms forced to the same charge
    restrained  atom indices with the hyperbolic restraint a*(sqrt(q^2+b^2)-b)
    """
    var_of, nvar = {}, 0
    for cls in equiv:
        for i in cls:
            assert i not in fixed, "equivalenced atom is fixed"
            var_of[i] = nvar
        nvar += 1
    for i in range(n):
        if i not in fixed and i not in var_of:
            var_of[i] = nvar
            nvar += 1
    M = np.zeros((n, nvar))
    for i, v in var_of.items():
        M[i, v] = 1.0
    q0 = np.zeros(n)
    for i, q in fixed.items():
        q0[i] = q

    A = np.zeros((nvar, nvar))
    B = np.zeros(nvar)
    for R, P, V in confs:
        D = 1.0 / np.linalg.norm(P[:, None, :] - R[None, :, :], axis=2)
        DM = D @ M
        A += DM.T @ DM
        B += DM.T @ (V - D @ q0)
    C = np.array([M[idx].sum(axis=0) for idx, _ in groups])
    d = np.array([tot - q0[idx].sum() for idx, tot in groups])
    nc = len(groups)

    x = np.zeros(nvar)
    for _ in range(maxiter):
        q = M @ x + q0
        Ar = A.copy()
        for i in restrained:
            if i in var_of:
                Ar[var_of[i], var_of[i]] += a / np.sqrt(q[i] ** 2 + b ** 2)
        K = np.block([[Ar, C.T], [C, np.zeros((nc, nc))]])
        x_new = np.linalg.solve(K, np.concatenate([B, d]))[:nvar]
        if np.max(np.abs(x_new - x)) < tol:
            x = x_new
            break
        x = x_new
    return M @ x + q0


def rrms(confs, q):
    out = []
    for R, P, V in confs:
        D = 1.0 / np.linalg.norm(P[:, None, :] - R[None, :, :], axis=2)
        out.append(np.sqrt(np.sum((V - D @ q) ** 2) / np.sum(V ** 2)))
    return out


def two_stage_resp(atoms, confs, symmetric=()):
    """Standard RESP-A1 style fit for an Ace-X-NMe dipeptide.

    symmetric: pairs of residue atom names fitted with one common charge.
    """
    n = len(atoms)
    idx = {(a[0], a[1]): i for i, a in enumerate(atoms)}
    res_atoms = [i for i, a in enumerate(atoms) if a[0] == "RES"]
    groups = [([i for i, a in enumerate(atoms) if a[0] == r], 0.0) for r in ("ACE", "RES", "NME")]
    fixed = {idx["RES", k]: v for k, v in BACKBONE.items()}
    heavy = [i for i, a in enumerate(atoms) if a[2] != "H"]
    equiv1 = [[idx["RES", x], idx["RES", y]] for x, y in symmetric]

    q1 = resp_fit(confs, n, fixed, groups, equiv1, heavy, a=0.0005)

    refit = {
        "RES": ["CB", "HB1", "HB2"],
        "ACE": ["CH3", "HH31", "HH32", "HH33"],
        "NME": ["CH3", "HH31", "HH32", "HH33"],
    }
    refit_idx = {idx[r, nm] for r, names in refit.items() for nm in names}
    fixed2 = {i: q1[i] for i in range(n) if i not in refit_idx}
    equiv2 = [[idx["RES", "HB1"], idx["RES", "HB2"]],
              [idx["ACE", f"HH3{k}"] for k in (1, 2, 3)],
              [idx["NME", f"HH3{k}"] for k in (1, 2, 3)]]
    restr2 = [idx["RES", "CB"], idx["ACE", "CH3"], idx["NME", "CH3"]]
    q2 = resp_fit(confs, n, fixed2, groups, equiv2, restr2, a=0.001)

    charges = {atoms[i][1]: q2[i] for i in res_atoms}
    return charges, rrms(confs, q2)


def round_neutral(charges, names, sink="CA"):
    """Round to 4 decimals and put the rounding residue on one atom (sum exactly 0)."""
    r = {k: round(charges[k], 4) for k in names}
    r[sink] = round(r[sink] - round(sum(r.values()), 4), 4)
    assert abs(sum(r.values())) < 1e-9, sum(r.values())
    return r


def load_confs(meta, names):
    """Optimized geometry (<name>.xyz) and ESP (<name>_vpot.out) next to <name>.json."""
    confs = []
    for name in names:
        d = meta[name]["dir"]
        elems, xyz = read_xyz(d / f"{name}.xyz")
        assert elems == [a[2] for a in meta[name]["atoms"]], f"atom order changed in {name}"
        P, V = read_vpot(d / f"{name}_vpot.out")
        confs.append((xyz * ANG2BOHR, P, V))
    return confs


def write_rtp_hdb(resname, charges, ffdir):
    """<resname>.rtp and <resname>.hdb (lower case) for appending to a GROMACS amber force field."""
    d = RESIDUES[resname]
    lines = ["", f"[ {resname} ]", f" ; {d['title']}, RESP HF/6-31G*, see dopa/README.md", " [ atoms ]"]
    for k, (name, atype) in enumerate(d["types"], start=1):
        q = charges[name] if abs(charges[name]) >= 5e-6 else 0.0
        lines.append(f"  {name:>5s}  {atype:>4s}  {q:10.5f}  {k:3d}")
    lines.append(" [ bonds ]")
    lines += [f"  {i:>5s}  {j:>5s}" for i, j in d["bonds"]]
    lines.append(" [ impropers ]")
    lines += ["  " + "  ".join(f"{x:>5s}" for x in imp) for imp in d["impropers"]]
    Path(ffdir, f"{resname.lower()}.rtp").write_text("\n".join(lines) + "\n")
    Path(ffdir, f"{resname.lower()}.hdb").write_text("\n".join([f"{resname}\t{len(d['hdb'])}"] + d["hdb"]) + "\n")


def cmd_list(args):
    for name, (resname, _, _) in CONFORMERS.items():
        if not args.res or resname == args.res.upper():
            print(name)


def cmd_fit(args):
    meta = {}
    for f in sorted(Path(args.workdir).glob("*/*.json")):
        m = json.loads(f.read_text())
        if not isinstance(m, dict) or m.get("name") not in CONFORMERS:
            continue
        m["dir"] = f.parent
        meta[m["name"]] = m
    # fit every residue whose conformers are all present; Tyr is required
    residues = []
    for res in ["TYR"] + list(RESIDUES):
        names = [n for n, v in CONFORMERS.items() if v[0] == res]
        have = [n for n in names if n in meta]
        if len(have) == len(names):
            residues.append(res)
        elif have or res == "TYR":
            raise SystemExit(f"ERROR: {res} conformers missing in {args.workdir}: {sorted(set(names) - set(have))}")
    if len(residues) < 2:
        raise SystemExit(f"ERROR: no modified residue with complete QM data in {args.workdir}")
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    report = [f"HF/6-31G* (Cartesian d), MK shells {SHELLS} x vdW, {DENSITY} points/A^2, "
              "two-stage RESP a = 0.0005/0.001, b = 0.1", ""]

    fitted = {}
    for resname in residues:
        names = [k for k, v in CONFORMERS.items() if v[0] == resname]
        atoms = [tuple(a) for a in meta[names[0]]["atoms"]]
        confs = load_confs(meta, names)
        fitted[resname], err = two_stage_resp(atoms, confs, SYMMETRY[resname])
        report.append(f"== {resname}: conformers {', '.join(names)}")
        report += [f"   RRMS {nm}: {e:.4f}  ({len(c[1])} points)" for nm, e, c in zip(names, err, confs)]

    tyr = round_neutral(fitted["TYR"], list(TYR_FF99SB))
    report.append("\n== Check: TYR refit vs ff99SB")
    report.append(f"   {'atom':5s} {'fit':>8s} {'ff99SB':>8s} {'diff':>8s}")
    for k, ref in TYR_FF99SB.items():
        report.append(f"   {k:5s} {tyr[k]:8.4f} {ref:8.4f} {tyr[k] - ref:8.4f}")
    dev = max(abs(tyr[k] - v) for k, v in TYR_FF99SB.items())
    report.append(f"   max |diff| = {dev:.4f}")
    (out / "tyr_refit_charges.json").write_text(json.dumps(tyr, indent=1))

    for resname in residues[1:]:
        # Tyr-referenced charges: ff99SB Tyr + (fit - Tyr fit) for the atoms both
        # residues share; the new O/H from the fit. The protocol difference to ff94
        # cancels, so WT (ff99SB Tyr) and the modified residue differ only by the
        # modification. For each replaced Tyr H, (Tyr fit - ff99SB) of that H is
        # taken off its carbon, which makes the residue exactly neutral again.
        d = RESIDUES[resname]
        names = [n for n, _ in d["types"]]
        direct = round_neutral(fitted[resname], names)
        referenced = {k: (TYR_FF99SB[k] + fitted[resname][k] - fitted["TYR"][k]) if k in TYR_FF99SB
                      else fitted[resname][k] for k in names}
        for h, c in d["replaced"].items():
            referenced[c] -= fitted["TYR"][h] - TYR_FF99SB[h]
        sink = next(iter(d["replaced"].values()))
        referenced[sink] -= sum(referenced.values())          # numerical rest only
        referenced = round_neutral(referenced, names, sink=sink)
        used = referenced if args.scheme == "referenced" else direct
        tag = resname.lower()

        report.append(f"\n== {resname} ({d['title']}) charges, scheme '{args.scheme}' is written to "
                      f"{tag}.rtp (residue sum = {sum(used.values()):.4f})")
        report.append(f"   {'atom':5s} {'type':4s} {'direct':>8s} {'tyr-ref':>8s} {'TYR ff':>8s}")
        for name, atype in d["types"]:
            ref = TYR_FF99SB.get(name)
            refs = f"{ref:8.4f}" if ref is not None else f"{'-':>8s}"
            report.append(f"   {name:5s} {atype:4s} {direct[name]:8.4f} {referenced[name]:8.4f} {refs}")
        (out / f"{tag}_charges.json").write_text(json.dumps(used, indent=1))
        (out / f"{tag}_charges_direct.json").write_text(json.dumps(direct, indent=1))
        (out / f"{tag}_charges_tyr_referenced.json").write_text(json.dumps(referenced, indent=1))
        write_rtp_hdb(resname, used, args.ffdir)

    text = "\n".join(report)
    print(text)
    (out / "resp_report.txt").write_text(text + "\n")
    print(f"\nwrote {out}/resp_report.txt and rtp/hdb for {', '.join(residues[1:])} in {args.ffdir}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("list")
    p.add_argument("--res", help="only conformers of this residue (TYR, DAH, TOP)")
    p = sub.add_parser("build")
    p.add_argument("--outdir", required=True)
    p.add_argument("--only", choices=list(CONFORMERS), help="build a single conformer")
    p = sub.add_parser("qm")
    p.add_argument("--name", choices=list(CONFORMERS), required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--memory-mb", type=int, default=40000)
    p.add_argument("--maxsteps", type=int, default=300)
    p = sub.add_parser("fit")
    p.add_argument("--workdir", required=True, help="directory with one subdirectory per conformer")
    p.add_argument("--outdir", required=True)
    p.add_argument("--ffdir", required=True)
    p.add_argument("--scheme", choices=("referenced", "direct"), default="referenced",
                   help="charges written to the rtp files: Tyr-referenced (default) or direct RESP fit")
    args = ap.parse_args()
    {"list": cmd_list, "build": cmd_build, "qm": cmd_qm, "fit": cmd_fit}[args.cmd](args)


if __name__ == "__main__":
    sys.exit(main())
