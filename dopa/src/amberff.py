"""Minimal reader for a GROMACS amber force field directory and an Amber energy
function for small molecules (used to fit the dityrosine bridge torsion and to
write explicit parameters for the cross-link topology).

Units as in GROMACS: nm, degrees, kJ/mol.
"""
import itertools
from pathlib import Path

import numpy as np

COULOMB = 138.935458          # kJ mol^-1 nm e^-2


def _sections(path):
    """Yield (section, fields) for every data line of an itp file."""
    sec = None
    for line in Path(path).read_text().split("\n"):
        line = line.split(";")[0].strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("["):
            sec = line.strip("[] ").lower()
            continue
        yield sec, line.split()


class AmberFF:
    """Bond/angle/dihedral types and LJ parameters from ffnonbonded.itp + ffbonded.itp
    (plus optional extra files with the same sections)."""

    def __init__(self, ffdir, extra=()):
        self.lj, self.bonds, self.angles, self.dihedrals = {}, {}, {}, []
        files = [Path(ffdir) / "ffnonbonded.itp", Path(ffdir) / "ffbonded.itp"] + [Path(p) for p in extra]
        for f in files:
            for sec, fl in _sections(f):
                if sec == "atomtypes" and len(fl) >= 7:
                    self.lj[fl[0]] = (float(fl[-2]), float(fl[-1]))           # sigma, epsilon
                elif sec == "bondtypes" and len(fl) >= 5:
                    p = (float(fl[3]), float(fl[4]))
                    self.bonds[(fl[0], fl[1])] = self.bonds[(fl[1], fl[0])] = p
                elif sec == "angletypes" and len(fl) >= 6:
                    p = (float(fl[4]), float(fl[5]))
                    self.angles[(fl[0], fl[1], fl[2])] = self.angles[(fl[2], fl[1], fl[0])] = p
                elif sec == "dihedraltypes" and len(fl) >= 8:
                    self.dihedrals.append((tuple(fl[:4]), int(fl[4]), float(fl[5]), float(fl[6]), int(fl[7])))

    def bond(self, a, b):
        return self.bonds[(a, b)]

    def angle(self, a, b, c):
        return self.angles[(a, b, c)]

    def dihedral(self, types, func):
        """All (phase, k, mult) terms of the most specific matching dihedral type."""
        best, terms = None, []
        for t, f, ph, k, n in self.dihedrals:
            if f != func:
                continue
            for cand in (types, types[::-1]):
                if all(x == "X" or x == y for x, y in zip(t, cand)):
                    nwild = sum(x == "X" for x in t)
                    if best is None or nwild < best[0]:
                        best, terms = (nwild, t), [(ph, k, n)]
                    elif nwild == best[0] and t == best[1]:
                        terms.append((ph, k, n))
                    break
        if not terms:
            raise KeyError(f"no dihedral type {func} for {types}")
        return terms


# ---------------------------------------------------------------------------
# geometry helpers (coordinates in nm or A, angles in degrees)
# ---------------------------------------------------------------------------
def angle(a, b, c):
    u, v = a - b, c - b
    return np.degrees(np.arccos(np.clip(u @ v / np.linalg.norm(u) / np.linalg.norm(v), -1, 1)))


def dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - (b0 @ b1) * b1
    w = b2 - (b2 @ b1) * b1
    return np.degrees(np.arctan2(np.cross(b1, v) @ w, v @ w))


def dihedrals_vec(p0, p1, p2, p3):
    """Vectorized dihedral angles (rad) for arrays of points (n, 3)."""
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / np.linalg.norm(b1, axis=1)[:, None]
    v = b0 - np.sum(b0 * b1, axis=1)[:, None] * b1
    w = b2 - np.sum(b2 * b1, axis=1)[:, None] * b1
    return np.arctan2(np.sum(np.cross(b1, v) * w, axis=1), np.sum(v * w, axis=1))


class CompiledMM:
    """Vectorized Amber energy of a Molecule (x in nm, kJ/mol), with optional harmonic
    dihedral restraints [(i, j, k, l, target_deg, k_kJ_per_rad2)]."""

    def __init__(self, mol, ff, skip_proper=None, restraints=()):
        T = mol.types
        b = np.array([(i, j, *ff.bond(T[i], T[j])) for i, j in mol.bonds])
        self.bI, self.bJ, self.b0, self.kb = b[:, 0].astype(int), b[:, 1].astype(int), b[:, 2], b[:, 3]
        a = np.array([(i, j, k, *ff.angle(T[i], T[j], T[k])) for i, j, k in mol.angles])
        self.aI, self.aJ, self.aK = (a[:, c].astype(int) for c in range(3))
        self.th0, self.kt = np.radians(a[:, 3]), a[:, 4]
        d = []
        for q in mol.propers:
            t = tuple(T[x] for x in q)
            if skip_proper and skip_proper(t):
                continue
            d += [(*q, np.radians(ph), k, n) for ph, k, n in ff.dihedral(t, 9)]
        for q in mol.impropers:
            d += [(*q, np.radians(ph), k, n) for ph, k, n in ff.dihedral(tuple(T[x] for x in q), 4)]
        d = np.array(d)
        self.dI, self.dJ, self.dK, self.dL = (d[:, c].astype(int) for c in range(4))
        self.dph, self.dk, self.dn = d[:, 4], d[:, 5], d[:, 6]
        pairs14 = set(mol.pairs14)
        P = []
        n = len(T)
        for i in range(n):
            for j in range(i + 1, n):
                key = frozenset((i, j))
                if key in mol.excluded and key not in pairs14:
                    continue
                si, ei = ff.lj[T[i]]
                sj, ej = ff.lj[T[j]]
                s, eps = 0.5 * (si + sj), np.sqrt(ei * ej)
                slj, sqq = (0.5, 1 / 1.2) if key in pairs14 else (1.0, 1.0)
                P.append((i, j, slj * 4 * eps * s ** 12, slj * 4 * eps * s ** 6, sqq * COULOMB * mol.q[i] * mol.q[j]))
        P = np.array(P)
        self.pI, self.pJ, self.A, self.B, self.C = P[:, 0].astype(int), P[:, 1].astype(int), P[:, 2], P[:, 3], P[:, 4]
        r = np.array(restraints, float).reshape(-1, 6)
        self.rI, self.rJ, self.rK, self.rL = (r[:, c].astype(int) for c in range(4))
        self.rt, self.rk = np.radians(r[:, 4]), r[:, 5]

    def energy(self, x, restraints=True):
        x = x.reshape(-1, 3)
        e = 0.5 * np.sum(self.kb * (np.linalg.norm(x[self.bI] - x[self.bJ], axis=1) - self.b0) ** 2)
        u, v = x[self.aI] - x[self.aJ], x[self.aK] - x[self.aJ]
        c = np.sum(u * v, axis=1) / np.linalg.norm(u, axis=1) / np.linalg.norm(v, axis=1)
        e += 0.5 * np.sum(self.kt * (np.arccos(np.clip(c, -1, 1)) - self.th0) ** 2)
        phi = dihedrals_vec(x[self.dI], x[self.dJ], x[self.dK], x[self.dL])
        e += np.sum(self.dk * (1 + np.cos(self.dn * phi - self.dph)))
        r = np.linalg.norm(x[self.pI] - x[self.pJ], axis=1)
        e += np.sum(self.A / r ** 12 - self.B / r ** 6 + self.C / r)
        if restraints and len(self.rk):
            dphi = dihedrals_vec(x[self.rI], x[self.rJ], x[self.rK], x[self.rL]) - self.rt
            dphi = (dphi + np.pi) % (2 * np.pi) - np.pi
            e += 0.5 * np.sum(self.rk * dphi ** 2)
        return e

    def gradient(self, x, h=1e-6):
        x = x.ravel()
        g = np.empty_like(x)
        for k in range(len(x)):
            xp, xm = x.copy(), x.copy()
            xp[k] += h
            xm[k] -= h
            g[k] = (self.energy(xp) - self.energy(xm)) / (2 * h)
        return g


def minimize(mm, x0, gtol=1.0, maxiter=2000, m=10):
    """L-BFGS with backtracking line search; returns (x (n, 3), energy)."""
    x = x0.ravel().copy()
    fx, g = mm.energy(x), mm.gradient(x)
    S, Y = [], []
    for _ in range(maxiter):
        if np.max(np.abs(g)) < gtol:
            break
        q, alphas = g.copy(), []
        for s, y in reversed(list(zip(S, Y))):
            rho = 1.0 / (y @ s)
            a = rho * (s @ q)
            alphas.append((rho, a))
            q -= a * y
        r = q * ((S[-1] @ Y[-1]) / (Y[-1] @ Y[-1]) if S else 1e-5)
        for (s, y), (rho, a) in zip(zip(S, Y), reversed(alphas)):
            r += s * (a - rho * (y @ r))
        d = -r
        if d @ g >= 0:
            d, S, Y = -1e-5 * g, [], []
        t = 1.0
        while True:
            xn = x + t * d
            fn = mm.energy(xn)
            if fn <= fx + 1e-4 * t * (g @ d) or t < 1e-12:
                break
            t *= 0.5
        gn = gn_ = mm.gradient(xn)
        s, y = xn - x, gn_ - g
        if y @ s > 1e-14:
            S.append(s)
            Y.append(y)
            if len(S) > m:
                S.pop(0)
                Y.pop(0)
        x, fx, g = xn, fn, gn
    return x.reshape(-1, 3), fx


class Molecule:
    """Small molecule with Amber atom types, charges, bonds and explicit impropers.
    Angles, proper dihedrals, 1-4 pairs and exclusions follow from the bonds."""

    def __init__(self, types, charges, bonds, impropers=()):
        self.types, self.q = list(types), np.asarray(charges, float)
        self.bonds = [tuple(b) for b in bonds]
        self.impropers = [tuple(i) for i in impropers]
        n = len(self.types)
        nb = [set() for _ in range(n)]
        for i, j in self.bonds:
            nb[i].add(j)
            nb[j].add(i)
        self.angles = sorted({(i, j, k) if i < k else (k, j, i)
                              for j in range(n) for i, k in itertools.combinations(sorted(nb[j]), 2)})
        props = set()
        for j, k in self.bonds:
            for i in nb[j] - {k}:
                for l in nb[k] - {j}:
                    if i != l:
                        props.add((i, j, k, l) if (i, j) < (l, k) else (l, k, j, i))
        self.propers = sorted(props)
        excl = {frozenset(b) for b in self.bonds} | {frozenset((a[0], a[2])) for a in self.angles}
        self.pairs14 = sorted({frozenset((d[0], d[3])) for d in self.propers} - excl, key=sorted)
        self.excluded = excl | set(self.pairs14)

    def energy(self, ff, x, skip_proper=None):
        """Amber energy (kJ/mol) at coordinates x (nm). skip_proper(types) -> True drops a proper term."""
        e = 0.0
        for i, j in self.bonds:
            b0, k = ff.bond(self.types[i], self.types[j])
            e += 0.5 * k * (np.linalg.norm(x[i] - x[j]) - b0) ** 2
        for i, j, k in self.angles:
            th0, kt = ff.angle(self.types[i], self.types[j], self.types[k])
            e += 0.5 * kt * np.radians(angle(x[i], x[j], x[k]) - th0) ** 2
        for d in self.propers:
            t = tuple(self.types[a] for a in d)
            if skip_proper and skip_proper(t):
                continue
            phi = dihedral(*x[list(d)])
            for ph, k, n in ff.dihedral(t, 9):
                e += k * (1 + np.cos(np.radians(n * phi - ph)))
        for d in self.impropers:
            t = tuple(self.types[a] for a in d)
            phi = dihedral(*x[list(d)])
            for ph, k, n in ff.dihedral(t, 4):
                e += k * (1 + np.cos(np.radians(n * phi - ph)))
        n = len(self.types)
        pairs14 = set(self.pairs14)
        for i in range(n):
            for j in range(i + 1, n):
                key = frozenset((i, j))
                if key in self.excluded and key not in pairs14:
                    continue
                r = np.linalg.norm(x[i] - x[j])
                si, ei = ff.lj[self.types[i]]
                sj, ej = ff.lj[self.types[j]]
                s, eps = 0.5 * (si + sj), np.sqrt(ei * ej)
                lj = 4 * eps * ((s / r) ** 12 - (s / r) ** 6)
                qq = COULOMB * self.q[i] * self.q[j] / r
                if key in pairs14:
                    lj, qq = 0.5 * lj, qq / 1.2
                e += lj + qq
        return e
