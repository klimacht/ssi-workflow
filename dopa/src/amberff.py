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
