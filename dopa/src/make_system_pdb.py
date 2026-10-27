#!/usr/bin/env python3
"""Cut a capped domain out of a protein model and optionally hydroxylate one Tyr.

  * keeps residues FIRST..LAST (heavy atoms only; pdb2gmx -ignh adds H)
  * residue FIRST-1 becomes an ACE cap (CA -> CH3, C, O)
  * residue LAST+1 becomes an NME cap (N, CA -> CH3)
  * residue numbers are kept (Y65 stays residue 65)

System names:
  wt                       unchanged
  y<N>_dah_ce1|ce2         Tyr N -> L-DOPA (DAH), new OH on the ring carbon named
                           CE1 / CE2 in the input
  y<N>_top_ce1|ce2         Tyr N -> 2,4,5-trihydroxyphenylalanine (TOPA, TOP):
                           DOPA plus an OH on the CD carbon para to the new 3-OH
  dah_ce1, dah_ce2         short for y<DEFAULT_MUT>_dah_ce1/ce2
  x<A>_<B>                 3,3'-dityrosine cross-link between Tyr A and Tyr B (both
                           residues -> DTY, bridge CE2(A)-CE2(B)); the side chains
                           get the chi1/chi2 rotamers that bring the two CE2 closest
                           with the best bridge alignment and no clashes; the link
                           data go to crosslink.json next to the output PDB
For ce1 the ring atom names are swapped (CD1<->CD2, CE1<->CE2) first, so the
residue templates (OH on CE2, and on CD1 for TOP) stay the same.

  make_system_pdb.py --parse SYSTEM --default-mut 65    -> "SITE RESNAME ORIENT"
"""
import argparse
import json
import re
from pathlib import Path

import numpy as np

RING = ("CG", "CD1", "CD2", "CE1", "CE2", "CZ")
C_OH = 1.364                          # parm99 CA-OH bond length (A)
EXTRA_OH = {"DAH": [("CE2", "OE2")], "TOP": [("CE2", "OE2"), ("CD1", "OD1")]}
ORDER = {
    "DAH": ["N", "CA", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "OE2", "CZ", "OH", "C", "O"],
    "TOP": ["N", "CA", "CB", "CG", "CD1", "OD1", "CD2", "CE1", "CE2", "OE2", "CZ", "OH", "C", "O"],
}


def parse_system(name, default_mut):
    """-> (site, resname, orient) or (None, None, None) for the wild type."""
    if name == "wt":
        return None, None, None
    m = re.fullmatch(r"dah_(ce1|ce2)", name)
    if m:
        return default_mut, "DAH", m.group(1)
    m = re.fullmatch(r"y(\d+)_(dah|top)_(ce1|ce2)", name)
    if m:
        return int(m.group(1)), m.group(2).upper(), m.group(3)
    m = re.fullmatch(r"x(\d+)_(\d+)", name)
    if m:
        return int(m.group(1)), "DTY", f"x{m.group(2)}"
    raise SystemExit(f"ERROR: unknown system name '{name}' (wt, dah_ce1, y112_dah_ce2, y65_top_ce1, x65_127, ...)")


# ---------------------------------------------------------------------------
# cross-link start geometry: chi1/chi2 search for the two tyrosines
# ---------------------------------------------------------------------------
SIDE = ["CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"]


def dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - (b0 @ b1) * b1
    w = b2 - (b2 @ b1) * b1
    return np.degrees(np.arctan2(np.cross(b1, v) @ w, v @ w))


def rotate(p, a, b, ang):
    k = (b - a) / np.linalg.norm(b - a)
    v = p - a
    c, s = np.cos(ang), np.sin(ang)
    return a + v * c + np.cross(k, v) * s + np.outer(v @ k, k) * (1 - c)


def rotamers(atoms, step=10):
    """All (chi1, chi2, side-chain coordinates in SIDE order) on a grid."""
    x1 = dihedral(atoms["N"], atoms["CA"], atoms["CB"], atoms["CG"])
    x2 = dihedral(atoms["CA"], atoms["CB"], atoms["CG"], atoms["CD1"])
    base = np.array([atoms[n] for n in SIDE])
    out = []
    for c1 in range(-180, 180, step):
        p1 = rotate(base, atoms["CA"], atoms["CB"], np.radians(c1 - x1))
        for c2 in range(-180, 180, step):
            p2 = p1.copy()
            p2[1:] = rotate(p1[1:], atoms["CB"], p1[0], np.radians(c2 - x2))
            out.append((c1, c2, p2))
    return out


def crosslink_rotamers(res, ra, rb, keep, min_clash=2.8, max_dist=4.0):
    """Rotamers of Tyr ra/rb that bring CE2(ra), CE2(rb) closest with the bridge pointing
    out of both rings, without clashes. Returns dict with chi values, coordinates, metrics."""
    env = [x for r in keep if r not in (ra, rb) for x in res[r]["atoms"].values()]
    env_a = np.array(env + [res[ra]["atoms"][n] for n in ("N", "C", "O")]
                     + [res[rb]["atoms"][n] for n in ("N", "CA", "C", "O", "CB")])
    env_b = np.array(env + [res[rb]["atoms"][n] for n in ("N", "C", "O")]
                     + [res[ra]["atoms"][n] for n in ("N", "CA", "C", "O", "CB")])
    ie, ring = SIDE.index("CE2"), [SIDE.index(n) for n in RING]

    def ok(rots, envx):
        return [(c1, c2, p) for c1, c2, p in rots
                if np.linalg.norm(p[1:, None] - envx[None], axis=2).min() > min_clash]

    best = None
    for c1a, c2a, p in ok(rotamers(res[ra]["atoms"]), env_a):
        ua = p[ie] - p[ring].mean(0)
        ua /= np.linalg.norm(ua)
        for c1b, c2b, q in ok(rotamers(res[rb]["atoms"]), env_b):
            v = q[ie] - p[ie]
            d = np.linalg.norm(v)
            if d > max_dist:
                continue
            dd = np.linalg.norm(p[:, None] - q[None], axis=2)
            dd[ie, ie] = 99.0
            if dd.min() < 2.4:
                continue
            ub = q[ie] - q[ring].mean(0)
            ub /= np.linalg.norm(ub)
            align = min(ua @ v / d, -(ub @ v) / d)
            key = (round(align, 2), -d)
            if best is None or key > best[0]:
                best = (key, dict(chi_a=[c1a, c2a], chi_b=[c1b, c2b], side_a=p, side_b=q,
                                  dist_A=float(d), align=float(align)))
    if best is None:
        raise SystemExit(f"ERROR: no clash-free rotamers bring CE2 of {ra} and {rb} within {max_dist} A")
    return best[1]


def read_pdb(path):
    residues = {}
    for line in open(path):
        if not line.startswith("ATOM"):
            continue
        name = line[12:16].strip()
        if name.startswith("H") or line[76:78].strip() == "H":
            continue
        resnr = int(line[22:26])
        res = residues.setdefault(resnr, {"resname": line[17:20].strip(), "atoms": {}})
        res["atoms"][name] = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return residues


def fmt(serial, name, resname, resnr, xyz):
    name = f" {name:<3s}" if len(name) < 4 else name
    elem = name.strip()[0]
    return (f"ATOM  {serial:5d} {name:4s} {resname:>3s} A{resnr:4d}    "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00          {elem:>2s}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--parse", metavar="SYSTEM", help="print SITE RESNAME ORIENT and exit")
    ap.add_argument("--input")
    ap.add_argument("--out")
    ap.add_argument("--first", type=int)
    ap.add_argument("--last", type=int)
    ap.add_argument("--default-mut", type=int, default=65)
    ap.add_argument("--system", help="system name, see above")
    args = ap.parse_args()

    if args.parse:
        site, resname, orient = parse_system(args.parse, args.default_mut)
        print(f"{site or 0} {resname or '-'} {orient or '-'}")
        return
    if not all([args.input, args.out, args.first, args.last, args.system]):
        ap.error("--input, --out, --first, --last and --system are required")
    site, newres, orient = parse_system(args.system, args.default_mut)

    res = read_pdb(args.input)
    for r in range(args.first - 1, args.last + 2):
        if r not in res:
            raise SystemExit(f"ERROR: residue {r} missing in {args.input}")
    if not {"CA", "C", "O"} <= res[args.first - 1]["atoms"].keys():
        raise SystemExit("ERROR: N-terminal neighbour lacks CA/C/O for the ACE cap")
    if not {"N", "CA"} <= res[args.last + 1]["atoms"].keys():
        raise SystemExit("ERROR: C-terminal neighbour lacks N/CA for the NME cap")

    if newres == "DTY":
        partner = int(orient[1:])
        for r in (site, partner):
            if res[r]["resname"] != "TYR":
                raise SystemExit(f"ERROR: residue {r} is {res[r]['resname']}, not TYR")
        link = crosslink_rotamers(res, site, partner, range(args.first, args.last + 1))
        for r, key in ((site, "side_a"), (partner, "side_b")):
            for n, x in zip(SIDE, link[key]):
                res[r]["atoms"][n] = x
            res[r]["resname"] = "DTY"
        meta = {"res_a": site, "res_b": partner, "atom": "CE2", "dist_A": round(link["dist_A"], 3),
                "align": round(link["align"], 3), "chi_a": link["chi_a"], "chi_b": link["chi_b"]}
        Path(args.out).with_name("crosslink.json").write_text(json.dumps(meta, indent=1))
        print(f"cross-link {site}-{partner}: rotamers chi1/chi2 {link['chi_a']} / {link['chi_b']}, "
              f"CE2-CE2 {link['dist_A']:.2f} A, bridge alignment {link['align']:.2f}")

    out = []
    prev = res[args.first - 1]["atoms"]
    out.append((args.first - 1, "ACE", [("CH3", prev["CA"]), ("C", prev["C"]), ("O", prev["O"])]))

    added = []
    for r in range(args.first, args.last + 1):
        resname = res[r]["resname"]
        atoms = dict(res[r]["atoms"])
        if r == site and newres != "DTY":
            if resname != "TYR":
                raise SystemExit(f"ERROR: residue {r} is {resname}, not TYR")
            if not set(RING) <= atoms.keys():
                raise SystemExit(f"ERROR: incomplete ring in residue {r}")
            if orient == "ce1":
                atoms["CD1"], atoms["CD2"] = atoms["CD2"], atoms["CD1"]
                atoms["CE1"], atoms["CE2"] = atoms["CE2"], atoms["CE1"]
            centroid = np.mean([atoms[a] for a in RING], axis=0)
            for carbon, oxygen in EXTRA_OH[newres]:
                u = atoms[carbon] - centroid
                atoms[oxygen] = atoms[carbon] + C_OH * u / np.linalg.norm(u)
                added.append(oxygen)
            atoms = {k: atoms[k] for k in ORDER[newres]}
            resname = newres
        out.append((r, resname, list(atoms.items())))

    nxt = res[args.last + 1]["atoms"]
    out.append((args.last + 1, "NME", [("N", nxt["N"]), ("CH3", nxt["CA"])]))

    serial = 0
    with open(args.out, "w") as fh:
        fh.write(f"REMARK  domain {args.first}-{args.last} of {args.input}, system {args.system}\n")
        for resnr, resname, atoms in out:
            for name, xyz in atoms:
                serial += 1
                fh.write(fmt(serial, name, resname, resnr, xyz))
        fh.write("TER\nEND\n")
    print(f"wrote {args.out}: {len(out)} residues (incl. caps), {serial} atoms, system {args.system}")

    # closest heavy atom of other residues to each added oxygen (start-structure clashes)
    for oxygen in added if newres != "DTY" else []:
        o = dict(next(a for r, n, a in out if r == site))[oxygen]
        best = min(((np.linalg.norm(x - o), f"{n}{r}:{an}") for r, n, a in out if r != site for an, x in a),
                   key=lambda t: t[0])
        print(f"  new {oxygen} of residue {site}: nearest other atom {best[1]} at {best[0]:.2f} A")


if __name__ == "__main__":
    main()
