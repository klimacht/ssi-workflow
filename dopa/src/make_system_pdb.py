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
For ce1 the ring atom names are swapped (CD1<->CD2, CE1<->CE2) first, so the
residue templates (OH on CE2, and on CD1 for TOP) stay the same.

  make_system_pdb.py --parse SYSTEM --default-mut 65    -> "SITE RESNAME ORIENT"
"""
import argparse
import re
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
    raise SystemExit(f"ERROR: unknown system name '{name}' (wt, dah_ce1, y112_dah_ce2, y65_top_ce1, ...)")


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

    out = []
    prev = res[args.first - 1]["atoms"]
    out.append((args.first - 1, "ACE", [("CH3", prev["CA"]), ("C", prev["C"]), ("O", prev["O"])]))

    added = []
    for r in range(args.first, args.last + 1):
        resname = res[r]["resname"]
        atoms = dict(res[r]["atoms"])
        if r == site:
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
    for oxygen in added:
        o = dict(next(a for r, n, a in out if r == site))[oxygen]
        best = min(((np.linalg.norm(x - o), f"{n}{r}:{an}") for r, n, a in out if r != site for an, x in a),
                   key=lambda t: t[0])
        print(f"  new {oxygen} of residue {site}: nearest other atom {best[1]} at {best[0]:.2f} A")


if __name__ == "__main__":
    main()
