#!/usr/bin/env python3
"""Cut a capped domain out of a protein model and optionally turn a Tyr into L-DOPA.

  * keeps residues FIRST..LAST (heavy atoms only; pdb2gmx -ignh adds H)
  * residue FIRST-1 becomes an ACE cap (CA -> CH3, C, O)
  * residue LAST+1 becomes an NME cap (N, CA -> CH3)
  * variant wt       : unchanged
    variant dah_ce2  : Tyr MUT -> DAH, 3-OH (OE2) on the ring carbon named CE2 in the input
    variant dah_ce1  : same, but on the carbon named CE1 (ring atom names swapped
                       CD1<->CD2, CE1<->CE2 first, so the DAH template stays the same)
Residue numbers are kept, so Y65 stays residue 65 in all systems.
"""
import argparse
import numpy as np

RING = ("CG", "CD1", "CD2", "CE1", "CE2", "CZ")
C_OH = 1.364                          # parm99 CA-OH bond length (A)


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
    ap.add_argument("--input", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--first", type=int, required=True)
    ap.add_argument("--last", type=int, required=True)
    ap.add_argument("--mut", type=int, required=True)
    ap.add_argument("--variant", choices=("wt", "dah_ce1", "dah_ce2"), required=True)
    args = ap.parse_args()

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

    for r in range(args.first, args.last + 1):
        resname = res[r]["resname"]
        atoms = dict(res[r]["atoms"])
        if r == args.mut and args.variant != "wt":
            if resname != "TYR":
                raise SystemExit(f"ERROR: residue {r} is {resname}, not TYR")
            if not set(RING) <= atoms.keys():
                raise SystemExit(f"ERROR: incomplete ring in residue {r}")
            if args.variant == "dah_ce1":
                atoms["CD1"], atoms["CD2"] = atoms["CD2"], atoms["CD1"]
                atoms["CE1"], atoms["CE2"] = atoms["CE2"], atoms["CE1"]
            centroid = np.mean([atoms[a] for a in RING], axis=0)
            u = atoms["CE2"] - centroid
            oe2 = atoms["CE2"] + C_OH * u / np.linalg.norm(u)
            order = ["N", "CA", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "OE2", "CZ", "OH", "C", "O"]
            atoms["OE2"] = oe2
            atoms = {k: atoms[k] for k in order}
            resname = "DAH"
        out.append((r, resname, list(atoms.items())))

    nxt = res[args.last + 1]["atoms"]
    out.append((args.last + 1, "NME", [("N", nxt["N"]), ("CH3", nxt["CA"])]))

    serial = 0
    with open(args.out, "w") as fh:
        fh.write(f"REMARK  domain {args.first}-{args.last} of {args.input}, variant {args.variant}\n")
        for resnr, resname, atoms in out:
            for name, xyz in atoms:
                serial += 1
                fh.write(fmt(serial, name, resname, resnr, xyz))
        fh.write("TER\nEND\n")
    n_res = len(out)
    print(f"wrote {args.out}: {n_res} residues (incl. caps), {serial} atoms, variant {args.variant}")


if __name__ == "__main__":
    main()
