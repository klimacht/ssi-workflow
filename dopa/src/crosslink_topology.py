#!/usr/bin/env python3
"""Prepare the topology of a 3,3'-dityrosine cross-link built by pdb2gmx (specbond.dat).

From topol.top (pdb2gmx + solvate + genion) and the solvated coordinates:
  topol.top        + the two improper dihedrals that keep the bridge in each ring
                     plane (pdb2gmx only takes impropers from the residue templates)
  topol_grow.top   every bonded term that spans the link gets explicit A/B
                     parameters: bond length from the start distance (A) to the
                     equilibrium length (B), angles/dihedrals/impropers switched on
                     from k = 0 (A) to full strength (B); with free-energy = yes and
                     delta-lambda the link closes slowly during MD (grow.mdp)
  posre_grow.itp   position restraints without the residues near the link, so the
                     cluster can rearrange while the bond forms (define -DPOSRES_GROW)

  crosslink_topology.py prepare --top topol.top --gro solv_ions.gro --link crosslink.json --ffdir FF
  crosslink_topology.py report  --gro em.gro --link crosslink.json      (bridge geometry)
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from amberff import AmberFF, angle, dihedral          # noqa: E402

RING_SC = ("CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH")


def read_gro(path):
    lines = Path(path).read_text().split("\n")
    n = int(lines[1])
    atoms = []
    for l in lines[2:2 + n]:
        atoms.append((int(l[0:5]), l[5:10].strip(), l[10:15].strip(),
                      np.array([float(l[20:28]), float(l[28:36]), float(l[36:44])])))
    return atoms


def gro_index(atoms, resnr, name):
    return next(i for i, a in enumerate(atoms) if a[0] == resnr and a[2] == name)


def link_atoms(link):
    return (link["res_a"], link["atom"]), (link["res_b"], link["atom"])


def cmd_prepare(args):
    link = json.loads(Path(args.link).read_text())
    (ra, na), (rb, nb) = link_atoms(link)
    top_lines = Path(args.top).read_text().split("\n")
    gro = read_gro(args.gro)
    ff = AmberFF(args.ffdir)

    # --- first moleculetype: atoms and bonded sections
    atoms, sec, mol = {}, None, 0
    for line in top_lines:
        s = line.split(";")[0].strip()
        if s.startswith("[") and "moleculetype" in s:
            mol += 1
            if mol == 2:
                break
        if s.startswith("["):
            sec = s.strip("[] ").lower()
            continue
        if mol == 1 and sec == "atoms" and s and not s.startswith("#"):
            f = s.split()
            atoms[int(f[0])] = (f[1], int(f[2]), f[3], f[4])
    by_name = {(v[1], v[3]): nr for nr, v in atoms.items()}
    ia, ib = by_name[(ra, na)], by_name[(rb, nb)]
    resn = lambda nr: atoms[nr][1]

    def spans(idx):
        return {ra, rb} <= {resn(i) for i in idx}

    grow, final = [], []
    sec, mol, n_span, found_bond, last_dih = None, 0, 0, False, None
    for k, line in enumerate(top_lines):
        s = line.split(";")[0].strip()
        if s.startswith("[") and "moleculetype" in s:
            mol += 1
        if s.startswith("["):
            sec = s.strip("[] ").lower()
        elif mol == 1 and s and not s.startswith("#") and sec in ("bonds", "angles", "dihedrals"):
            f = s.split()
            na_ = {"bonds": 2, "angles": 3, "dihedrals": 4}[sec]
            idx, func = [int(x) for x in f[:na_]], int(f[na_])
            if sec == "dihedrals":
                last_dih = len(final)
            if spans(idx):
                n_span += 1
                types = tuple(atoms[i][0] for i in idx)
                lead = " ".join(f"{i:5d}" for i in idx)
                if sec == "bonds":
                    found_bond |= set(idx) == {ia, ib}
                    b0, kb = ff.bond(*types)
                    d0 = np.linalg.norm(gro[ia - 1][3] - gro[ib - 1][3])
                    grow.append(f"{lead} {func} {d0:.5f} {kb:.1f} {b0:.5f} {kb:.1f} ; cross-link, grows")
                elif sec == "angles":
                    th, kt = ff.angle(*types)
                    grow.append(f"{lead} {func} {th:.3f} 0.0 {th:.3f} {kt:.3f} ; cross-link")
                else:
                    for ph, kd, n in ff.dihedral(types, func):
                        grow.append(f"{lead} {func} {ph:.2f} 0.0 {n} {ph:.2f} {kd:.5f} {n} ; cross-link")
                final.append(line)
                continue
        grow.append(line)
        final.append(line)
    if not found_bond:
        raise SystemExit(f"ERROR: no bond between {ra}:{na} and {rb}:{nb} in {args.top} (specbond.dat?)")

    # --- bridge impropers (rtp templates cannot reference the other residue)
    imp = []
    for r, o in ((ra, rb), (rb, ra)):
        q = [by_name[(r, "CZ")], by_name[(r, "CD2")], by_name[(r, "CE2")], by_name[(o, "CE2")]]
        types = tuple(atoms[i][0] for i in q)
        lead = " ".join(f"{i:5d}" for i in q)
        imp.append((f"{lead} 4 ; cross-link bridge improper",
                    [f"{lead} 4 {ph:.2f} 0.0 {n} {ph:.2f} {kd:.5f} {n} ; cross-link bridge improper"
                     for ph, kd, n in ff.dihedral(types, 4)]))
    pos_final = last_dih + 1
    final[pos_final:pos_final] = [x[0] for x in imp]
    # in topol_grow.top the impropers go after the last dihedral line of the protein as well
    grow_dih_end = None
    sec, mol = None, 0
    for i, l in enumerate(grow):
        s = l.split(";")[0].strip()
        if s.startswith("[") and "moleculetype" in s:
            mol += 1
        if s.startswith("["):
            sec = s.strip("[] ").lower()
        elif mol == 1 and sec == "dihedrals" and s and not s.startswith("#"):
            grow_dih_end = i + 1
    grow[grow_dih_end:grow_dih_end] = [y for x in imp for y in x[1]]

    # --- posre_grow include after the POSRES block
    for i, l in enumerate(grow):
        if '#include "posre.itp"' in l:
            j = next(k for k in range(i, len(grow)) if grow[k].strip().startswith("#endif"))
            grow[j + 1:j + 1] = ["", "; position restraints away from the cross-link (growth)",
                                 "#ifdef POSRES_GROW", '#include "posre_grow.itp"', "#endif"]
            break
    else:
        raise SystemExit("ERROR: no posre.itp include in the topology")

    # --- posre_grow.itp: drop residues near either ring
    ring = np.array([gro[by_name[(r, n)] - 1][3] for r in (ra, rb) for n in RING_SC])
    prot = [a for a in gro[:len(atoms)]]
    near = {a[0] for a in prot if np.linalg.norm(ring - a[3], axis=1).min() < args.radius}
    out, kept, dropped = [], 0, 0
    for l in Path(args.posre).read_text().split("\n"):
        f = l.split(";")[0].split()
        if len(f) >= 5 and f[0].isdigit():
            if resn(int(f[0])) in near:
                dropped += 1
                continue
            kept += 1
        out.append(l)
    Path(args.top).with_name("posre_grow.itp").write_text("\n".join(out))

    Path(args.top).with_name("topol_nolinkimp.top").write_text("\n".join(top_lines))
    Path(args.top).write_text("\n".join(final))
    Path(args.top).with_name("topol_grow.top").write_text("\n".join(grow))
    d0 = np.linalg.norm(gro[ia - 1][3] - gro[ib - 1][3])
    print(f"cross-link {ra}:{na} (atom {ia}) - {rb}:{nb} (atom {ib}): start distance {d0:.3f} nm, "
          f"{n_span} bonded terms span the link (+2 bridge impropers)")
    print(f"posre_grow.itp: {kept} restrained atoms, {dropped} released "
          f"(residues within {args.radius} nm of the rings: {sorted(near)})")


def cmd_report(args):
    link = json.loads(Path(args.link).read_text())
    (ra, na), (rb, nb) = link_atoms(link)
    gro = read_gro(args.gro)
    x = lambda r, n: gro[gro_index(gro, r, n)][3]
    d = np.linalg.norm(x(ra, na) - x(rb, nb))
    angs = [angle(x(r, nn), x(r, "CE2"), x(o, "CE2")) for r, o in ((ra, rb), (rb, ra)) for nn in ("CD2", "CZ")]
    phi = dihedral(x(ra, "CD2"), x(ra, "CE2"), x(rb, "CE2"), x(rb, "CD2"))
    print(f"bridge {ra}:{na}-{rb}:{nb}: {d:.4f} nm; angles CD2/CZ-CE2-CE2' "
          + " ".join(f"{a:.1f}" for a in angs) + f" deg; torsion CD2-CE2-CE2'-CD2' {phi:.1f} deg")
    if abs(d - args.expect) > 0.01:
        raise SystemExit(f"ERROR: bridge length {d:.4f} nm differs from {args.expect} nm by more than 0.01 nm")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("prepare")
    p.add_argument("--top", required=True)
    p.add_argument("--gro", required=True)
    p.add_argument("--link", required=True)
    p.add_argument("--ffdir", required=True)
    p.add_argument("--posre", default="posre.itp")
    p.add_argument("--radius", type=float, default=0.6, help="nm around the rings released during growth")
    p = sub.add_parser("report")
    p.add_argument("--gro", required=True)
    p.add_argument("--link", required=True)
    p.add_argument("--expect", type=float, default=0.149)
    args = ap.parse_args()
    {"prepare": cmd_prepare, "report": cmd_report}[args.cmd](args)


if __name__ == "__main__":
    main()
