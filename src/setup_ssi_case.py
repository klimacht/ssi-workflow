#!/usr/bin/env python3
import argparse
import json
import os
import subprocess
from pathlib import Path

import numpy as np
import MDAnalysis as mda


def run(cmd):
    print("RUN:", " ".join(str(x) for x in cmd), flush=True)
    subprocess.check_call([str(x) for x in cmd])


def phi_dir_name(phi):
    whole = int(phi)
    frac = int(round((phi - whole) * 10))
    return f"phi_{whole:02d}p{frac}"


def compress_indices(indices):
    """
    Compress sorted 1-based atom indices into ranges.
    For water oxygens in GROMACS this should usually become:
      first-last:3
    """
    indices = sorted(int(x) for x in indices)
    if not indices:
        raise RuntimeError("No indices supplied to compress_indices.")

    runs = []
    i = 0
    n = len(indices)

    while i < n:
        start = indices[i]

        if i == n - 1:
            runs.append((start, start, 1))
            break

        step = indices[i + 1] - indices[i]
        j = i + 1

        while j + 1 < n and indices[j + 1] - indices[j] == step:
            j += 1

        end = indices[j]
        runs.append((start, end, step))
        i = j + 1

    parts = []
    for start, end, step in runs:
        if start == end:
            parts.append(str(start))
        elif step == 1:
            parts.append(f"{start}-{end}")
        else:
            parts.append(f"{start}-{end}:{step}")

    return " ".join(parts)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdbid", required=True)
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--work-dir", required=True)
    ap.add_argument("--gmx", required=True)
    ap.add_argument("--phi", required=True, type=float)
    ap.add_argument("--radius-nm", required=True, type=float)
    ap.add_argument("--cell-size-nm", required=True, type=float)
    ap.add_argument("--union-mode", required=True, choices=["product", "nearest"])
    ap.add_argument("--bias-mode", required=True, choices=["plumed-restraint", "internal-indus"])
    ap.add_argument("--mdp", required=True)
    args = ap.parse_args()

    run_dir = Path(args.run_dir)
    work_dir = Path(args.work_dir)
    unbiased_dir = run_dir / "unbiased"

    work_dir.mkdir(parents=True, exist_ok=True)
    os.chdir(work_dir)

    required = [
        unbiased_dir / "prod_unbiased.gro",
        unbiased_dir / "prod_unbiased.cpt",
        unbiased_dir / "topol.top",

    ]

    for f in required:
        if not f.is_file() or f.stat().st_size == 0:
            raise RuntimeError(f"Missing or empty required file: {f}")

    # Copy topology files locally; handle single or multi-chain posre files.
    posre_files = sorted(unbiased_dir.glob("posre*.itp"))
    if not posre_files:
        raise RuntimeError(f"No posre*.itp files found in {unbiased_dir}")
    topol_chain_files = sorted(unbiased_dir.glob("topol_*.itp"))
    for src in [unbiased_dir / "topol.top"] + posre_files + topol_chain_files:
        (work_dir / src.name).write_bytes(src.read_bytes())

    # Use symlinks for large coordinate/checkpoint inputs.
    for src_name in ["prod_unbiased.gro", "prod_unbiased.cpt"]:
        src = unbiased_dir / src_name
        dst = work_dir / src_name
        if dst.exists() or dst.is_symlink():
            dst.unlink()
        dst.symlink_to(src)

    # Generate biased TPR.
    run([
        "mpirun", "-np", "1", args.gmx, "grompp",
        "-f", args.mdp,
        "-c", "prod_unbiased.gro",
        "-t", "prod_unbiased.cpt",
        "-r", "prod_unbiased.gro",
        "-p", "topol.top",
        "-o", "ssi_phi.tpr",
        "-maxwarn", "1",
    ])

    if not Path("ssi_phi.tpr").is_file() or Path("ssi_phi.tpr").stat().st_size == 0:
        raise RuntimeError("ssi_phi.tpr was not created.")

    # MDAnalysis may not expose box dimensions from GROMACS TPR files.
    # Therefore, use the final unbiased GRO file for coordinates, atom indices,
    # atom names/residue names, and box dimensions. GRO/PDB-like coordinates in
    # MDAnalysis are in Angstrom, so convert Angstrom -> nm by multiplying by 0.1.
    u = mda.Universe("prod_unbiased.gro")

    if u.dimensions is None:
        raise RuntimeError("MDAnalysis could not read box dimensions from prod_unbiased.gro")

    protein_heavy = u.select_atoms("protein and not name H*")
    if len(protein_heavy) == 0:
        raise RuntimeError("No protein heavy atoms found.")

    water_oxygen = u.select_atoms(
        "(resname SOL or resname WAT or resname HOH) and "
        "(name OW or name OH2 or name O)"
    )
    if len(water_oxygen) == 0:
        raise RuntimeError("No water oxygen atoms found.")

    centers_nm = protein_heavy.positions * 0.1
    np.savetxt("heavy_atom_centers_fixed.dat", centers_nm, fmt="%.6f")

    box_nm = np.asarray(u.dimensions[:3], dtype=float) * 0.1
    if np.any(box_nm <= 0.0):
        raise RuntimeError(f"Invalid box dimensions read from prod_unbiased.gro: {box_nm}")

    # PLUMED/INDUS atom_index is treated as 1-based.
    water_oxygen_indices_1based = [int(a.index) + 1 for a in water_oxygen]
    target_expression = compress_indices(water_oxygen_indices_1based)

    # Write INDUS input.
    indus_lines = []
    indus_lines.append("# Autogenerated INDUS input for SSI hydrophobicity mapping")
    indus_lines.append("# Units: lengths in nm, phi in kJ mol^-1 particle^-1")
    indus_lines.append("")
    indus_lines.append(f"Target = [ atom_index {target_expression} ]")
    indus_lines.append("")
    indus_lines.append("ProbeVolume = {")
    indus_lines.append("  type          = union_spheres")
    indus_lines.append("  centers_file  = heavy_atom_centers_fixed.dat")
    indus_lines.append(f"  box_lengths   = [ {box_nm[0]:.5f} {box_nm[1]:.5f} {box_nm[2]:.5f} ]")
    indus_lines.append(f"  r_max         = {args.radius_nm:.5f}")
    indus_lines.append(f"  cell_size     = {args.cell_size_nm:.5f}")
    indus_lines.append(f"  union_mode    = {args.union_mode}")
    indus_lines.append("")
    indus_lines.append("  sigma         = 0.01")
    indus_lines.append("  alpha_c       = 0.02")
    indus_lines.append("}")
    indus_lines.append("")

    if args.bias_mode == "internal-indus":
        indus_lines.append("Bias = {")
        indus_lines.append("  order_parameter = ntilde")
        indus_lines.append(f"  phi = {args.phi:.6f}")
        indus_lines.append("}")
        indus_lines.append("")

    Path("indus_union_runtime.input").write_text("\n".join(indus_lines))

    # Write PLUMED input.
    if args.bias_mode == "plumed-restraint":
        plumed = f"""# Autogenerated PLUMED input for SSI / INDUS
# Bias mode: PLUMED RESTRAINT
# Linear unfavorable potential U = phi * Ntilde
# No internal INDUS Bias block is used in this mode.

indus: INDUS INPUTFILE=indus_union_runtime.input

restraint: RESTRAINT ARG=indus.ntilde AT=0.0 KAPPA=0.0 SLOPE={args.phi:.6f}

PRINT ...
  ARG=indus.n,indus.ntilde,restraint.bias
  STRIDE=500
  FILE=COLVAR
... PRINT
"""
    else:
        plumed = """# Autogenerated PLUMED input for SSI / INDUS
# Bias mode: internal INDUS bias
# Do not add a PLUMED RESTRAINT in this mode.

indus: INDUS INPUTFILE=indus_union_runtime.input NO_SHARE_ALL_DERIVATIVES

PRINT ...
  ARG=indus.n,indus.ntilde,indus.ubias
  STRIDE=500
  FILE=COLVAR
... PRINT
"""
    Path("plumed_runtime.dat").write_text(plumed)

    # Summary for inspection.
    summary = {
        "pdbid": args.pdbid,
        "work_dir": str(work_dir),
        "phi": args.phi,
        "radius_nm": args.radius_nm,
        "cell_size_nm": args.cell_size_nm,
        "union_mode": args.union_mode,
        "bias_mode": args.bias_mode,
        "n_atoms_total": int(u.atoms.n_atoms),
        "n_protein_atoms": int(len(u.select_atoms("protein"))),
        "n_protein_heavy_atoms": int(len(protein_heavy)),
        "n_water_oxygen_atoms": int(len(water_oxygen)),
        "first_water_oxygen_index_1based": int(water_oxygen_indices_1based[0]),
        "last_water_oxygen_index_1based": int(water_oxygen_indices_1based[-1]),
        "target_expression": target_expression,
        "box_nm": [float(x) for x in box_nm],
        "centers_file": "heavy_atom_centers_fixed.dat",
        "indus_input": "indus_union_runtime.input",
        "plumed_input": "plumed_runtime.dat",
    }

    Path("ssi_setup_summary.json").write_text(json.dumps(summary, indent=2))

    with open("ssi_setup_summary.txt", "w") as fh:
        fh.write("SSI setup summary\n")
        fh.write("=================\n")
        for k, v in summary.items():
            fh.write(f"{k}: {v}\n")

        fh.write("\nFirst 5 heavy-atom centers in nm:\n")
        for row in centers_nm[:5]:
            fh.write(f"{row[0]:.6f} {row[1]:.6f} {row[2]:.6f}\n")

        fh.write("\nLast 5 heavy-atom centers in nm:\n")
        for row in centers_nm[-5:]:
            fh.write(f"{row[0]:.6f} {row[1]:.6f} {row[2]:.6f}\n")

    print("============================================================")
    print("SSI setup completed")
    print("Summary:")
    print(Path("ssi_setup_summary.txt").read_text())
    print("============================================================")


if __name__ == "__main__":
    main()
