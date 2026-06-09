"""
render_maps_pymol.py — STEP 1: render raw PyMOL images for all proteins/maps.

Run via:
  source /mpcdf/soft/RHEL_9/packages/x86_64/anaconda/3/2023.03/etc/profile.d/conda.sh
  conda activate /u/tklimach/conda-envs/pymol_env
  pymol -cq render_maps_pymol.py

Then add colorbars in step 2:
  module purge
  module load python-waterboa/2025.06
  python3 render_maps_compose.py
"""

import os
import glob
import csv
from pathlib import Path

ROOT = Path("/ptmp/tklimach/ssi_workflow")
RAWDIR = Path("/ptmp/tklimach/ssi-workflow/examples/figures/_raw")
RAWDIR.mkdir(parents=True, exist_ok=True)

PROTEINS = ["1DPX", "1HRC", "1FS3", "1RBB", "3LDJ", "4IBA", "PHLP5"]
RENDER_WIDTH = 500
RENDER_HEIGHT = 500
PHISTAR_RANGE = (-3.0, 3.0)
ETA_RANGE = (0.0, 1.0)
VIEWS = [("front", (0,0,0)), ("side", (0,90,0)), ("back", (0,180,0))]


def find_eta_pdb(protein):
    pat = str(ROOT / "runs" / protein / "dewetting" / "dewetting_phi_*.pdb")
    files = [f for f in glob.glob(pat) if "phistar" not in f]
    tsv = ROOT / "runs" / protein / "analysis" / "N_vs_phi.tsv"
    phistar_dir = None
    if tsv.is_file():
        with open(tsv) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                if r.get("phi_star", "0") == "1":
                    phistar_dir = r["dirname"]; break
    if phistar_dir:
        for f in files:
            if phistar_dir in f:
                return f, phistar_dir
    return (files[0] if files else None, phistar_dir or "phi_05p0")


def find_phistar_pdb(protein):
    f = ROOT / "runs" / protein / "dewetting" / "dewetting_phistar.pdb"
    return str(f) if f.is_file() else None


def render_pymol(pdb_path, b_range, is_phistar):
    cmd.delete("all")
    cmd.load(pdb_path, "prot")
    if is_phistar:
        cmd.hide("everything", "prot")
        cmd.show("surface", "prot and b > -7")
    else:
        cmd.hide("everything", "prot")
        cmd.show("surface", "prot and b > -0.5")
    lo, hi = b_range
    cmd.spectrum("b", "red_white_blue", "prot", minimum=lo, maximum=hi)
    cmd.bg_color("white")
    cmd.set("ray_shadows", 0)
    cmd.set("ray_opaque_background", 0)
    cmd.set("surface_quality", 1)
    cmd.set("ambient", 0.3)


def render_view(out_path, rx, ry, rz, base_view):
    cmd.set_view(base_view)
    if rx: cmd.rotate("x", rx, "prot")
    if ry: cmd.rotate("y", ry, "prot")
    if rz: cmd.rotate("z", rz, "prot")
    cmd.png(str(out_path), width=RENDER_WIDTH, height=RENDER_HEIGHT, dpi=100, ray=1)


def main():
    for protein in PROTEINS:
        print(f"=== {protein} ===")

        eta_pdb, phi_star_dir = find_eta_pdb(protein)
        if eta_pdb:
            render_pymol(eta_pdb, ETA_RANGE, is_phistar=False)
            cmd.orient("prot"); cmd.zoom("prot", buffer=2.0)
            base = cmd.get_view()
            for name, (rx, ry, rz) in VIEWS:
                out = RAWDIR / f"{protein}_eta_{name}.png"
                render_view(out, rx, ry, rz, base)
                print(f"  wrote {out.name}")

        phistar_pdb = find_phistar_pdb(protein)
        if phistar_pdb:
            render_pymol(phistar_pdb, PHISTAR_RANGE, is_phistar=True)
            cmd.orient("prot"); cmd.zoom("prot", buffer=2.0)
            base = cmd.get_view()
            for name, (rx, ry, rz) in VIEWS:
                out = RAWDIR / f"{protein}_phistar_{name}.png"
                render_view(out, rx, ry, rz, base)
                print(f"  wrote {out.name}")

    print(f"\nRaw images in {RAWDIR}")


main()
