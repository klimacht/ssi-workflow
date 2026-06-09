"""
render_maps.py — render dewetting maps as PNG via PyMOL.

For each protein, produces three orthogonal views of:
  - the single-window eta map  (dewetting_phi_05p0.pdb or similar)
  - the per-atom phi*_i map    (dewetting_phistar.pdb)

Run via:
  module purge
  module load anaconda/3/2023.03
  conda activate /u/tklimach/conda-envs/pymol_env
  python render_maps.py

Output: examples/figures/<PDBID>_<map>_<view>.png
"""

import os
import glob
from pathlib import Path

# ---- config ----
ROOT = Path("/ptmp/tklimach/ssi_workflow")
OUTDIR = Path("/ptmp/tklimach/ssi-workflow/examples/figures")
OUTDIR.mkdir(parents=True, exist_ok=True)

PROTEINS = ["1DPX", "1HRC", "1FS3", "1RBB", "3LDJ", "4IBA", "PHLP5"]

# Per-map color range. phi*_i map is centered around 0 (median mode), eta is 0..1.
PHISTAR_RANGE = (-3.0, 3.0)      # symmetric around 0
ETA_RANGE = (0.0, 1.0)

# Three orthogonal rotations
VIEWS = [
    ("front",  ( 0,   0, 0)),
    ("side",   ( 0,  90, 0)),
    ("back",   ( 0, 180, 0)),
]


def find_eta_pdb(protein):
    """Find the single-window eta PDB (file name depends on phi*)."""
    pat = str(ROOT / "runs" / protein / "dewetting" / "dewetting_phi_*.pdb")
    files = [f for f in glob.glob(pat) if "phistar" not in f]
    return files[0] if files else None


def find_phistar_pdb(protein):
    f = ROOT / "runs" / protein / "dewetting" / "dewetting_phistar.pdb"
    return str(f) if f.is_file() else None


def render_map(pdb_path, out_basename, b_range, reverse=False):
    """Render a PDB colored by B-factor, three views."""
    cmd.delete("all")
    cmd.load(pdb_path, "prot")

    # Hide buried-atom sentinel and any negative B for eta map
    if "phistar" in str(pdb_path):
        cmd.hide("everything", "prot")
        cmd.show("surface", "prot and b > -7.5")
    else:
        cmd.hide("everything", "prot")
        cmd.show("surface", "prot and b > -0.5")

    # Color by B-factor
    lo, hi = b_range
    if reverse:
        cmd.spectrum("b", "blue_white_red", "prot", minimum=lo, maximum=hi)
    else:
        cmd.spectrum("b", "red_white_blue", "prot", minimum=lo, maximum=hi)

    # Common view settings
    cmd.bg_color("white")
    cmd.set("ray_shadows", 0)
    cmd.set("ray_opaque_background", 0)
    cmd.set("surface_quality", 2)
    cmd.set("ambient", 0.3)

    cmd.orient("prot")
    cmd.zoom("prot", buffer=2.0)

    base_view = cmd.get_view()

    for name, (rx, ry, rz) in VIEWS:
        cmd.set_view(base_view)
        if rx: cmd.rotate("x", rx, "prot")
        if ry: cmd.rotate("y", ry, "prot")
        if rz: cmd.rotate("z", rz, "prot")
        out = OUTDIR / f"{out_basename}_{name}.png"
        cmd.png(str(out), width=800, height=800, dpi=150, ray=1)
        print(f"  wrote {out.name}")


def main():
    for protein in PROTEINS:
        print(f"=== {protein} ===")

        eta_pdb = find_eta_pdb(protein)
        if eta_pdb:
            print(f"  eta:    {os.path.basename(eta_pdb)}")
            render_map(eta_pdb, f"{protein}_eta", ETA_RANGE, reverse=False)
        else:
            print("  eta:    (no single-window map found)")

        phistar_pdb = find_phistar_pdb(protein)
        if phistar_pdb:
            print(f"  phi*_i: dewetting_phistar.pdb")
            render_map(phistar_pdb, f"{protein}_phistar", PHISTAR_RANGE, reverse=False)
        else:
            print("  phi*_i: (no per-atom map found)")

    print(f"\nAll figures in {OUTDIR}")


main()
