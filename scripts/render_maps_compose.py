#!/usr/bin/env python3
"""
render_maps_compose.py — STEP 2: add matplotlib colorbar + labels to raw
PyMOL images and write the final PNGs.

Run via:
  module purge
  module load python-waterboa/2025.06
  python3 render_maps_compose.py
"""

import csv
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable

ROOT = Path("/ptmp/tklimach/ssi_workflow")
RAWDIR = Path("/ptmp/tklimach/ssi-workflow/examples/figures/_raw")
OUTDIR = Path("/ptmp/tklimach/ssi-workflow/examples/figures")

PROTEINS = ["1DPX", "1HRC", "1FS3", "1RBB", "3LDJ", "4IBA", "PHLP5"]
PHISTAR_RANGE = (-3.0, 3.0)
ETA_RANGE = (0.0, 1.0)
VIEWS = ["front", "side", "back"]

RWB = LinearSegmentedColormap.from_list(
    "rwb_custom",
    [(0.0, "#c0202c"), (0.5, "#ffffff"), (1.0, "#2050a0")])


def get_phistar_dir(protein):
    tsv = ROOT / "runs" / protein / "analysis" / "N_vs_phi.tsv"
    if not tsv.is_file(): return None
    with open(tsv) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r.get("phi_star", "0") == "1":
                return r["dirname"]
    return None


def get_center_info(protein):
    """Read centering metadata from the sidecar .center file."""
    f = ROOT / "runs" / protein / "dewetting" / "dewetting_phistar.center"
    if not f.is_file():
        return None
    meta = {}
    with open(f) as fh:
        for line in fh:
            k, v = line.strip().split("	", 1)
            meta[k] = v
    return meta


def compose(raw_png, final_png, b_range, label, subtitle):
    img = mpimg.imread(str(raw_png))
    fig = plt.figure(figsize=(5.0, 6.2), dpi=100)
    ax_img = fig.add_axes([0.0, 0.20, 1.0, 0.78])
    ax_img.imshow(img); ax_img.axis("off")
    ax_cb = fig.add_axes([0.15, 0.11, 0.70, 0.030])
    cb = fig.colorbar(ScalarMappable(norm=Normalize(*b_range), cmap=RWB),
                      cax=ax_cb, orientation="horizontal")
    cb.set_label(label, fontsize=9)
    cb.ax.tick_params(labelsize=8)
    if subtitle:
        fig.text(0.5, 0.025, subtitle, ha="center", fontsize=7,
                 color="#444", style="italic")
    fig.savefig(str(final_png), dpi=100)
    plt.close(fig)


def main():
    for protein in PROTEINS:
        print(f"=== {protein} ===")
        phistar_dir = get_phistar_dir(protein)
        phi_label = phistar_dir.replace("phi_", "").replace("p", ".") \
                    if phistar_dir else "?"
        center_info = get_center_info(protein)

        for view in VIEWS:
            # eta map
            raw = RAWDIR / f"{protein}_eta_{view}.png"
            if raw.is_file():
                out = OUTDIR / f"{protein}_eta_{view}.png"
                compose(raw, out, ETA_RANGE,
                        f"η = ⟨n_i⟩_φ* / ⟨n_i⟩₀  at φ* = {phi_label} kJ/mol",
                        "single-window dewetting map")
                print(f"  wrote {out.name}")

            # phistar map
            raw = RAWDIR / f"{protein}_phistar_{view}.png"
            if raw.is_file() and center_info is not None:
                out = OUTDIR / f"{protein}_phistar_{view}.png"
                label = center_info.get("center_label", "?")
                compose(raw, out, PHISTAR_RANGE,
                        "Δφ*_i = φ*_i − φ*_m  [kJ/mol]",
                        f"reference: {label} kJ/mol")
                print(f"  wrote {out.name}")


main()
