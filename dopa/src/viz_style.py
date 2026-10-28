"""Colour scales, ranges and labels shared by render_dopa_pymol.py (PyMOL) and
render_dopa_compose.py (matplotlib), so that the colorbars match the renderings.

Scales are evenly spaced colour stops: PyMOL interpolates `spectrum` between them,
matplotlib builds the same LinearSegmentedColormap.
"""

# solvent exposure (A^2, fraction of time): one hue, light = buried, dark = exposed
ASA_STOPS = ["#fbe6d3", "#f5b783", "#e8803f", "#bd5120", "#7a2d0c"]
ASA_RANGE = (0.0, 30.0)
ACC_RANGE = (0.0, 1.0)

# Calpha RMSF of the wild type (A): one hue
RMSF_STOPS = ["#e3eef8", "#a7c9e6", "#5f9ccf", "#2f68a6", "#16386b"]
RMSF_RANGE = (0.5, 3.0)

# RMSF(system) - RMSF(wt) (A): blue = stiffer, red = more flexible, neutral grey middle
# (poles as in the SSI maps, scripts/render_maps_compose.py)
DRMSF_STOPS = ["#2050a0", "#8ea6cf", "#e6e6e6", "#df8f95", "#c0202c"]
DRMSF_RANGE = (-1.5, 1.5)

# structure overlays: wild type vs. modified system
WT_CARBON = "#8c8c8c"
MOD_CARBON = "#e8803f"
CONTACT = "#e8803f"          # distance at or below the coupling cutoff
NO_CONTACT = "#8c8c8c"

TRACE = "#2f68a6"            # movie side panels
TRACE_FUTURE = "#c9d6e6"
INK = "#222222"
MUTED = "#666666"

STILL_PX = 600               # raw renderings: width; cluster views are square,
DOMAIN_ASPECT = 0.6          # whole-domain views height / width
MOVIE_PX = 720               # square movie renderings (composed to 1280 x 720)
MOVIE_FPS = 25

MODIFIED = {"DAH": "DOPA", "TOP": "TOPA", "DTY": "dityrosine"}
LABEL = {"DAH": "DOPA", "TOP": "TOPA", "DTY": "Y", "TYR": "Y"}    # residue labels in the images
REACTIVE = ("OH", "CE1", "CE2", "CD1", "CD2")
PAPER_NAME = {"OH": "O", "CE1": "CE1", "CE2": "CE2", "CD1": "CD1", "CD2": "CD2"}
VIEWS = (("front", 0), ("side", 90), ("back", 180))


def hex_rgb(h):
    h = h.lstrip("#")
    return [int(h[i:i + 2], 16) / 255.0 for i in (0, 2, 4)]
