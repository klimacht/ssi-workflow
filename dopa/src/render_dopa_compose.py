#!/usr/bin/env python3
"""Final figures and movies from the raw PyMOL images of render_dopa_pymol.py.

  python3 render_dopa_compose.py --viz ANALYSIS/viz --raw ANALYSIS/viz/raw --out ANALYSIS/figures

Figures (PNG)
  exposure_<view>.png       Table 1 | MD mean | fraction of time accessible, Tyr cluster of the WT
  rmsf_wt.png               WT Calpha RMSF, three views
  drmsf.png                 RMSF(system) - RMSF(wt) for every modified system
  cluster.png               hydroxylated systems over the WT (representative structures)
  dityrosine_<system>.png   cross-linked systems over the WT
Movies (MP4, H.264, 1280 x 720)
  movie_<name>.mp4          rendering + time series with a cursor (see viz_prep.py)

ffmpeg: $FFMPEG, else ffmpeg on the PATH (Viper: module load gcc/15 ffmpeg/7.1), else imageio-ffmpeg.
"""
import argparse
import json
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.image as mpimg  # noqa: E402
from matplotlib.cm import ScalarMappable  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, Normalize  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

import viz_style as S  # noqa: E402

plt.rcParams.update({"font.size": 9, "text.color": S.INK, "axes.labelcolor": S.INK,
                     "xtick.color": S.MUTED, "ytick.color": S.MUTED, "axes.edgecolor": S.MUTED})


def cmap(name, stops):
    return LinearSegmentedColormap.from_list(name, stops)


ASA, RMSF, DRMSF = cmap("asa", S.ASA_STOPS), cmap("rmsf", S.RMSF_STOPS), cmap("drmsf", S.DRMSF_STOPS)


def colorbar(fig, rect, cm, rng, label, extend="neither"):
    ax = fig.add_axes(rect)
    cb = fig.colorbar(ScalarMappable(norm=Normalize(*rng), cmap=cm), cax=ax, orientation="horizontal",
                      extend=extend)
    cb.set_label(label, fontsize=9, labelpad=4)
    cb.ax.xaxis.set_label_position("top")
    cb.ax.tick_params(labelsize=8)
    cb.outline.set_visible(False)
    return cb


def image(fig, rect, path, title=None):
    ax = fig.add_axes(rect)
    ax.imshow(mpimg.imread(str(path)))
    ax.axis("off")
    if title:
        ax.set_title(title, fontsize=9.5)
    return ax


def system_label(system, meta):
    if system == "wt":
        return "WT"
    mods = meta["modified"]
    if all(n == "DTY" for _, n in mods):
        return "dityrosine " + "–".join(f"Y{r}" for r, _ in mods)
    (r, n), = mods
    side = system.rsplit("_", 1)[-1].upper()
    return f"{S.MODIFIED[n]}{r}, 3-OH on {side}"


def footnote(fig, text, y=0.015):
    fig.text(0.5, y, text, ha="center", fontsize=7.5, color=S.MUTED, style="italic")


def fig_exposure(raw, out, info):
    for vname, _ in S.VIEWS:
        paths = [raw / f"exposure_{t}_{vname}.png" for t in ("paper", "md", "acc")]
        if not all(p.exists() for p in paths):
            continue
        fig = plt.figure(figsize=(9.0, 4.1), dpi=200)
        titles = ("Table 1 (static, relaxed AF3 model)", f"MD mean, WT (t ≥ {info['skip_ns']:g} ns)",
                  "MD fraction of time accessible")
        for i, (p, t) in enumerate(zip(paths, titles)):
            image(fig, [i / 3, 0.22, 1 / 3, 0.71], p, t)
        colorbar(fig, [0.06, 0.1, 0.54, 0.03], ASA, S.ASA_RANGE, "atom ASA (Å²)", extend="max")
        colorbar(fig, [0.72, 0.1, 0.23, 0.03], ASA, S.ACC_RANGE, "fraction of frames with ASA > 1 Å²")
        footnote(fig, "spheres: phenol O, ortho C (CE1, CE2), meta C (CD1, CD2); "
                      "representative WT structure; heavy-atom ASA, probe 1.4 Å")
        fig.savefig(out / f"exposure_{vname}.png")
        plt.close(fig)
        print(f"  wrote exposure_{vname}.png")


def fig_rmsf(raw, out):
    paths = [raw / f"rmsf_wt_{v}.png" for v, _ in S.VIEWS]
    if not all(p.exists() for p in paths):
        return
    fig = plt.figure(figsize=(9.0, 3.0), dpi=200)
    for i, (p, (v, _)) in enumerate(zip(paths, S.VIEWS)):
        image(fig, [i / 3, 0.3, 1 / 3, 0.62], p, v)
    colorbar(fig, [0.3, 0.14, 0.4, 0.04], RMSF, S.RMSF_RANGE, "Cα RMSF (Å)", extend="both")
    footnote(fig, "WT, cartoon width and colour = Cα RMSF (mean over replicas); grey sticks: Y65, Y108, Y112, Y124, Y127")
    fig.savefig(out / "rmsf_wt.png")
    plt.close(fig)
    print("  wrote rmsf_wt.png")


def grid(n, ncol):
    nrow = int(np.ceil(n / ncol))
    return nrow, ncol


def fig_panels(raw, out, info, prefix, systems, fname, title_fn, legend=None, cbar=None, note="", ncol=4):
    paths = [(s, raw / f"{prefix}_{s}_front.png") for s in systems]
    paths = [(s, p) for s, p in paths if p.exists()]
    if not paths:
        return
    ncol = min(ncol, len(paths))
    nrow, ncol = grid(len(paths), ncol)
    h, w_px = mpimg.imread(str(paths[0][1])).shape[:2]
    w = 3.2 if ncol <= 3 else 2.4
    h_img, h_bottom = w * h / w_px + 0.3, 0.8
    fig = plt.figure(figsize=(w * ncol, h_img * nrow + h_bottom), dpi=200)
    H = h_img * nrow + h_bottom
    for k, (s, p) in enumerate(paths):
        r, c = divmod(k, ncol)
        image(fig, [c / ncol, (h_bottom + (nrow - 1 - r) * h_img) / H, 1 / ncol, (h_img - 0.25) / H], p,
              title_fn(s))
    if cbar:
        cm, rng, label, extend = cbar
        colorbar(fig, [0.3, 0.4 * h_bottom / H, 0.4, 0.1 / H], cm, rng, label, extend)
    if legend:
        fig.legend(handles=legend, loc="lower center", bbox_to_anchor=(0.5, 0.42 * h_bottom / H),
                   ncol=len(legend), frameon=False, fontsize=8.5)
    footnote(fig, note, y=0.08 * h_bottom / H)
    fig.savefig(out / fname)
    plt.close(fig)
    print(f"  wrote {fname}")


OVERLAY_LEGEND = [Line2D([], [], color=S.WT_CARBON, lw=3, alpha=0.6, label="WT"),
                  Line2D([], [], color=S.MOD_CARBON, lw=3, label="modified system"),
                  Line2D([], [], color="#4d4d4d", lw=1.5, ls=(0, (2, 2)), label="polar contact of a modified-residue O")]


def fig_dityrosine(raw, out, info):
    for system, meta in info["systems"].items():
        paths = [raw / f"dityrosine_{system}_{v}.png" for v, _ in S.VIEWS[:2]]
        if not all(p.exists() for p in paths):
            continue
        fig = plt.figure(figsize=(6.0, 4.2), dpi=200)
        for i, (p, (v, _)) in enumerate(zip(paths, S.VIEWS[:2])):
            image(fig, [i / 2, 0.15, 1 / 2, 0.7], p, v)
        fig.suptitle(system_label(system, meta) + " over the WT", fontsize=10.5, y=0.975)
        fig.legend(handles=OVERLAY_LEGEND[:2], loc="lower center", bbox_to_anchor=(0.5, 0.07), ncol=2,
                   frameon=False, fontsize=8.5)
        footnote(fig, "representative structures, fitted on the core Cα (37–140)")
        fig.savefig(out / f"dityrosine_{system}.png")
        plt.close(fig)
        print(f"  wrote dityrosine_{system}.png")


# ---------------------------------------------------------------------------
# movies
# ---------------------------------------------------------------------------

def read_rows(path):
    lines = [l.rstrip("\n").split("\t") for l in open(path) if l.strip()]
    return lines[0], lines[1:]


def ffmpeg_exe():
    exe = os.environ.get("FFMPEG") or shutil.which("ffmpeg")
    if exe:
        return exe
    try:
        import imageio_ffmpeg
        return imageio_ffmpeg.get_ffmpeg_exe()
    except ImportError:
        return None


def trace_axes(fig, rect, t, y, ylabel, cutoff=None):
    ax = fig.add_axes(rect)
    ax.plot(t, y, color=S.TRACE_FUTURE, lw=1.2)
    past, = ax.plot([], [], color=S.TRACE, lw=1.6)
    cursor = ax.axvline(t[0], color=S.INK, lw=0.8)
    dot, = ax.plot([], [], "o", color=S.TRACE, ms=6, mec="white", mew=1.0)
    if cutoff is not None:
        ax.axhline(cutoff, color=S.MUTED, lw=0.8, ls=(0, (3, 3)))
        ax.text(t[-1], cutoff, f"  {cutoff:g} Å", va="center", ha="left", fontsize=8, color=S.MUTED)
    ax.set_xlim(t[0], t[-1])
    ax.set_ylabel(ylabel)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.grid(axis="y", color="#e6e6e6", lw=0.6)
    return lambda k: (past.set_data(t[:k + 1], y[:k + 1]), cursor.set_xdata([t[k], t[k]]),
                      dot.set_data([t[k]], [y[k]]))


def movie(raw, out, viz, name, m):
    fdir = raw / f"movie_{name}"
    frames = sorted(fdir.glob("f_*.png"))
    if not frames:
        return
    head, rows = read_rows(viz / f"movie_{name}.tsv")
    col = {h: i for i, h in enumerate(head)}
    idx = [int(f.stem[2:]) for f in frames]
    t = np.array([float(r[0]) for r in rows])
    a, b = m["pair"]
    fig = plt.figure(figsize=(12.8, 7.2), dpi=100)
    img = fig.add_axes([0.0, 0.0, 0.5625, 1.0])
    im = img.imshow(mpimg.imread(str(frames[0])))
    img.axis("off")
    x0, w = 0.63, 0.33
    if m["kind"] == "exposure":
        tyr = m["exposure_tyr"]
        title = f"Y{tyr} opens transiently (WT)"
        y1 = np.array([float(r[col[f"Y{tyr}:O"]]) for r in rows])
        y2 = np.array([float(r[col[f"Y{tyr}:CE1"]]) + float(r[col[f"Y{tyr}:CE2"]]) for r in rows])
        updates = [trace_axes(fig, [x0, 0.56, w, 0.22], t, y1, f"Y{tyr} phenol O\nASA (Å²)"),
                   trace_axes(fig, [x0, 0.27, w, 0.22], t, y2, f"Y{tyr} ortho C\n(CE1+CE2) ASA (Å²)")]
        note = "spheres: reactive Tyr atoms, coloured by their ASA in this frame"
    else:
        title = f"Y{a}–Y{b} stay close (WT)"
        y1 = np.array([float(r[col["pair_d"]]) for r in rows])
        updates = [trace_axes(fig, [x0, 0.30, w, 0.45], t, y1, f"shortest ortho C···C\nY{a}–Y{b} (Å)",
                              cutoff=m["cutoff_A"])]
        note = f"dashed line: closest ortho carbons, orange when ≤ {m['cutoff_A']:g} Å"
    fig.axes[-1].set_xlabel("time (ns)")
    fig.text(x0, 0.92, title, fontsize=15, weight="bold")
    fig.text(x0, 0.875, f"replica {m['rep']}, coordinates averaged over {m['smooth_ns']:g} ns", fontsize=9.5,
             color=S.MUTED)
    tlabel = fig.text(x0, 0.81, "", fontsize=13)
    colorbar(fig, [x0, 0.105, w, 0.022], ASA, S.ASA_RANGE, "atom ASA (Å²)", extend="max")
    fig.text(x0, 0.03, note, fontsize=8.5, color=S.MUTED)
    cdir = raw / f"movie_{name}_composed"
    cdir.mkdir(exist_ok=True)
    for old in cdir.glob("c_*.png"):
        old.unlink()
    for n, (k, f) in enumerate(zip(idx, frames)):
        im.set_data(mpimg.imread(str(f)))
        for u in updates:
            u(k)
        tlabel.set_text(f"t = {t[k]:6.1f} ns")
        fig.savefig(cdir / f"c_{n:04d}.png", dpi=100)
    plt.close(fig)
    exe = ffmpeg_exe()
    if not exe:
        print(f"  movie {name}: {len(frames)} composed frames in {cdir}, no ffmpeg for the MP4")
        return
    mp4 = out / f"movie_{name}.mp4"
    subprocess.run([exe, "-y", "-loglevel", "error", "-framerate", str(m["fps"]), "-i", str(cdir / "c_%04d.png"),
                    "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "20", "-movflags", "+faststart",
                    str(mp4)], check=True)
    print(f"  wrote {mp4.name} ({len(frames)} frames, {len(frames) / m['fps']:.1f} s)")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--viz", required=True, help="output directory of viz_prep.py")
    ap.add_argument("--raw", required=True, help="output directory of render_dopa_pymol.py")
    ap.add_argument("--out", required=True)
    ap.add_argument("--no-movies", action="store_true")
    args = ap.parse_args()

    viz, raw, out = Path(args.viz), Path(args.raw), Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    info = json.loads((viz / "viz_info.json").read_text())
    systems = info["systems"]
    fig_exposure(raw, out, info)
    fig_rmsf(raw, out)
    mod = [s for s in systems if s != "wt"]
    fig_panels(raw, out, info, "drmsf", mod, "drmsf.png", lambda s: system_label(s, systems[s]),
               cbar=(DRMSF, S.DRMSF_RANGE, "RMSF − RMSF(WT) (Å)", "both"), ncol=3,
               note="representative structure of each system; orange sticks: modified residues")
    hyd = [s for s in mod if any(n in ("DAH", "TOP") for _, n in systems[s]["modified"])]
    fig_panels(raw, out, info, "cluster", hyd, "cluster.png", lambda s: system_label(s, systems[s]),
               legend=OVERLAY_LEGEND, note="representative structures, fitted on the core Cα (37–140)")
    fig_dityrosine(raw, out, info)
    if not args.no_movies:
        for name, m in info["movies"].items():
            movie(raw, out, viz, name, m)
    print(f"figures in {out}")


if __name__ == "__main__":
    main()
