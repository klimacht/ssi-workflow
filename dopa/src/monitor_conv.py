#!/usr/bin/env python3
"""Convergence numbers and terminal plots for one DOPA system (used by dopa_monitor.sh).

Per replica:
  energies   from the GROMACS logs (md.log, else npt.log / nvt.log):
             temperature, pressure, volume (box), potential energy
  structure  from rep*/md.xtc (protein only), every --stride-th frame:
             Calpha RMSD of the helical core vs. the first production frame,
             radius of gyration (Calpha), ring-centroid distance MUT-PARTNER
Structure values are cached in rep*/.monitor/ so a refresh only reads new frames.
Plots are drawn with gnuplot's dumb terminal (replicas as 1, 2, 3, 4).
"""
import argparse
import json
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np

RING = ("CG", "CD1", "CD2", "CE1", "CE2", "CZ")


# ---------------------------------------------------------------------------
# GROMACS log energies
# ---------------------------------------------------------------------------
def parse_log(path):
    """List of dicts (time in ps + all energy terms + Volume) from a GROMACS log."""
    rows, t = [], None
    lines = Path(path).read_text(errors="replace").split("\n")
    i = 0
    while i < len(lines):
        line = lines[i]
        if "A V E R A G E S" in line:
            break
        if line.strip().startswith("Step") and "Time" in line and i + 1 < len(lines):
            f = lines[i + 1].split()
            try:
                t = float(f[1])
            except (IndexError, ValueError):
                pass
            i += 2
            continue
        if "Energies (kJ/mol)" in line and t is not None:
            vals, j = {}, i + 1
            while j + 1 < len(lines) and lines[j].strip():
                names = [lines[j][k:k + 15].strip() for k in range(0, len(lines[j]), 15)]
                nums = lines[j + 1].split()
                for nm, v in zip([n for n in names if n], nums):
                    try:
                        vals[nm] = float(v)
                    except ValueError:
                        pass
                j += 2
            if j + 2 < len(lines) and "Box-X" in lines[j + 1]:
                box = [float(v) for v in lines[j + 2].split()[:3]]
                vals["Volume"] = box[0] * box[1] * box[2]
                j += 3
            vals["time"] = t
            rows.append(vals)
            i = j
            continue
        i += 1
    return rows


# ---------------------------------------------------------------------------
# structure from the trajectory (cached, incremental)
# ---------------------------------------------------------------------------
def update_structure(rep, sysdir, args):
    """Append new frames to rep/.monitor/struct.tsv; return array (t_ns, rmsd, rg, dist)."""
    xtc = rep / "md.xtc"
    cache = rep / ".monitor"
    table = cache / "struct.tsv"
    state_f = cache / "state.json"
    ref_f = cache / "ref_core.npy"
    old = np.loadtxt(table, ndmin=2) if table.exists() and table.stat().st_size else np.empty((0, 4))
    if not xtc.exists() or xtc.stat().st_size == 0:
        return old
    state = json.loads(state_f.read_text()) if state_f.exists() else {}
    if state.get("xtc_size") == xtc.stat().st_size:
        return old

    import MDAnalysis as mda
    from MDAnalysis.analysis.align import rotation_matrix
    from MDAnalysis.lib.distances import minimize_vectors

    cache.mkdir(exist_ok=True)
    try:
        u = mda.Universe(str(sysdir / "processed.gro"), str(xtc))
    except Exception as exc:                      # file being written, offsets, ...
        print(f"      ({rep.name}: trajectory not readable yet: {exc.__class__.__name__})")
        return old
    ca = u.select_atoms(f"name CA and resid {args.first}:{args.last}")
    core_mask = np.isin(ca.resids, np.arange(args.core[0], args.core[1] + 1))
    rings = {}
    for r in (args.mut, args.partner):
        g = u.select_atoms(f"resid {r} and name " + " ".join(RING))
        if len(g) == 6:
            rings[r] = (g, ca.resids.tolist().index(r))

    def frame_values(ts):
        box = ts.dimensions
        raw = ca.positions
        steps = minimize_vectors(np.diff(raw, axis=0), box)
        pos = np.vstack([raw[:1], raw[0] + np.cumsum(steps, axis=0)])   # whole Calpha trace
        cent = {}
        for r, (g, k) in rings.items():
            cent[r] = pos[k] + minimize_vectors(g.positions - raw[k], box).mean(axis=0)
        return pos, cent

    if ref_f.exists():
        ref = np.load(ref_f)
    else:
        u.trajectory[0]
        pos, _ = frame_values(u.trajectory.ts)
        ref = pos[core_mask] - pos[core_mask].mean(axis=0)
        np.save(ref_f, ref)

    t_first = state.get("t_first")
    start = state.get("next_frame", 0)
    rows = []
    try:
        n = len(u.trajectory)
        for ts in u.trajectory[start:n:args.stride]:
            if t_first is None:
                t_first = ts.time
            pos, cent = frame_values(ts)
            core = pos[core_mask] - pos[core_mask].mean(axis=0)
            _, rmsd = rotation_matrix(core, ref)
            rg = np.sqrt(((pos - pos.mean(axis=0)) ** 2).sum(axis=1).mean())
            dist = (np.linalg.norm(cent[args.mut] - cent[args.partner])
                    if args.mut in cent and args.partner in cent else np.nan)
            rows.append(((ts.time - t_first) / 1000.0, rmsd, rg, dist))
            start = ts.frame + args.stride
    except Exception as exc:
        print(f"      ({rep.name}: stopped reading at frame {start}: {exc.__class__.__name__})")
    new = np.vstack([old, np.array(rows)]) if rows else old
    tmp = cache / "struct.tsv.tmp"
    np.savetxt(tmp, new, fmt="%.4f", header="t_ns rmsd_A rg_A dist_A")
    os.replace(tmp, table)
    state.update(next_frame=start, t_first=t_first, xtc_size=xtc.stat().st_size)
    state_f.write_text(json.dumps(state))
    return new


# ---------------------------------------------------------------------------
# output
# ---------------------------------------------------------------------------
def tail_stats(t, y, window):
    if len(t) == 0:
        return None
    m = t >= t[-1] - window
    yy = y[m][np.isfinite(y[m])]
    if len(yy) == 0:
        return None
    return yy.mean(), yy.std()


def fmt(stat, f="{:.2f}"):
    if stat is None:
        return "-"
    return (f + "±" + f).format(*stat)


def gnuplot(series, title, ylabel, width, height):
    """series: list of (label, x, y). Returns the dumb-terminal plot as text."""
    series = [(lab, x, y) for lab, x, y in series if len(x) > 1]
    if not series or not shutil.which("gnuplot"):
        return ""
    data = []
    plots = []
    for k, (lab, x, y) in enumerate(series):
        block = "\n".join(f"{a:.4f} {b:.5g}" for a, b in zip(x, y) if np.isfinite(b))
        data.append(block + "\ne")
        plots.append(f"'-' using 1:2 with points pt '{lab}' notitle")
    script = (f"set terminal dumb size {width},{height}\n"
              f"set title '{title}'\nset ylabel '{ylabel}'\nset xlabel 'ns'\n"
              "set tics nomirror\nset border 3\n"
              "plot " + ", ".join(plots) + "\n" + "\n".join(data) + "\n")
    try:
        out = subprocess.run(["gnuplot"], input=script, capture_output=True, text=True, timeout=30).stdout
    except Exception:
        return ""
    lines = [l.rstrip() for l in out.split("\n")]
    while lines and not lines[0].strip():
        lines.pop(0)
    while lines and not lines[-1].strip():
        lines.pop()
    return "\n".join("    " + l for l in lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sysdir", required=True)
    ap.add_argument("--nrep", type=int, default=4)
    ap.add_argument("--first", type=int, default=34)
    ap.add_argument("--last", type=int, default=144)
    ap.add_argument("--mut", type=int, default=65)
    ap.add_argument("--partner", type=int, default=112)
    ap.add_argument("--core", default="37:140")
    ap.add_argument("--stride", type=int, default=10, help="every n-th frame (10 ps each)")
    ap.add_argument("--window", type=float, default=10.0, help="ns for the averages")
    ap.add_argument("--plot", choices=("none", "compact", "full"), default="compact")
    ap.add_argument("--width", type=int, default=0)
    args = ap.parse_args()
    args.core = tuple(int(v) for v in args.core.split(":"))

    sysdir = Path(args.sysdir)
    cols = args.width or int(os.environ.get("COLUMNS", 0)) or shutil.get_terminal_size((110, 40)).columns
    width = max(60, min(cols - 6, 140))

    reps = [sysdir / f"rep{r}" for r in range(1, args.nrep + 1)]
    energy, struct, phase = {}, {}, {}
    for rep in reps:
        for name in ("md", "npt", "nvt"):
            log = rep / f"{name}.log"
            if log.exists() and log.stat().st_size:
                energy[rep.name] = parse_log(log)
                phase[rep.name] = name
                break
        if (rep / "md.xtc").exists():
            struct[rep.name] = update_structure(rep, sysdir, args)

    win_ps = args.window * 1000.0
    for rep in reps:
        e = energy.get(rep.name, [])
        parts = []
        if e:
            t = np.array([r["time"] for r in e])
            def col(key):
                return np.array([r.get(key, np.nan) for r in e])
            parts.append(f"{phase[rep.name].upper():3s} T {fmt(tail_stats(t, col('Temperature'), win_ps), '{:.1f}')} K")
            parts.append(f"P {fmt(tail_stats(t, col('Pressure (bar)'), win_ps), '{:.0f}')} bar")
            if "Volume" in e[-1]:
                parts.append(f"V {fmt(tail_stats(t, col('Volume'), win_ps), '{:.2f}')} nm3")
        s = struct.get(rep.name)
        if s is not None and len(s):
            parts.append(f"RMSD {fmt(tail_stats(s[:, 0], s[:, 1], args.window))} A")
            parts.append(f"Rg {fmt(tail_stats(s[:, 0], s[:, 2], args.window))} A")
            parts.append(f"d{args.mut}-{args.partner} {fmt(tail_stats(s[:, 0], s[:, 3], args.window))} A")
        if parts:
            print(f"    {rep.name:5s} " + "   ".join(parts))
    if any(energy.values()) or struct:
        print(f"          (mean±SD over the last {args.window:g} ns of each replica)")

    if args.plot == "none":
        return
    h_small, h_big = 13, 16

    def struct_series(k):
        return [(str(i + 1), struct[r.name][:, 0], struct[r.name][:, k])
                for i, r in enumerate(reps) if r.name in struct and len(struct[r.name])]

    def energy_series(key, scale=1.0):
        out = []
        for i, r in enumerate(reps):
            e = energy.get(r.name, [])
            if e:
                out.append((str(i + 1), np.array([x["time"] for x in e]) / 1000.0,
                            np.array([x.get(key, np.nan) for x in e]) * scale))
        return out

    plots = []
    if struct:
        plots.append((struct_series(1), f"Calpha RMSD core {args.core[0]}-{args.core[1]} (A), replicas 1-4", "A"))
        if args.plot == "full":
            plots.append((struct_series(3), f"ring distance {args.mut}-{args.partner} (A)", "A"))
            plots.append((struct_series(2), "radius of gyration Calpha (A)", "A"))
    if args.plot == "full" or not struct:
        if any(phase.get(r.name) in ("npt", "md") for r in reps):
            plots.append((energy_series("Volume"), f"box volume (nm3), {'/'.join(sorted(set(phase.values())))}", "nm3"))
        if args.plot == "full":
            plots.append((energy_series("Temperature"), "temperature (K)", "K"))
            plots.append((energy_series("Potential", 1e-3), "potential energy (MJ/mol)", "MJ/mol"))
    for series, title, ylab in plots:
        txt = gnuplot(series, title, ylab, width, h_big if args.plot == "full" else h_small)
        if txt:
            print(txt)


if __name__ == "__main__":
    main()
