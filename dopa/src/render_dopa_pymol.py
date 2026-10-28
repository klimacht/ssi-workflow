#!/usr/bin/env python3
"""Raw PyMOL images and movie frames from the output of viz_prep.py.

  python3 render_dopa_pymol.py --viz ANALYSIS/viz --out ANALYSIS/viz/raw [--scenes ...]

Run with the Python of a PyMOL installation (Viper: /u/tklimach/conda-envs/pymol_env, see
dopa/slurm/16_render.slurm). Writes PNGs with a transparent background; render_dopa_compose.py
then adds colorbars, titles and the movie side panels (as render_maps_pymol.py /
render_maps_compose.py do for the SSI maps).

Scenes
  exposure    reactive Tyr atoms (O, CE1, CE2, CD1, CD2) as spheres coloured by ASA: Table 1
              (exposure_paper.pdb), MD mean and fraction of time accessible (exposure_md.pdb);
              Tyr cluster in three views
  flex        WT putty cartoon, width and colour = Calpha RMSF (three views); per modified system
              the cartoon coloured by RMSF - RMSF(wt), modified residues as sticks
  cluster     per hydroxylated system: Tyr cluster over the WT (grey), polar contacts of the
              oxygens of the modified residue
  dityrosine  cross-linked systems over the WT
  movies      frames of every movie in viz_info.json, reactive atoms coloured by their ASA in
              that frame; for the pair movie the closest ortho C...C distance as a dashed line
"""
import argparse
import json
import os
from pathlib import Path

import pymol
pymol.finish_launching(["pymol", "-qc"])
from pymol import cmd, util  # noqa: E402

import viz_style as S  # noqa: E402

CL = "65+108+112+124+127"
REACT = "+".join(S.REACTIVE)
SIDE = "(sidechain or name CA) and not hydro"


def palette(prefix, stops):
    names = []
    for i, h in enumerate(stops):
        cmd.set_color(f"{prefix}{i}", S.hex_rgb(h))
        names.append(f"{prefix}{i}")
    return " ".join(names)


def setup():
    cmd.reinitialize()
    cmd.feedback("disable", "all", "actions details")
    cmd.bg_color("white")
    cmd.set("max_threads", len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else os.cpu_count())
    for key, value in (("ray_opaque_background", 0), ("ray_shadows", 0), ("ambient", 0.35),
                       ("specular", 0.15), ("antialias", 2), ("stick_radius", 0.17),
                       ("sphere_scale", 0.42), ("cartoon_transparency", 0.55),
                       ("label_size", 17), ("label_font_id", 7), ("label_position", (0.0, 1.2, 10.0)),     # larger z ends up in front of the clipping plane
                       ("float_labels", 1), ("dash_gap", 0.25), ("dash_radius", 0.07),
                       ("dash_color", "grey30"), ("depth_cue", 0), ("orthoscopic", 1),
                       ("two_sided_lighting", 1)):
        cmd.set(key, value)
    cmd.set_color("wtc", S.hex_rgb(S.WT_CARBON))
    cmd.set_color("modc", S.hex_rgb(S.MOD_CARBON))
    cmd.set_color("contact", S.hex_rgb(S.CONTACT))
    cmd.set_color("nocontact", S.hex_rgb(S.NO_CONTACT))
    cmd.set_color("ink", S.hex_rgb(S.INK))
    cmd.set("label_color", "ink")
    return {"asa": palette("asa", S.ASA_STOPS), "rmsf": palette("rmsf", S.RMSF_STOPS),
            "drmsf": palette("drm", S.DRMSF_STOPS)}


def label_cluster(obj, resnames=None):
    resnames = resnames or {}
    for resi in CL.split("+"):
        text = S.LABEL[resnames.get(int(resi), "TYR")] + resi
        cmd.label(f"{obj} and resi {resi} and name CA", repr(text))


def cluster_view(obj, buffer=3.0):
    sel = f"{obj} and resi {CL} and {SIDE}"
    cmd.orient(sel)
    cmd.zoom(sel, buffer)
    return cmd.get_view()


def closer(view, factor):
    """Move the camera towards the origin of rotation (zoom in), keeping the clipping planes in place."""
    v = list(view)
    shift = abs(v[11]) * (1.0 - factor)
    v[11] *= factor
    v[15], v[16] = max(1.0, v[15] - shift), v[16] - shift
    return v


def png(path, view, angle=0, px=S.STILL_PX, aspect=1.0):
    cmd.set_view(view)
    if angle:
        cmd.turn("y", angle)
    cmd.png(str(path), width=px, height=int(px * aspect) // 2 * 2, dpi=150, ray=1)
    print(f"  wrote {Path(path).name}")


def base_style(obj, cartoon="grey90"):
    cmd.hide("everything", obj)
    cmd.show("cartoon", obj)
    cmd.set("cartoon_color", cartoon, obj)
    cmd.show("sticks", f"{obj} and resi {CL} and {SIDE}")
    cmd.color("wtc", f"{obj} and resi {CL} and {SIDE}")


def scene_exposure(viz, out, px):
    pal = setup()
    for obj in ("paper", "md"):
        cmd.load(str(viz / f"exposure_{obj}.pdb"), obj)
        base_style(obj)
        cmd.show("spheres", f"{obj} and resi {CL} and name {REACT}")
        label_cluster(obj)
    view = cluster_view("md")
    for tag, obj, prop, rng in (("paper", "paper", "b", S.ASA_RANGE), ("md", "md", "b", S.ASA_RANGE),
                                ("acc", "md", "q", S.ACC_RANGE)):
        cmd.disable("all")
        cmd.enable(obj)
        cmd.spectrum(prop, pal["asa"], f"{obj} and resi {CL} and name {REACT}", minimum=rng[0], maximum=rng[1])
        for vname, angle in S.VIEWS:
            png(out / f"exposure_{tag}_{vname}.png", view, angle, px)


def scene_flex(viz, info, out, px):
    pal = setup()
    cmd.load(str(viz / "rep_wt.pdb"), "wt")
    cmd.hide("everything")
    cmd.show("cartoon", "wt")
    cmd.cartoon("putty", "wt")
    cmd.set("cartoon_transparency", 0.0)
    cmd.spectrum("b", pal["rmsf"], "wt", minimum=S.RMSF_RANGE[0], maximum=S.RMSF_RANGE[1])
    cmd.show("sticks", f"wt and resi {CL} and {SIDE}")
    cmd.color("wtc", f"wt and resi {CL} and {SIDE}")
    cmd.viewport(px, int(px * S.DOMAIN_ASPECT))          # zoom for the wide frame
    cmd.orient("wt and name CA")
    cmd.zoom("wt", 2.0, complete=1)
    view = closer(cmd.get_view(), 0.62)                  # zoom fits the bounding sphere, too loose here
    for vname, angle in S.VIEWS:
        png(out / f"rmsf_wt_{vname}.png", view, angle, px, S.DOMAIN_ASPECT)
    for system, meta in info["systems"].items():
        if system == "wt" or not meta.get("rmsf"):
            continue
        pal = setup()
        cmd.load(str(viz / f"rep_{system}.pdb"), "s")
        cmd.hide("everything")
        cmd.show("cartoon", "s")
        cmd.set("cartoon_transparency", 0.0)
        cmd.spectrum("q", pal["drmsf"], "s", minimum=S.DRMSF_RANGE[0], maximum=S.DRMSF_RANGE[1])
        mod = "+".join(str(r) for r, _ in meta["modified"])
        cmd.show("sticks", f"s and resi {mod} and {SIDE}")
        cmd.color("modc", f"s and resi {mod} and sidechain and elem C")      # CA keeps the cartoon colour
        util.cnc(f"s and resi {mod} and sidechain")
        cmd.viewport(px, int(px * S.DOMAIN_ASPECT))
        png(out / f"drmsf_{system}_front.png", view, 0, px, S.DOMAIN_ASPECT)


def overlay(viz, system, meta):
    setup()
    cmd.load(str(viz / "rep_wt.pdb"), "wt")
    cmd.load(str(viz / f"rep_{system}.pdb"), "mod")
    cmd.hide("everything")
    cmd.show("cartoon", "mod")
    cmd.set("cartoon_color", "grey90", "mod")
    cmd.show("sticks", f"wt and resi {CL} and {SIDE}")
    cmd.color("wtc", "wt")
    cmd.set("stick_transparency", 0.45, "wt")
    mod_res = "+".join(str(r) for r, _ in meta["modified"])
    cmd.show("sticks", f"mod and resi {CL} and {SIDE}")
    cmd.show("sticks", f"mod and resi {mod_res} and hydro and neighbor elem O")
    cmd.color("modc", f"mod and elem C")
    util.cnc("mod")
    label_cluster("mod", {r: n for r, n in meta["modified"]})
    return mod_res


def scene_cluster(viz, info, out, px):
    for system, meta in info["systems"].items():
        if not any(n in ("DAH", "TOP") for _, n in meta["modified"]):
            continue
        mod_res = overlay(viz, system, meta)
        cmd.distance("hb", f"mod and resi {mod_res} and elem O and not name O",
                     f"mod and not resi {mod_res} and (elem N or elem O)", mode=2)
        cmd.hide("labels", "hb")
        view = cluster_view("wt")
        png(out / f"cluster_{system}_front.png", view, 0, px)


def scene_dityrosine(viz, info, out, px):
    for system, meta in info["systems"].items():
        if not any(n == "DTY" for _, n in meta["modified"]):
            continue
        mod_res = overlay(viz, system, meta)
        sel = f"mod and resi {mod_res} and {SIDE}"
        cmd.orient(sel)
        cmd.zoom(sel, 8.0)
        view = cmd.get_view()
        for vname, angle in S.VIEWS[:2]:
            png(out / f"dityrosine_{system}_{vname}.png", view, angle, px)


def read_rows(path):
    lines = [l.rstrip("\n").split("\t") for l in open(path) if l.strip()]
    return [dict(zip(lines[0], l)) for l in lines[1:]]


def scene_movies(viz, info, out, px, every):
    """Each frame is a single-state copy of that trajectory state ("cur"), styled like the stills."""
    for name, m in info["movies"].items():
        pal = setup()
        cmd.load(str(viz / f"movie_{name}.pdb"), "m")
        cmd.load_traj(str(viz / f"movie_{name}.xtc"), "m", state=1)
        cmd.disable("m")
        rows = read_rows(viz / f"movie_{name}.tsv")
        n = min(cmd.count_states("m"), len(rows))
        a, b = m["pair"]
        react = f"cur and resi {CL} and name {REACT}"
        fdir = out / f"movie_{name}"
        fdir.mkdir(parents=True, exist_ok=True)
        view = None
        for k in range(0, n, every):
            row = rows[k]
            vals = {}
            for t in CL.split("+"):
                for atom in S.REACTIVE:
                    key = f"Y{t}:{S.PAPER_NAME[atom]}"
                    if key in row:
                        vals[(t, atom)] = float(row[key])
            cmd.delete("cur")
            cmd.create("cur", "m", k + 1, 1)
            base_style("cur")
            cmd.show("spheres", react)
            if m["kind"] == "pair":
                cmd.set("stick_radius", 0.24, "cur")
            cmd.alter(react, "b = vals.get((resi, name), 0.0)", space={"vals": vals})
            cmd.spectrum("b", pal["asa"], react, minimum=S.ASA_RANGE[0], maximum=S.ASA_RANGE[1])
            label_cluster("cur")
            if m["kind"] == "pair":
                cmd.delete("dd")
                cmd.distance("dd", f"cur and resi {a} and name {row['pair_a']}",
                             f"cur and resi {b} and name {row['pair_b']}")
                cmd.hide("labels", "dd")
                cmd.set("dash_radius", 0.11, "dd")
                cmd.set("dash_color", "contact" if float(row["pair_d"]) <= m["cutoff_A"] else "nocontact", "dd")
            if view is None:                       # fixed camera: the cluster in the first frame
                view = cluster_view("cur")
            cmd.set_view(view)
            cmd.png(str(fdir / f"f_{k:04d}.png"), width=px, height=px, dpi=150, ray=1)
        print(f"  movie {name}: {len(range(0, n, every))} frames in {fdir}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--viz", required=True, help="output directory of viz_prep.py")
    ap.add_argument("--out", required=True)
    ap.add_argument("--scenes", nargs="+", default=["exposure", "flex", "cluster", "dityrosine", "movies"])
    ap.add_argument("--scale", type=float, default=1.0, help="image size factor (e.g. 0.5 for quick tests)")
    ap.add_argument("--movie-every", type=int, default=1, help="render every n-th movie frame")
    args = ap.parse_args()

    viz, out = Path(args.viz), Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    info = json.loads((viz / "viz_info.json").read_text())
    still = int(S.STILL_PX * args.scale)
    for scene in args.scenes:
        print(f"=== {scene}")
        if scene == "exposure":
            scene_exposure(viz, out, still)
        elif scene == "flex":
            scene_flex(viz, info, out, still)
        elif scene == "cluster":
            scene_cluster(viz, info, out, still)
        elif scene == "dityrosine":
            scene_dityrosine(viz, info, out, still)
        elif scene == "movies":
            scene_movies(viz, info, out, int(S.MOVIE_PX * args.scale) // 2 * 2, args.movie_every)
        else:
            raise SystemExit(f"unknown scene {scene}")


if __name__ == "__main__":
    main()
