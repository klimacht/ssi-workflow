"""Why is Y65 the preferred site? QM/MM single points on MD frames (PySCF, electrostatic embedding).

Per frame (written by qmmm_prep.py) and per Tyr site, the side chain from CB is the QM region
(p-cresol): cut at CA–CB, link H on the CB→CA bond at 1.09 A. The MM region is every other
atom as a force-field point charge, with CA and HA of the cut residue set to zero (they sit next
to the link H) and all charges Gaussian-smeared (--mm-radius) against over-polarisation.

States of each site, all at the MD geometry (vertical):
    PhOH     phenol, singlet             PhOH+    radical cation, doublet
    PhO-     phenolate (HH removed)      PhO.     phenoxyl radical, doublet
    IE_phenol     = E(PhOH+) - E(PhOH)   how easily the ring gives up an electron
    IE_phenolate  = E(PhO.)  - E(PhO-)   electron detachment from the phenolate (electron transfer to O3)
    E_dep         = E(PhO-)  - E(PhOH)   deprotonation; only differences between sites mean anything
                                         (proton energy and missing relaxation are the same for all)
Triad (default Y65, Y127, Y112): the three side chains as one QM region, neutral and radical
cation; the Loewdin spin population summed per ring shows where the hole sits.

Every quantity is computed twice: 'mm' (protein charges) and 'gas' (same QM coordinates, no
charges). mm - gas is the effect of the protein environment, separated from geometry.
Each radical starts from the orbitals of its closed-shell partner (same atoms) with one electron
taken out of the HOMO, instead of an atomic guess.

    python qmmm_tyr.py check   --qmmm runs/dopa/qmmm/wt                      no PySCF needed: QM regions, MM charges
    python qmmm_tyr.py run     --qmmm runs/dopa/qmmm/wt --index 7            one frame: all sites + triad
    python qmmm_tyr.py run     --qmmm runs/dopa/qmmm/wt --index 7 --part triad   only the triad (own array task)
    python qmmm_tyr.py summary --qmmm runs/dopa/qmmm/wt                      mean +- SEM per site, triad spin per ring

Results: QMMM/results/<frame>_sites.json and <frame>_triad.json (one file per part, so both
parts of a frame can run at the same time; restartable: finished states are skipped),
QMMM/summary_sites.tsv and QMMM/summary_triad.tsv.
"""
import argparse
import json
import time
from pathlib import Path

import numpy as np

HARTREE_EV = 27.211386245988
BACKBONE = {"N", "H", "CA", "HA", "C", "O", "H1", "H2", "H3", "OC1", "OC2"}
STATES = {"PhOH": (0, 0, False), "PhOH+": (1, 1, False), "PhO-": (-1, 0, True), "PhO.": (0, 1, True)}


# --------------------------------------------------------------------------- QM/MM regions

def load(path):
    f = np.load(path)
    return {k: f[k] for k in f.files}


def side_chain(fr, resid):
    """Indices of the side-chain atoms of one residue, and of its CA and HA."""
    res = np.flatnonzero(fr["resids"] == resid)
    names = fr["names"][res]
    side = res[[n not in BACKBONE for n in names]]
    ca = res[names == "CA"]
    ha = res[np.isin(names, ["HA", "HA1", "HA2"])]
    cb = res[names == "CB"]
    if len(ca) != 1 or len(cb) != 1:
        raise ValueError(f"residue {resid}: no single CA/CB")
    return side, ca[0], cb[0], ha


def region(fr, resids, deprotonate=False):
    """QM atoms (element, xyz) with link H per residue, fragment label per QM atom, MM coords/charges.

    The neutral QM side chain replaces side chain, CA and HA of each cut residue; the charge
    left on the rest of that residue (N, H, C, O: -0.114 e for Tyr) is spread evenly over those
    atoms so that the MM part of the residue keeps the residue's total charge.
    """
    elem, xyz, frag, qm_idx, zero = [], [], [], [], []
    q = fr["charges"].copy()
    for r in resids:
        side, ca, cb, ha = side_chain(fr, r)
        for i in side:
            if deprotonate and fr["names"][i] == "HH":
                continue
            elem.append(str(fr["elements"][i]))
            xyz.append(fr["xyz"][i])
            frag.append(r)
        v = fr["xyz"][ca] - fr["xyz"][cb]
        elem.append("H")
        xyz.append(fr["xyz"][cb] + 1.09 * v / np.linalg.norm(v))
        frag.append(r)
        qm_idx.extend(side)
        zero.extend([ca, *ha])
        res = np.flatnonzero(fr["resids"] == r)
        rest = np.setdiff1d(res, np.concatenate([side, [ca], ha]))
        q[rest] += (fr["charges"][res].sum() - fr["charges"][rest].sum()) / len(rest)
    q[zero] = 0.0
    mm = np.setdiff1d(np.arange(len(fr["xyz"])), qm_idx + zero)   # HH stays out of MM when removed
    return elem, np.array(xyz), np.array(frag), fr["xyz"][mm], q[mm]


# --------------------------------------------------------------------------- PySCF

def ionised_guess(mf):
    """Start density for the radical of the same atoms: the closed-shell orbitals of mf with one
    electron taken out of the HOMO (beta spin)."""
    c = mf.mo_coeff[:, mf.mo_occ > 0]
    return np.array([c @ c.T, c[:, :-1] @ c[:, :-1].T])


def scf_energy(elem, xyz, charge, spin, mm, a, frag=None, dm0=None):
    """Energy (and spin per fragment) of one state; also returns the converged SCF object."""
    from pyscf import gto, dft, qmmm, lib
    lib.num_threads(a.threads)
    t0 = time.time()
    mol = gto.M(atom=list(zip(elem, xyz.tolist())), basis=a.basis, charge=charge, spin=spin,
                unit="Angstrom", verbose=0, max_memory=a.memory_mb)
    mf = (dft.UKS(mol) if spin else dft.RKS(mol)).density_fit()
    mf.xc, mf.conv_tol, mf.max_cycle = a.xc, 1e-8, 150
    if mm is not None:
        coords, q = mm
        mf = qmmm.mm_charge(mf, coords, q, radii=np.full(len(q), a.mm_radius), unit="Angstrom")
    e = mf.kernel(dm0)
    cycles = mf.cycles
    if not mf.converged:                                       # second-order fallback
        mf = mf.newton()
        e = mf.kernel(mf.make_rdm1())
        cycles += mf.cycles
    out = {"E_h": float(e), "converged": bool(mf.converged), "cycles": int(cycles),
           "guess": "ionised" if dm0 is not None else "atomic", "time_s": round(time.time() - t0, 1)}
    if spin and frag is not None:
        da, db = mf.make_rdm1()
        s = mol.intor_symmetric("int1e_ovlp")
        w, v = np.linalg.eigh(s)
        sh = v @ np.diag(np.sqrt(w)) @ v.T
        pop = np.diag(sh @ (da - db) @ sh)
        per_atom = np.array([pop[p0:p1].sum() for _, _, p0, p1 in mol.aoslice_by_atom()])
        out["spin"] = {str(r): float(per_atom[frag == r].sum()) for r in np.unique(frag)}
    return out, mf


def run(a):
    qd = Path(a.qmmm)
    frames = [l.split("\t") for l in (qd / "frames.tsv").read_text().splitlines()[1:]]
    _, fname, system, rep, t_ps = frames[a.index][:5]
    fr = load(qd / "frames" / fname)
    res_dir = qd / "results"
    res_dir.mkdir(exist_ok=True)
    parts = [p for p in (["sites", "triad"] if a.part == "all" else [a.part]) if p == "sites" or a.triad]

    for part in parts:
        rpath = res_dir / fname.replace(".npz", f"_{part}.json")
        res = json.loads(rpath.read_text()) if rpath.exists() else {
            "frame": fname, "system": system, "rep": int(rep), "time_ps": float(t_ps),
            "settings": {"xc": a.xc, "basis": a.basis, "mm_radius": a.mm_radius}, part: {}}

        def save():
            rpath.write_text(json.dumps(res, indent=1))

        def state(e, label, st, elem, xyz, charge, spin, mm, frag, closed):
            """Compute one state unless it is in e already; returns the SCF object for a closed shell."""
            if st in e:
                return None
            dm0 = ionised_guess(closed) if spin and closed is not None else None
            e[st], mf = scf_energy(elem, xyz, charge, spin, mm, a, frag=frag, dm0=dm0)
            print(f"{fname} {label} {st:7s} {e[st]['E_h']:.6f} Eh  {e[st]['cycles']:3d} cycles "
                  f"({e[st]['guess']})  {e[st]['time_s']:5.0f} s", flush=True)
            save()
            return None if spin else mf

        for site in a.sites if part == "sites" else []:
            d = res["sites"].setdefault(str(site), {})
            for env in ("mm", "gas"):
                e = d.setdefault(env, {})
                closed = None                                  # PhOH -> PhOH+, PhO- -> PhO.
                for st, (charge, spin, deprot) in STATES.items():
                    elem, xyz, _, mmc, mmq = region(fr, [site], deprotonate=deprot)
                    closed = state(e, f"Y{site} {env:3s}", st, elem, xyz, charge, spin,
                                   (mmc, mmq) if env == "mm" else None, None, closed)
                if all(s in e for s in STATES):
                    h = lambda s: e[s]["E_h"] * HARTREE_EV
                    e["IE_phenol_eV"] = h("PhOH+") - h("PhOH")
                    e["IE_phenolate_eV"] = h("PhO.") - h("PhO-")
                    e["E_dep_eV"] = h("PhO-") - h("PhOH")
                    save()

        if part == "triad":
            elem, xyz, frag, mmc, mmq = region(fr, a.triad)
            for env in ("mm", "gas"):
                e = res["triad"].setdefault(env, {})
                closed = None
                for st, (charge, spin) in {"neutral": (0, 0), "cation": (1, 1)}.items():
                    closed = state(e, f"triad {env:3s}", st, elem, xyz, charge, spin,
                                   (mmc, mmq) if env == "mm" else None, frag, closed)
                if "neutral" in e and "cation" in e:
                    e["IE_eV"] = (e["cation"]["E_h"] - e["neutral"]["E_h"]) * HARTREE_EV
                    save()
    print(f"done {fname} ({', '.join(parts)})")


# --------------------------------------------------------------------------- check / summary

def check(a):
    qd = Path(a.qmmm)
    frames = [l.split("\t") for l in (qd / "frames.tsv").read_text().splitlines()[1:]]
    print(f"{len(frames)} frames in {qd / 'frames.tsv'}")
    fr = load(qd / "frames" / frames[0][1])
    print(f"first frame {frames[0][1]}: {len(fr['xyz'])} atoms, total charge {fr['charges'].sum():+.3f}")
    for resids, label in [([s], f"Y{s}") for s in a.sites] + ([(a.triad, "triad")] if a.triad else []):
        for deprot in (False, True) if len(resids) == 1 else (False,):
            elem, xyz, frag, mmc, mmq = region(fr, resids, deprotonate=deprot)
            dmin = np.linalg.norm(mmc[:, None] - xyz[None], axis=-1).min()
            near = mmq[np.linalg.norm(mmc[:, None] - xyz[None], axis=-1).min(1) < 6.0].sum()
            formula = "".join(f"{x}{elem.count(x)}" for x in sorted(set(elem)))
            print(f"  {label:6s} {'PhO-' if deprot else 'PhOH':5s} QM {formula:10s} ({len(elem)} atoms incl. "
                  f"{len(resids)} link H)  MM {len(mmq)} charges, sum {mmq.sum():+.3f}, within 6 A {near:+.3f}, "
                  f"closest {dmin:.2f} A")
    side, ca, cb, _ = side_chain(fr, a.sites[0])
    link = region(fr, [a.sites[0]])[1][-1]
    print(f"  link H of Y{a.sites[0]}: CB-H {np.linalg.norm(link - fr['xyz'][cb]):.2f} A, "
          f"on CB->CA: {np.degrees(np.arccos(np.dot(link - fr['xyz'][cb], fr['xyz'][ca] - fr['xyz'][cb]) / np.linalg.norm(link - fr['xyz'][cb]) / np.linalg.norm(fr['xyz'][ca] - fr['xyz'][cb]))):.1f} deg off")


def summary(a):
    qd = Path(a.qmmm)
    merged = {}                                                  # sites and triad files of one frame
    for p in sorted((qd / "results").glob("*.json")):
        r = json.loads(p.read_text())
        m = merged.setdefault(r["frame"], {"sites": {}, "triad": {}})
        m["sites"].update(r.get("sites", {}))
        m["triad"].update(r.get("triad", {}))
    res = list(merged.values())
    keys =["IE_phenol_eV", "IE_phenolate_eV", "E_dep_eV"]
    lines = ["site\tenv\tquantity\tmean\tsem\tn\trelative_to_site_mean"]
    table = {}
    for site in a.sites:
        for env in ("mm", "gas"):
            for k in keys:
                v = np.array([r["sites"][str(site)][env][k] for r in res
                              if k in r["sites"].get(str(site), {}).get(env, {})])
                if len(v):
                    table[(site, env, k)] = (v.mean(), v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else 0.0, len(v))
        for k in keys:                                           # protein effect = mm - gas, frame by frame
            v = np.array([r["sites"][str(site)]["mm"][k] - r["sites"][str(site)]["gas"][k] for r in res
                          if k in r["sites"].get(str(site), {}).get("mm", {}) and k in r["sites"][str(site)].get("gas", {})])
            if len(v):
                table[(site, "mm-gas", k)] = (v.mean(), v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else 0.0, len(v))
    for env in ("mm", "gas", "mm-gas"):
        for k in keys:
            vals = [table[(s, env, k)][0] for s in a.sites if (s, env, k) in table]
            ref = np.mean(vals) if vals else 0.0
            for s in a.sites:
                if (s, env, k) in table:
                    m, e, n = table[(s, env, k)]
                    lines.append(f"Y{s}\t{env}\t{k}\t{m:.3f}\t{e:.3f}\t{n}\t{m - ref:+.3f}")
    (qd / "summary_sites.tsv").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))

    if a.triad:
        tl = ["env\tquantity\tmean\tsem\tn"]
        for env in ("mm", "gas"):
            rs = [r["triad"][env] for r in res if "IE_eV" in r.get("triad", {}).get(env, {})]
            if not rs:
                continue
            ie = np.array([x["IE_eV"] for x in rs])
            tl.append(f"{env}\tIE_eV\t{ie.mean():.3f}\t{ie.std(ddof=1) / np.sqrt(len(ie)) if len(ie) > 1 else 0:.3f}\t{len(ie)}")
            spin = np.array([[x["cation"]["spin"][str(t)] for t in a.triad] for x in rs])
            for j, t in enumerate(a.triad):
                sv = spin[:, j]
                top = (spin.argmax(1) == j).mean()
                tl.append(f"{env}\tspin_Y{t}\t{sv.mean():.3f}\t{sv.std(ddof=1) / np.sqrt(len(sv)) if len(sv) > 1 else 0:.3f}\t{len(sv)}")
                tl.append(f"{env}\tfraction_hole_mainly_on_Y{t}\t{top:.3f}\t\t{len(sv)}")
        (qd / "summary_triad.tsv").write_text("\n".join(tl) + "\n", encoding="utf-8")
        print("\n".join(tl))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mode", choices=["check", "run", "summary"])
    ap.add_argument("--qmmm", required=True, help="directory written by qmmm_prep.py")
    ap.add_argument("--index", type=int, default=0, help="run: line of frames.tsv (array task id)")
    ap.add_argument("--part", choices=["all", "sites", "triad"], default="all",
                    help="run: single sites, the triad, or both")
    ap.add_argument("--sites", type=int, nargs="+", default=[65, 108, 112, 124, 127])
    ap.add_argument("--triad", type=int, nargs="*", default=[65, 127, 112], help="empty: no triad")
    ap.add_argument("--xc", default="wb97x")
    ap.add_argument("--basis", default="6-31+g*")
    ap.add_argument("--mm-radius", type=float, default=0.8, help="Gaussian width of the MM charges (A)")
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--memory-mb", type=int, default=8000)
    a = ap.parse_args()
    {"check": check, "run": run, "summary": summary}[a.mode](a)


if __name__ == "__main__":
    main()
