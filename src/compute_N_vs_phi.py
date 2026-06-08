#!/usr/bin/env python3
"""
Compute <N>_phi and susceptibility from COLVAR files.
Produces the N vs phi curve and identifies phi*.

Usage:
    python3 compute_N_vs_phi.py --pdbid 1DPX
                                [--root __SSI_ROOT__]
                                [--t-start 1000] [--phi-list "0 1 2 3 4 5 6 7 8 9 10 11 12"]
"""
import argparse
from pathlib import Path

import numpy as np


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdbid", required=True)
    ap.add_argument("--root", default="__SSI_ROOT__")
    ap.add_argument("--t-start", type=float, default=1000.0,
                    help="Discard data before this time in ps (default: 1000)")
    ap.add_argument("--phi-list", default="0 1 2 3 4 5 6 7 8 9 10 11 12")
    return ap.parse_args()


def phi_to_dirname(phi):
    whole = int(phi)
    frac  = int(round((phi - whole) * 10))
    return f"phi_{whole:02d}p{frac}"


def read_colvar(colvar_path, t_start):
    """Read COLVAR file, return arrays of (time, N, Ntilde, bias)."""
    times, N, Nt, bias = [], [], [], []
    with open(colvar_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            t = float(parts[0])
            if t < t_start:
                continue
            times.append(t)
            N.append(float(parts[1]))
            Nt.append(float(parts[2]))
            bias.append(float(parts[3]))
    return (np.array(times), np.array(N),
            np.array(Nt), np.array(bias))


def main():
    args    = parse_args()
    root    = Path(args.root)
    pdbid   = args.pdbid.upper()
    run_dir = root / "runs" / pdbid
    ana_dir = run_dir / "analysis"
    ana_dir.mkdir(parents=True, exist_ok=True)

    phi_list = [float(x) for x in args.phi_list.split()]

    print(f"=== N vs phi analysis for {pdbid} ===")
    print(f"t_start = {args.t_start} ps\n")

    results = []

    for phi in phi_list:
        dirname  = phi_to_dirname(phi)
        colvar   = run_dir / "ssi_production" / dirname / "COLVAR"

        if not colvar.is_file():
            print(f"  phi={phi:5.1f}  {dirname}  COLVAR missing — skipping")
            continue

        times, N, Nt, bias_arr = read_colvar(colvar, args.t_start)

        if len(N) == 0:
            print(f"  phi={phi:5.1f}  {dirname}  no frames after t_start — skipping")
            continue

        N_mean   = N.mean()
        Nt_mean  = Nt.mean()
        N_var    = N.var()
        n_frames = len(N)

        results.append({
            "phi":     phi,
            "dirname": dirname,
            "N_mean":  N_mean,
            "Nt_mean": Nt_mean,
            "N_var":   N_var,
            "n_frames": n_frames,
        })

        print(f"  phi={phi:5.1f}  {dirname}  "
              f"<N>={N_mean:8.2f}  <Nt>={Nt_mean:8.2f}  "
              f"var(N)={N_var:8.2f}  frames={n_frames}")

    if not results:
        print("No results found.")
        return

    # Find phi* = the susceptibility maximum within the collective dewetting
    # transition region.
    #
    # Naively taking the global max of var(N) can land on a secondary
    # transition at high phi; naively taking the first local max can land on
    # a tiny noise bump at low phi (before any dewetting). We therefore
    # restrict the search to the transition window where <N> has dropped to
    # between 20% and 90% of N0 -- i.e. where collective dewetting is actually
    # happening -- and take the maximum of var(N) there.
    vars_ = [r["N_var"] for r in results]
    means_ = [r["N_mean"] for r in results]
    N0 = results[0]["N_mean"] if results[0]["phi"] == 0.0 else means_[0]

    lo, hi = 0.20 * N0, 0.90 * N0
    candidates = [i for i, m in enumerate(means_) if lo <= m <= hi]

    if candidates:
        # max var(N) within the transition window
        phi_star_idx = max(candidates, key=lambda i: vars_[i])
    else:
        # fallback: global maximum of susceptibility
        phi_star_idx = int(np.argmax(vars_))

    phi_star = results[phi_star_idx]["phi"]

    print(f"\nphi*  = {phi_star:.1f} kJ/mol  "
          f"(max susceptibility var(N) = {results[phi_star_idx]['N_var']:.2f})")
    if N0 is not None:
        norm_susc = results[phi_star_idx]["N_var"] / N0
        print(f"N0    = {N0:.2f}")
        print(f"Normalized susceptibility = {norm_susc:.4f}")

    # Write summary TSV — mark phi* row
    out_tsv = ana_dir / "N_vs_phi.tsv"
    with open(out_tsv, "w") as fh:
        fh.write("phi\tdirname\tN_mean\tNt_mean\tN_var\tn_frames\tphi_star\n")
        for r in results:
            is_star = "1" if r["phi"] == phi_star else "0"
            fh.write(f"{r['phi']}\t{r['dirname']}\t"
                     f"{r['N_mean']:.6f}\t{r['Nt_mean']:.6f}\t"
                     f"{r['N_var']:.6f}\t{r['n_frames']}\t{is_star}\n")
    print(f"\nWrote {out_tsv}")
    print(f"\nUse phi* = {phi_star:.1f} for dewetting map.")
    print(f"  python3 compute_nS.py --pdbid {pdbid} --phi {phi_to_dirname(phi_star)}")
    print(f"  python3 combine_dewetting_map.py --pdbid {pdbid} --phi {phi_to_dirname(phi_star)}")


if __name__ == "__main__":
    main()
