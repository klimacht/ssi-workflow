# Example data

Complete results for the six benchmark proteins from
[Sinha, Garde & Cramer (2023)](https://doi.org/10.1021/acs.jpcb.3c04902)
plus PhLP5 (the protein this pipeline was originally built for), as produced
by the workflow in this repository. Use these to validate a local
installation without rerunning 13 × 3 ns of biased MD per protein.

## Layout

    examples/
    ├── README.md
    ├── analysis_comparison.tsv         cross-protein summary table
    ├── figures/                        PyMOL renderings (PNG)
    │   └── <PDBID>_<map>_<view>.png
    └── <PDBID>/
        ├── <PDBID>_input.pdb           cleaned input structure
        ├── prod_unbiased.gro           equilibrated reference frame
        ├── analysis/
        │   ├── N_vs_phi.tsv            ⟨N⟩_φ, var(N), φ* table
        │   └── n0_per_atom.npy         ⟨n_i⟩₀ array
        └── dewetting/
            ├── dewetting_phi_XX.pdb    single-window η-map (B-factor = η)
            ├── dewetting_phi_XX.dat    same data in tabular form
            ├── dewetting_phistar.pdb   per-atom φ*_i map (centered, all φ used)
            └── dewetting_phistar.dat   full curves ⟨n_i⟩_φ for every atom

## Cross-protein summary

`analysis_comparison.tsv` lists for every protein:

| Column | Meaning |
|--------|---------|
| `N0` | mean water count in the 0.6 nm hydration shell at φ = 0 |
| `phi_star` | collective dewetting potential (susceptibility max, kJ/mol) |
| `gamma` | var(N)_φ* / N0 — overall hydrophobicity (higher = more hydrophobic) |
| `patch_contrast_kJmol` | p95 − p05 of Δφ_i — surface heterogeneity |
| `n_surface` | number of surface atoms (⟨n_i⟩₀ > 4) |

Cross-protein ranking matches the reference paper: 4IBA most hydrophobic
overall (γ ≈ 1.27), 1HRC least (γ ≈ 0.42); the strongest patch contrasts
appear in 1RBB / 1HRC / 1DPX (proteins with clearly differentiated
hydrophobic spots on a hydrophilic background), while 4IBA and 3LDJ are
homogeneous (uniformly hydrophobic).

## Viewing a map in VMD

```tcl
mol new examples/1DPX/dewetting/dewetting_phistar.pdb
mol modselect 0 top "beta > -7"        # hide buried-atom sentinel
mol modstyle  0 top Surf 1.4
mol modcolor  0 top Beta
color scale method RWB                 # red=hydrophobic, blue=hydrophilic
mol scaleminmax top 0 -3 3              # symmetric kJ/mol range
```

## Reproducing from scratch

```bash
cd <your-workflow-checkout>
bash bin/ssi_run_full.sh 1DPX
sbatch slurm/06_full_analysis.slurm 1DPX
```

Then `runs/1DPX/dewetting/dewetting_phistar.pdb` should match the example.
The atom-by-atom values will differ slightly between runs (different random
seed in equilibration), but η-ranking and the global metrics (N0, γ,
patch-contrast) are reproducible within ~5%.
