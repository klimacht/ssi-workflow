# SSI Workflow

Automated, cluster-ready pipeline for **Sparse Sampling INDUS (SSI)** protein
surface hydrophobicity mapping, following the method of
[Sinha, Garde & Cramer, *J. Phys. Chem. B* 2023](https://doi.org/10.1021/acs.jpcb.3c04902)
and the original SSI approach of
[Rego, Xi & Patel, *PNAS* 2021](https://doi.org/10.1073/pnas.2018234118).

Given a protein structure (RCSB ID or local PDB), the pipeline runs the full
chain automatically and produces a per-atom dewetting map:

    PDB → preparation → equilibration → unbiased production
    → biased SSI array (φ = 0–12 kJ/mol) → susceptibility
    → single-window η-map AND per-atom φ*_i-map

## Requirements & Upstream Source

- **GROMACS 2023.x**: [Source / Downloads](https://manual.gromacs.org/documentation/2023/download.html)
- **PLUMED 2.9.0**: [GitHub Release](https://github.com/plumed/plumed2/releases/tag/v2.9.0) — downloaded and patched automatically by `setup.sh`
- **INDUS** (PLUMED plugin): [Patel Lab](https://github.com/patellab511/indus) — patched locally for performance, see `indus-patches/`
- Python ≥ 3.9 with `numpy` and `MDAnalysis`
- Slurm workload manager

## Quick start

    git clone https://github.com/klimacht/ssi-workflow.git ssi_workflow
    cd ssi_workflow

    # Replaces placeholders, creates config/workflow.env, fetches and patches PLUMED.
    ./setup.sh /abs/path/to/ssi_workflow /abs/path/to/software

    # Build PLUMED (it has been configured by setup.sh):
    cd /abs/path/to/software/src/plumed-2.9.0 && make -j 4 && make install

    # Back in the workflow dir: review config/workflow.env, then run a protein
    bash bin/ssi_run_full.sh 1DPX

    # Monitor progress (works for one or several proteins)
    watch -n 30 bash bin/monitor_pipeline.sh 1DPX

## Repository layout

    ssi_workflow/
    ├── README.md
    ├── LICENSE                      MIT
    ├── CITATION.cff                 cite-this-repo metadata
    ├── setup.sh                     paths + PLUMED download/patch + config
    ├── config/
    │   └── workflow.env.example     all tunable parameters
    ├── slurm/
    │   ├── 00_prepare.slurm         download/clean PDB, pdb2gmx, solvate, ions, EM
    │   ├── 01_equilibrate.slurm     NVT + NPT equilibration
    │   ├── 02_unbiased.slurm        2 ns unbiased NPT production
    │   ├── 03_ssi_preflight.slurm   geometry + INDUS/PLUMED setup, grompp check
    │   ├── 04_ssi_production.slurm  biased array job (one task per φ window)
    │   └── 06_full_analysis.slurm   per-atom analysis on a full node, all φ in parallel
    ├── mdp/                         GROMACS .mdp parameter files
    ├── src/
    │   ├── setup_ssi_case.py        per-φ TPR + INDUS input + PLUMED input
    │   ├── compute_n0.py            ⟨n_i⟩₀ from unbiased trajectory
    │   ├── compute_nS.py            ⟨n_i⟩_φ from a biased trajectory
    │   ├── compute_N_vs_phi.py      ⟨N⟩_φ, susceptibility, picks φ*
    │   ├── combine_dewetting_map.py            single-window η-map at φ*
    │   └── combine_dewetting_phistar_map.py    per-atom φ*_i-map using all φ
    ├── bin/
    │   ├── ssi_lib.sh               shared helper functions
    │   ├── ssi_run_full.sh          submit the whole pipeline with dependencies
    │   ├── ssi_rerun_chain.sh       rerun keeping a single chain (multimers)
    │   ├── run_analysis.sh          run analysis steps after production
    │   ├── evaluate_all.sh          submit full analysis for many proteins
    │   ├── collect_comparison.sh    cross-protein comparison table
    │   ├── monitor_pipeline.sh      multi-protein stage + φ-window monitor
    │   ├── monitor_ssi_production.sh per-φ live progress table
    │   ├── ssi_doctor.sh            environment sanity checks
    │   ├── load_env.sh              source config + modules interactively
    │   └── ssi_submit_*.sh          submit individual stages
    ├── docs/
    │   └── methods.md               concise methods text for a manuscript
    └── indus-patches/
        ├── README.md
        └── 0001-*.patch             INDUS modifications

## Method summary

| Stage | Duration | Ensemble | Notes |
|-------|----------|----------|-------|
| Energy minimization | — | — | steepest descent |
| NVT | — | 300 K | V-rescale, τ=1 ps |
| NPT | — | 300 K, 1 bar | Parrinello–Rahman, τ=2 ps |
| Unbiased production | 2 ns | NPT | starting structure for biased runs |
| Biased production | 3 ns × 13 | NPT | φ = 0–12 kJ/mol, protein heavy atoms restrained |

The observation volume is the union of spheres (radius **0.6 nm**) centered on
every protein heavy atom. A linear biasing potential *U = φÑ* is applied
through PLUMED to progressively dewet the protein surface. The susceptibility
⟨δN²⟩_φ peaks at the characteristic dewetting potential **φ\***. Surface atoms
are defined by ⟨n_i⟩₀ > 4 waters within the 0.6 nm shell.

## Two flavors of dewetting map

Both maps store one scalar per atom in the PDB B-factor column but use the
simulation data differently.

**Single-window map** (`combine_dewetting_map.py`) — the original method.
Uses one window, φ*, and reports

    η_i = ⟨n_i⟩_φ* / ⟨n_i⟩₀          (0 = hydrophobic, 1 = hydrophilic)

Matches the reference paper but discards 12 of the 13 simulations.

**Per-atom transition map** (`combine_dewetting_phistar_map.py`) — uses every
window. For each atom it builds the full hydration curve ⟨n_i⟩_φ and finds
the continuous potential **φ*_i** at which the atom loses half its hydration
(linear interpolation between bracketing windows). Low φ*_i = dewets early =
hydrophobic; high φ*_i = holds water = hydrophilic. Continuous, uses all data,
robust against noise in any single window.

By default the values are *centered*: the B-factor is **Δφ_i = φ*_i − φ\*_global**,
with 0 at the protein average (white), negative for hydrophobic patches (red),
positive for hydrophilic regions (blue). Alternatives:

    --center-mode phistar_global   subtract global φ* from susceptibility (default)
    --center-mode median           subtract the per-atom median of φ*_i
    --no-center                    write raw φ*_i in kJ/mol

Both maps are produced together by a single Slurm job that recomputes all 13
nS windows in parallel on a full node:

    sbatch slurm/06_full_analysis.slurm 1DPX

For a batch over many proteins:

    bash bin/evaluate_all.sh                       # submit
    bash bin/collect_comparison.sh                 # after all finish

Outputs (per protein, in `runs/PDBID/`):

    analysis/N_vs_phi.tsv              ⟨N⟩_φ, var(N), φ*
    dewetting/dewetting_phi_XX.pdb     single-window η-map at φ*
    dewetting/dewetting_phistar.pdb    per-atom φ*_i-map (all windows)

The cross-protein table `analysis_comparison.tsv` lists, per protein:

- **γ = var(N)_φ\* / N0** — overall hydrophobicity (high = more hydrophobic)
- **patch contrast (p95 − p05 of Δφ_i)** — surface heterogeneity (high =
  sharper hydrophobic/hydrophilic patches on the same protein)

The two quantities are independent: 4IBA is overall the most hydrophobic
(γ ≈ 1.27) but has a fairly uniform surface (low contrast), while 1HRC is
relatively hydrophilic on average (γ ≈ 0.42) yet has the strongest patch
contrast — exactly the cytochrome C behaviour highlighted in the reference
paper.

## Coloring the maps

    # VMD: red = hydrophobic, blue = hydrophilic
    mol new dewetting_phistar.pdb
    mol modselect 0 top "beta > -7"        # hide buried-atom sentinel
    mol modstyle  0 top VDW 1.0 12.0
    mol modcolor  0 top Beta
    color scale method RWB
    mol scaleminmax top 0 -3 3              # symmetric kJ/mol range

For the single-window η-map use the same setup with `mol scaleminmax top 0 0 1.25`.

## Handling special cases

**Custom (non-RCSB) structures** — place the file at
`runs/PDBID/pdb/PDBID.pdb` before submitting; `00_prepare.slurm` detects it
and skips the download.

**Multimeric crystal structures** — many crystal structures contain several
identical chains. To reproduce a monomer measurement, keep a single chain:

    bash bin/ssi_rerun_chain.sh 3LDJ A     # keep only chain A

or pass the chain to prepare directly:

    sbatch slurm/00_prepare.slurm 3LDJ A

**Incomplete termini** — residues missing backbone atoms (N, CA, C) at the
C-terminus are trimmed automatically during PDB cleaning.

## INDUS modifications

The INDUS PLUMED plugin was extended (see `indus-patches/`):

1. **`union_mode = nearest`** — approximate union-of-spheres indicator using
   the nearest-center distance; faster for large observation volumes.
2. **`NO_SHARE_ALL_DERIVATIVES`** — suppresses MPI broadcast of derivatives
   when biasing through an external PLUMED RESTRAINT.

## Citation

Please cite this repository together with the original method papers:

> Sinha, I.; Garde, S.; Cramer, S. M. A Comparative Analysis of Protein Surface
> Hydrophobicity Maps Determined by Sparse Sampling INDUS and Spatial
> Aggregation Propensity. *J. Phys. Chem. B* **2023**.
> [doi:10.1021/acs.jpcb.3c04902](https://doi.org/10.1021/acs.jpcb.3c04902)

> Rego, N. B.; Xi, E.; Patel, A. J. Identifying hydrophobic protein patches to
> inform protein interaction interfaces. *PNAS* **2021**, *118*, e2018234118.
> [doi:10.1073/pnas.2018234118](https://doi.org/10.1073/pnas.2018234118)

See `CITATION.cff` for machine-readable citation metadata.
