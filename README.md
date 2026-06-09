# SSI Workflow

Automated, cluster-ready pipeline for **Sparse Sampling INDUS (SSI)** protein
surface hydrophobicity mapping, following the method of
[Sinha, Garde & Cramer, *J. Phys. Chem. B* 2023](https://doi.org/10.1021/acs.jpcb.3c04902)
and the original SSI approach of
[Rego, Xi & Patel, *PNAS* 2021](https://doi.org/10.1073/pnas.2018234118).

Given a protein structure (RCSB ID or local PDB), the pipeline runs the full
chain automatically and produces a per-atom dewetting map:

```
PDB → preparation → equilibration → unbiased production
→ biased SSI array (φ = 0–12 kJ/mol) → susceptibility → dewetting map
```

## Requirements

- GROMACS 2023.x patched with PLUMED 2.9.0 + INDUS (see `indus-patches/`)
- Python ≥ 3.9 with `numpy` and `MDAnalysis`
- Slurm workload manager

## Quick start

```bash
git clone <repo-url> ssi_workflow
cd ssi_workflow

# Configure paths (replaces placeholders, creates config/workflow.env)
./setup.sh /abs/path/to/ssi_workflow /abs/path/to/software

# Edit config/workflow.env for your cluster (module names etc.)

# Run a protein end-to-end (submits all stages with Slurm dependencies)
bash bin/ssi_run_full.sh 1DPX

# Monitor
watch -n 30 bash bin/monitor_pipeline.sh 1DPX
```

## Repository layout

```
ssi_workflow/
├── README.md
├── setup.sh                     configure placeholders + create config
├── config/
│   └── workflow.env.example     all tunable parameters (copy to workflow.env)
├── slurm/
│   ├── 00_prepare.slurm         download/clean PDB, pdb2gmx, solvate, ions, EM
│   ├── 01_equilibrate.slurm     NVT + NPT equilibration
│   ├── 02_unbiased.slurm        2 ns unbiased NPT production
│   ├── 03_ssi_preflight.slurm   geometry + INDUS/PLUMED input setup, grompp check
│   └── 04_ssi_production.slurm  biased array job (one task per φ)
├── mdp/                         GROMACS .mdp parameter files
├── src/
│   ├── setup_ssi_case.py        per-φ TPR + INDUS input + PLUMED input
│   ├── compute_n0.py            ⟨n_i⟩₀ from unbiased trajectory
│   ├── compute_nS.py            ⟨n_i⟩_φ from a biased trajectory
│   ├── compute_N_vs_phi.py      ⟨N⟩_φ and susceptibility, picks φ*
│   └── combine_dewetting_map.py η_i = nS/n0 → .dat + .pdb (B-factor = η)
├── bin/
│   ├── ssi_lib.sh               shared helper functions (sourced by all stages)
│   ├── ssi_run_full.sh          submit the whole pipeline with dependencies
│   ├── ssi_rerun_chain.sh       rerun keeping a single chain (multimeric crystals)
│   ├── run_analysis.sh          run all analysis steps after production
│   ├── monitor_pipeline.sh      multi-protein stage + φ-window progress monitor
│   ├── monitor_ssi_production.sh per-φ live progress table
│   ├── ssi_doctor.sh            environment sanity checks
│   ├── load_env.sh              source config + modules interactively
│   └── ssi_submit_*.sh          submit individual stages
├── docs/
│   └── methods.md               concise methods text for a manuscript
└── indus-patches/
    ├── README.md
    └── 0001-*.patch             INDUS modifications used here
```

## Method summary

| Stage | Duration | Ensemble | Notes |
|-------|----------|----------|-------|
| Energy minimization | — | — | steepest descent |
| NVT | — | 300 K | V-rescale, τ=1 ps |
| NPT | — | 300 K, 1 bar | Parrinello–Rahman, τ=2 ps |
| Unbiased production | 2 ns | NPT | provides the starting structure for biased runs |
| Biased production | 3 ns × 13 | NPT | φ = 0–12 kJ/mol, protein heavy atoms restrained |

The observation volume is the union of spheres (radius **0.6 nm**) centered on
every protein heavy atom. A linear biasing potential *U = φÑ* is applied through
PLUMED to progressively dewet the protein surface. The susceptibility
⟨δN²⟩_φ peaks at the characteristic dewetting potential **φ\***. For each
surface atom (⟨n_i⟩₀ > 4), the dewetting parameter

```
η_i = ⟨n_i⟩_φ* / ⟨n_i⟩₀
```

is written to the B-factor column of an output PDB (η = 0 hydrophobic →
η = 1 hydrophilic), ready to color in VMD or ChimeraX.

## Two flavors of dewetting map

Two per-atom maps can be produced; both write a single scalar per atom into
the PDB B-factor column, but use the simulation data differently.

**Single-window map** (`combine_dewetting_map.py`) — the original method.
Uses one window, phi*, and reports eta_i = <n_i>_phi* / <n_i>_0. Fast, matches
the reference paper, but discards 12 of the 13 simulations.

**Per-atom transition map** (`combine_dewetting_phistar_map.py`) — uses every
window. For each atom it builds the full hydration curve <n_i>_phi and finds
the continuous potential phi*_i at which the atom loses half its hydration
(linear interpolation between bracketing windows). Low phi*_i = dewets early =
hydrophobic; high phi*_i = holds water = hydrophilic. This uses all the data,
is continuous rather than tied to the phi grid, and is less sensitive to noise
in any single window.

```bash
# compute nS for ALL phi windows (needed once), then build the richer map
bash bin/compute_all_nS.sh PHLP5
python3 src/combine_dewetting_phistar_map.py --pdbid PHLP5
```

Output: `runs/PDBID/dewetting/dewetting_phistar.pdb` (B-factor = phi*_i).

## Example data

The repository ships with a complete set of pre-computed results for the six
benchmark proteins from [Sinha et al.][paper-jpcb] plus PhLP5, under
[`examples/`](examples/). For each protein you get the cleaned input
structure, ⟨N⟩_φ table, both dewetting maps (η at φ* and the per-atom φ*_i),
and the per-atom hydration curves. Renderings are in `examples/figures/`.

Use these to validate a local installation without rerunning 13 × 3 ns of
biased MD. See [`examples/README.md`](examples/README.md) for details.

The cross-protein summary table `examples/analysis_comparison.tsv`
reproduces the ranking from the reference paper: 4IBA most hydrophobic
overall (γ ≈ 1.27), 1HRC least (γ ≈ 0.42); strongest patch contrasts in
1RBB / 1HRC / 1DPX.

[paper-jpcb]: https://doi.org/10.1021/acs.jpcb.3c04902


## Gallery

Side-by-side comparison of the two map types for selected proteins.
The single-window η-map (left) shows hydration at the single bias φ\*; the
per-atom φ\*_i map (right) integrates information from all 13 biased
simulations into one continuous scalar per atom, centered on the protein
median. See [`examples/figures/`](examples/figures/) for all proteins and views.

| Protein | single-window η at φ\* | per-atom φ\*_i (centered) |
|---------|------------------------|-----------------------------|
| **PhLP5** — strongest patch contrast (4.78 kJ/mol) | ![](examples/figures/PHLP5_eta_front.png) | ![](examples/figures/PHLP5_phistar_front.png) |
| **1HRC** (cytochrome C) — hydrophilic on average but sharp hotspots | ![](examples/figures/1HRC_eta_front.png) | ![](examples/figures/1HRC_phistar_front.png) |
| **4IBA** (β-lactoglobulin) — most hydrophobic overall, homogeneous | ![](examples/figures/4IBA_eta_front.png) | ![](examples/figures/4IBA_phistar_front.png) |

Color scale: **red = hydrophobic**, **blue = hydrophilic**. The η-map runs
0–1 (no hydration → full hydration), the φ\*_i map is symmetric around 0 in
kJ/mol (negative = dewets earlier than the protein average).

## Handling special cases

**Custom (non-RCSB) structures** — place the file at
`runs/PDBID/pdb/PDBID.pdb` before submitting; `00_prepare.slurm` detects it
and skips the download.

**Multimeric crystal structures** — many crystal structures contain several
identical chains. To reproduce a monomer measurement, keep a single chain:

```bash
bash bin/ssi_rerun_chain.sh 3LDJ A     # keep only chain A
```

or pass the chain to prepare directly:

```bash
sbatch slurm/00_prepare.slurm 3LDJ A
```

**Incomplete termini** — residues missing backbone atoms (N, CA, C) at the
C-terminus are trimmed automatically during PDB cleaning.

## Coloring the map

```tcl
# VMD: red = hydrophobic (η→0), blue = hydrophilic (η→1)
mol modselect 0 top "beta >= 0"
mol modstyle  0 top VDW 1.0 12.0
mol modcolor  0 top Beta
color scale method RWB
color scale midpoint 0.5
```

## INDUS modifications

The INDUS PLUMED plugin was extended (see `indus-patches/`):

1. **`union_mode = nearest`** — approximate union-of-spheres indicator using
   the nearest-center distance; faster for large observation volumes.
2. **`NO_SHARE_ALL_DERIVATIVES`** — suppresses MPI broadcast of derivatives
   when biasing through an external PLUMED RESTRAINT.

## Citation

> Sinha, I.; Garde, S.; Cramer, S. M. A Comparative Analysis of Protein Surface
> Hydrophobicity Maps Determined by Sparse Sampling INDUS and Spatial
> Aggregation Propensity. *J. Phys. Chem. B* **2023**.

> Rego, N. B.; Xi, E.; Patel, A. J. Identifying hydrophobic protein patches to
> inform protein interaction interfaces. *PNAS* **2021**, *118*, e2018234118.


## scripts/ (author tools)

The `scripts/` directory contains tools used to maintain this repository
(copying results into `examples/`, rendering figures via PyMOL). They are
hardcoded to the author's cluster paths and not meant for general use.
