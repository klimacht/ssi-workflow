# Phl p 5a Y65 → L-DOPA: free MD of WT vs. mutant

Unrestrained MD to see how hydroxylation of Tyr65 to L-DOPA
(3,4-dihydroxyphenylalanine, PDB residue code `DAH`) changes the structure of
the grass pollen allergen Phl p 5a. Context: site-selective Tyr hydroxylation
by peroxynitrite (Bothen et al., *Free Radic. Res.* 2026,
doi:10.1080/10715762.2026.2752377).

Unlike the SSI runs in the main workflow (heavy atoms restrained), everything
here is free MD, including the wild-type reference.

## Design

| Choice | Value | Why |
|--------|-------|-----|
| Structure | AF3 model (`examples/PHLP5/PHLP5_input.pdb`) | no complete experimental structure |
| System | helical domain 1, residues 34–144, ACE/NME caps | termini and linker have pLDDT < 50; full length would be ~800k atoms in a cube vs. ~35–40k; every residue within 6 Å of Y65 lies in 61–127, the domain touches the rest of the model at one residue (222) |
| Variants | `wt`, `dah_ce1`, `dah_ce2` | DOPA is not symmetric: the new 3-OH sits on either ring edge; buried rings flip on µs–s time scales, so both orientations are simulated |
| Force field | amber99sb-ildn + TIP3P | validated for folded proteins (Lindorff-Larsen 2010, 2012); shipped with GROMACS |
| Ions | neutralize + 0.15 M NaCl | Y65 neighbours D63, K64, R66, E116, K123 |
| Box | rhombic dodecahedron, 1.2 nm | ~71 % of a cube's volume |
| Equilibration | EM, 100 ps NVT, 1 ns NPT, heavy atoms restrained | |
| Production | 4 replicas × 300 ns per system, NPT 300 K / 1 bar | V-rescale, C-rescale, PME, 2 fs, h-bond constraints; 4 × 32 cores = one Viper node (`mdrun -multidir`) |

Note from building the variants: in `dah_ce2` the new O starts 2.5 Å from the
ring of Y127 (in `dah_ce1` the nearest atom is 4.7 Å away), i.e. that ortho
carbon is sterically blocked in the model. Energy minimization removes the
clash; compare with the atom-level ASA of CE1/CE2 in Table 1 of the paper.

## DOPA parameters

Bonded terms: all from parm99 (Wang et al. 2000). Both hydroxylated ring
carbons are typed `CA` (Tyr uses `C` for CZ, but `C–CA–OH` is not defined),
so `CA–OH`, `CA–CA–OH`, `CA–OH–HO`, `X–CA–OH–X` and the improper
`CA–CA–CA–OH` cover the catechol.

Charges: RESP with the ff94/ff99SB recipe (Bayly 1993; Cornell 1995;
Cieplak 1995), `params/resp_dopa.py`:

- Ace-DAH-NMe, backbone α (φ, ψ = −60°, −40°) and β (−120°, 140°), φ/ψ
  constrained, both intramolecular catechol H-bond patterns → 4 conformers
- HF/6-31G* (Cartesian d functions, as Gaussian's 6-31G* used for ff94)
  optimization and ESP with PySCF (Sun et al. 2020), φ/ψ frozen with
  geomeTRIC (Wang & Song 2016); venv: `params/make_qm_venv.sh`
- Merz–Kollman shells 1.4/1.6/1.8/2.0 × vdW radius, 6 points/Å²
- two-stage RESP (a = 0.0005 / 0.001, b = 0.1), backbone N, H, C, O fixed to
  ff99SB values, residue and caps neutral
- check: the same protocol applied to Tyr is compared with the ff99SB Tyr
  charges (`params/results/resp_report.txt`)

The fit writes `ff/dah.rtp` and `ff/dah.hdb`; `ff/make_ff.sh` appends them to a
copy of the GROMACS force field (`amber99sb-ildn-dopa.ff`) and adds `DAH` to a
local `residuetypes.dat` in each run directory.

## Running

Configured like the main workflow (`setup.sh` replaces `__SSI_ROOT__`;
`config/workflow.env` provides `GMX` and modules). Settings: `dopa/dopa.env`.

```bash
bash dopa/params/make_qm_venv.sh /ptmp/$USER/venvs/qm   # once (DOPA_QM_VENV)
bash dopa/bin/dopa_submit.sh resp            # QM array + fit (once; results are committed)
bash dopa/bin/dopa_submit.sh setup wt        # build + equilibrate one system
bash dopa/bin/dopa_submit.sh bench wt        # 20k-step benchmark -> rep*/bench.log
bash dopa/bin/dopa_submit.sh all             # setup -> production for all systems
bash dopa/bin/dopa_submit.sh analysis        # after production
```

Production (`13_production.slurm`) runs 12 h jobs (cluster walltime) with
`-maxh 11.5` and resubmits itself from the checkpoints until all replicas
are done (at most `DOPA_MAX_RESUBMIT` times).

Runs go to `runs/dopa/<system>/rep<i>/`, the comparison to `runs/dopa/analysis/`:

| File | Content |
|------|---------|
| `summary.tsv` | per metric and system: mean ± SEM over replicas, Welch p vs. WT, per-replica values |
| `per_residue.tsv` | Cα RMSF and DSSP helicity per residue |
| `contacts.tsv` | side-chain contact occupancy of residue 65 (≥ 5 % in any system) |
| `hbonds.tsv` | H-bonds of the OH groups of residue 65 (donor and acceptor) |
| `rmsf_helicity.png`, `rmsd_time.png`, `chi1_chi2.png` | overview plots |

Metrics for residue 65: χ1 rotamer populations, ring-flip fraction
(|Δχ2| > 90° from the start), ring-centroid distances to Y112, Y124, Y127,
side-chain SASA, contacts, H-bonds. The first 50 ns of each replica are
discarded (`DOPA_SKIP_NS`).

## Limitations

- one domain only: possible transient contacts with the N-terminus or domain 2 are missing
- 300 ns covers local rearrangements (side chains, H-bond networks, helix ends), not unfolding;
  a stability change (ΔΔG) would need alchemical free-energy calculations (e.g. Gapsys et al. 2016)
- DOPA is modelled in its reduced, neutral form (catechol pKa ≈ 9–10); oxidation to dopaquinone is not included
- n = 4 replicas per system: p-values are indicative only

## References

- Bayly, Cieplak, Cornell, Kollman, *J. Phys. Chem.* 1993, 97, 10269 (RESP)
- Cornell et al., *J. Am. Chem. Soc.* 1995, 117, 5179 (ff94)
- Cieplak, Cornell, Bayly, Kollman, *J. Comput. Chem.* 1995, 16, 1357 (multi-conformation RESP)
- Wang, Cieplak, Kollman, *J. Comput. Chem.* 2000, 21, 1049 (parm99)
- Sun et al., *J. Chem. Phys.* 2020, 153, 024109 (PySCF)
- Wang, Song, *J. Chem. Phys.* 2016, 144, 214108 (geomeTRIC)
- Hornak et al., *Proteins* 2006, 65, 712 (ff99SB)
- Lindorff-Larsen et al., *Proteins* 2010, 78, 1950 (ILDN); *PLoS ONE* 2012, 7, e32131
- Bussi, Donadio, Parrinello, *J. Chem. Phys.* 2007, 126, 014101; Bernetti, Bussi, *J. Chem. Phys.* 2020, 153, 114107
- Knapp, Ospina, Deane, *J. Chem. Theory Comput.* 2018, 14, 6127 (replicas)
- Kabsch, Sander, *Biopolymers* 1983, 22, 2577 (DSSP)
- Weininger et al., *Biochemistry* 2014, 53, 4519 (ring flips)
- Gapsys et al., *Angew. Chem. Int. Ed.* 2016, 55, 7364 (ΔΔG by free-energy calculations)
