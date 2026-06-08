# Methods (concise — SSI map as a supporting analysis)

> Drop-in paragraph for a manuscript whose main subject is something else
> (e.g. PhLP5 derivatization) and that uses the SSI map only as one input.

The surface hydrophobicity map of PhLP5 was computed using sparse sampling
INDUS (SSI) as described by Sinha et al. Briefly, the protein was solvated in
a cubic box (minimum protein–wall distance 1.6 nm), neutralized with
counter-ions, energy-minimized, and equilibrated in the NVT (300 K) and NPT
(1 bar) ensembles using the Amber99SB force field and the SPC/E water model.
Following a 2 ns unbiased production run, thirteen biased simulations (3 ns
each) were performed at φ = 0–12 kJ mol⁻¹ particle⁻¹, applying an unfavorable
linear potential *U = φÑ* to the water count within a union of 0.6 nm spheres
centered on each protein heavy atom, with protein heavy atoms harmonically
restrained (k = 1000 kJ mol⁻¹ nm⁻²). The first 1 ns of each biased trajectory
was discarded. The characteristic dewetting potential φ\* was identified from
the maximum of the water-number susceptibility ⟨δN²⟩_φ. For each surface atom
(defined by ⟨n_i⟩₀ > 4 waters within the 0.6 nm shell), the dewetting
parameter η_i = ⟨n_i⟩_φ\* / ⟨n_i⟩₀ was mapped onto the structure (η = 0,
hydrophobic; η = 1, hydrophilic). Simulations used GROMACS 2023.3 with PLUMED
2.9.0 incorporating the INDUS collective variable; the analysis pipeline is
available at [repository URL].

---

## Optional sentence if the modified INDUS code is released

> Add to the paragraph above if you publish the patched plugin.

The INDUS union-of-spheres implementation was extended with a nearest-center
evaluation mode for improved performance on large observation volumes; the
modified source and the complete pipeline are provided at [repository URL].

---

## Longer version (if a full Methods subsection is wanted)

See the upstream method papers for the full protocol. The key parameters used
here were: Amber99SB / SPC/E; cubic box with 1.6 nm padding; PME electrostatics
with 1.0 nm real-space cutoff; 1.0 nm van der Waals cutoff; V-rescale
thermostat (τ = 1 ps) and Parrinello–Rahman barostat (τ = 2 ps) at 300 K /
1 bar; 2 fs timestep with LINCS-constrained H-bonds. The observation volume was
the union of 0.6 nm spheres on all protein heavy atoms; bias windows spanned
φ = 0–12 kJ mol⁻¹ in 1 kJ mol⁻¹ steps, 3 ns each, with the first 1 ns
discarded. Protein heavy atoms were position-restrained (k = 1000
kJ mol⁻¹ nm⁻²) throughout the biased runs.
