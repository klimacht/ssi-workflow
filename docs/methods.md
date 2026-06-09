# Methods (concise — SSI map as a supporting analysis)

> Drop-in paragraph for a manuscript whose main subject is something else
> (e.g. PhLP5 derivatization) and that uses the SSI map only as one input.

## Surface hydrophobicity mapping by sparse sampling INDUS

The surface hydrophobicity map of PhLP5 was obtained using sparse sampling
INDUS (SSI) as described by Sinha et al. and the original approach of
Rego et al. Briefly, the protein was solvated in a cubic box (minimum
protein–wall distance 1.6 nm), neutralized with counter-ions,
energy-minimized, and equilibrated in the NVT (300 K) and NPT (1 bar)
ensembles using the Amber99SB force field and the SPC/E water model.
Following a 2 ns unbiased production run, thirteen biased simulations of
3 ns each were performed at φ = 0–12 kJ mol⁻¹ particle⁻¹, applying an
unfavourable linear potential *U = φÑ* through PLUMED to the water count
inside a union of 0.6 nm spheres centred on every protein heavy atom.
Protein heavy atoms were harmonically restrained
(k = 1000 kJ mol⁻¹ nm⁻²) and the first 1 ns of every biased trajectory
was discarded.

Rather than reporting a single dewetting value at one bias window, we
extracted the **complete hydration curve** ⟨n_i⟩_φ for every protein heavy
atom from all thirteen simulations and defined a per-atom half-dewetting
potential **φ\*_i** by cubic interpolation: the value of φ at which the
atom loses half of its unbiased hydration. The same definition applied
to the full hydration shell ⟨N⟩_φ yields the molecule-wide
half-dewetting potential **φ\*_m**. Mapping the centred quantity
**Δφ\*_i = φ\*_i − φ\*_m** onto the structure gives a continuous,
zero-centred surface descriptor: atoms whose hydration is destabilised by
less than the molecular average appear as Δφ\*_i < 0, atoms that retain
water more tightly as Δφ\*_i > 0. Surface atoms were identified by
⟨n_i⟩₀ > 4 waters within the 0.6 nm shell; buried atoms and atoms whose
hydration curve does not cross the half-dehydration level within
φ ≤ 12 kJ mol⁻¹ were excluded from the colour scale and rendered in grey.
Simulations used GROMACS 2023.3 with PLUMED 2.9.0 incorporating the INDUS
collective variable. The complete pipeline, the modified INDUS plugin
and the validation against the reference proteins of Sinha et al.
are available at [repository URL].

---

## Figure caption (suggested)

> **Figure X. Surface hydrophobicity map of PhLP5.** Per-atom dewetting
> map of PhLP5 obtained by sparse sampling INDUS, coloured by
> **Δφ\*_i = φ\*_i − φ\*_m** (kJ mol⁻¹), where φ\*_i is the bias potential
> at which a given surface atom loses half of its unbiased hydration and
> φ\*_m = 4.77 kJ mol⁻¹ is the corresponding molecule-wide value.
> **Red regions** dewet at lower bias than the protein average — water at
> these sites is weakly bound and easily displaced. **White regions**
> dewet at approximately the molecular half-dehydration potential and
> represent the typical hydration behaviour of the surface. **Blue
> regions** retain water under stronger bias than the average and
> correspond to strongly bound, structurally engaged hydration. Grey
> patches mark atoms whose hydration curve does not cross the half-
> dehydration level within the sampled φ range (extreme hydrophilicity or
> sterically protected water). The map was computed from thirteen 3 ns
> biased simulations at φ = 0–12 kJ mol⁻¹ particle⁻¹ on a 0.6 nm union of
> heavy-atom spheres. PhLP5 displays the strongest surface heterogeneity
> of all proteins examined here (Δφ\*_i p95–p05 = 4.79 kJ mol⁻¹),
> indicating well-defined hydrophobic hotspots against a comparatively
> hydrophilic background — the hallmark of a surface suitable for
> selective chemical addressing.

---

## Optional sentence on the modified INDUS code

> Add to the methods paragraph above if you publish the patched plugin.

The INDUS union-of-spheres implementation was extended with a
nearest-centre evaluation mode and configurable derivative sharing for
improved performance on large observation volumes; the modified source
and the full pipeline are provided at [repository URL].

---

## Longer Methods version (if a full Methods subsection is wanted)

See the upstream method papers for the full protocol. The key parameters
used here were: Amber99SB / SPC/E; cubic box with 1.6 nm padding; PME
electrostatics with 1.0 nm real-space cutoff; 1.0 nm van der Waals
cutoff; V-rescale thermostat (τ = 1 ps) and Parrinello–Rahman barostat
(τ = 2 ps) at 300 K / 1 bar; 2 fs timestep with LINCS-constrained
H-bonds. The observation volume was the union of 0.6 nm spheres on all
protein heavy atoms; bias windows spanned φ = 0–12 kJ mol⁻¹ in 1 kJ
mol⁻¹ steps, 3 ns each, with the first 1 ns discarded.

For each atom i and each bias φ, the time-averaged water count
⟨n_i⟩_φ within the 0.6 nm shell was obtained from the biased
trajectories. The half-dewetting potential φ\*_i was defined as the φ at
which ⟨n_i⟩_φ = ½⟨n_i⟩₀, obtained by cubic interpolation through the
four φ-values bracketing the crossing (linear interpolation at the
boundaries). The same definition applied to the full hydration shell
⟨N⟩_φ yields the molecule-wide value φ\*_m. The reported per-atom
quantity is Δφ\*_i = φ\*_i − φ\*_m so that all maps are centred at zero
and directly comparable between proteins.

Validation against the six benchmark proteins from Sinha et al. (1DPX,
1HRC, 1FS3, 1RBB, 3LDJ, 4IBA) reproduced the published ranking by the
normalised susceptibility γ = ⟨δN²⟩_φ* / N0, with 4IBA most hydrophobic
(γ ≈ 1.27) and 1HRC least (γ ≈ 0.42). The patch-contrast metric
introduced here (Δφ\*_i p95–p05) is highest for 1RBB (4.53), 1HRC (4.47)
and 1DPX (4.16) — proteins with sharply defined hydrophobic spots on
hydrophilic backgrounds — and lowest for the surface-uniform 4IBA (2.72)
and 3LDJ (2.93).
