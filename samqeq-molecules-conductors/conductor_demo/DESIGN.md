# Central result: molecules and a conductor in one solve (design, 2026-09-27)

Tod: "design a clean system to demonstrate the self-energy conductor limit. this will be the central result in
the paper." Claim to demonstrate: **conventional QEq is well posed for neither limit; samQEq is well posed for
both, in one cell and one solve.** Everything is a single-point solve on a fixed geometry (no dynamics), one
rank, the release binary (`~/codes/lammps/lammps-release/build/lmp`), so every number is deterministic and
reproducible from this directory.

## Model (no new parameters)
- **Metal:** the Au(111) slab of the paper's §2.4 study (`ms_alphasigma_overleaf/conductor_slab`): 5 layers of
  48, 17.651 × 20.381 Å, layers at z = 10.0 … 19.608 Å, Gaussian kernel R_c = 1.618 Å, χ = 5.31 eV,
  η = 5.172 eV (physical Thomas–Fermi value), self-energy width w = 0.5 Å (Scalfi).
- **Molecules:** TIP4P-FQ water, the Gaussian-kernel port shipped in `examples/tip4pfq_liquid`
  (`water4_liquid_gauss.param`; O uncharged, H and M fluctuate), geometry cut from its equilibrated 256-molecule
  liquid (`tip4pfq_liquid.data`).
- One kernel (Gaussian) for every site, as the §2.4 positivity argument requires. `boundary p p f`,
  `kspace_modify slab 3.0`, `pppm/samqeq`, solve tolerance 1e-8.

## Experiment I — one cell, both limits (the 2×2)
Au slab + a water film (≈ 7 Å thick, ~80 molecules) on the top face, closest O ≥ 3.0 Å above the top Au plane,
molecules overlapping across the new lateral boundaries removed. No applied field: the metal responds to the
water's own field.

| arm | constraint | self-energy on Au |
|---|---|---|
| **QEq** (conventional) | one global neutrality over metal + water | off |
| QEq + self-energy | one global | on |
| per-molecule, no self-energy | one per water molecule + one for the metal | off |
| **samQEq** | one per water molecule + one for the metal | on |

Scan η_Au ∈ {5.172, 2.5, 1.0, 0.3, 0.05} (the parser requires η > 0; 0.05 eV stands for the perfect-conductor limit).
Observables per arm and η:
1. λ_min of the constrained operator (code's 200-step Lanczos, as in §2.4);
2. charge transferred between water and metal, Σq(water) — zero by construction under per-molecule constraints,
   free under the global one ("molecules treated as metals");
3. spread of per-molecule net charge (max |Σq_m|) and max |q_i|;
4. mean water dipole relative to the same film with no metal (the molecular limit must be undisturbed);
5. layer charges of the slab (the image response to the water film).

Expected: only the samQEq arm is positive definite at every η **and** keeps every molecule neutral; the global
arms move charge between water and metal; the arms without the self-energy go indefinite below η ≈ 3.5 eV.

## Experiment II — a conductor must obey the image law
Au slab replicated 3 × 3 laterally (2160 atoms, 52.95 × 61.14 Å) so that single-image physics holds out to
z ≈ 10 Å. A fixed +1 point charge above the top face and a fixed −1 below the bottom face (neutral cell; both
outside the solve group, entering through the fixed-charge field), at heights z = 2 … 12 Å from the outer
planes. Arms: self-energy on / off, η_Au ∈ {5.172, 0.3}.
Observables: charge induced on the near face (→ −q for a conductor), and the interaction energy
E(z) − E(z_ref) against the classical image law −k_e q² / [4 (z − z₀)], z₀ fitted (image-plane position).

## Deliverables
`build.py` (geometry + decks), `run_*.sh`, `analyze.py` → `RESULT_conductor_demo_*.md` + a figure script. The
paper's central figure: (a) λ_min vs η for the four arms; (b) water→metal charge transfer vs η; (c) image-law
energy vs height with and without the self-energy.
