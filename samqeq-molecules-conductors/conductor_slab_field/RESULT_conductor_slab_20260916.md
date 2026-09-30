# RESULT — conductor limit on the Au(111) slab, current binary (2026-09-16)

For α·σ highlight H1 ("molecules, insulators and conductors"). Plan: `../PLAN_restructure_highlights_20260916.md` §8.

## System and protocol

Regression case `~/codes/samQEq/internal/tests/cases/metal_slab_gself`, copied here unchanged: 240-atom
orthorhombic Au(111), 5 layers of 48, box 17.651 × 20.381 × 60 Å, `boundary p p f`, PQEq Gaussian kernel
(R_c = 1.618 Å, λ = 0.462770), `pppm/samqeq 1e-6` with `kspace_modify slab 3.0 mesh 24 28 80 gewald 0.35`,
`fix efield 0 0 0.5` (V/Å), global neutrality (one fragment), solve tolerance 1e-8, `run 0`, one rank.
χ_Au = 5.31 eV. The arm with the self-energy uses `fix_modify chg gself on width 0.50` (16.248 eV/e² added
to η).
Binary: `~/codes/lammps/lammps-dev2/build/lmp_parallel` (2026-09-12). Its `src/SAMQEQ` is byte-identical to
samQEq `deb4a08` on the five core files. **Control:** the unmodified deck reproduces the committed golden
exactly (pe 2.35068898730608, max per-atom charge difference 0.0), in `golden_check/`.

Commands: `bash run_sweep.sh && python3 analyze.py` (→ `layers.csv`); `bash run_lmin.sh` (→ `lmin.txt`).

## Layer charges (e), Gauss limit ε₀EA = 0.9940 e per face

| gself | η (eV) | L1 | L2 | L3 | L4 | L5 | face/Gauss | λ_min (eV) |
|---|---|---|---|---|---|---|---|---|
| off | 5.172 | −2.1707 | +1.4896 | −0.0224 | −1.4479 | +2.1514 | 2.174 | +1.699 |
| off | 2.5 | −9.1927 | +13.5358 | +0.3830 | −14.0396 | +9.3134 | 9.309 | −0.9733 |
| off | 1.0 | +2.3986 | −6.5464 | −0.0645 | +6.6176 | −2.4053 | 2.416 | −2.473 |
| off | 0.3 | +1.0284 | −4.1856 | −0.0310 | +4.2173 | −1.0291 | 1.035 | −3.173 |
| on | 5.172 | −1.0562 | −0.1690 | +0.0001 | +0.1703 | +1.0548 | 1.062 | +17.95 |
| on | 2.5 | −1.1048 | −0.1270 | +0.0001 | +0.1285 | +1.1032 | 1.111 | +15.27 |
| on | 1.0 | −1.1371 | −0.0954 | +0.0000 | +0.0972 | +1.1352 | 1.143 | +13.77 |
| on | 0.3 | −1.1538 | −0.0781 | −0.0000 | +0.0800 | +1.1518 | 1.160 | +13.07 |

The total charge is zero to ≤ 6e-12 in every arm. λ_min is the code's own 200-step Lanczos estimate of the
projected operator (the `ridge eig 1000 200` diagnostic runs in `runs/*/lmin/`; their charges are not used).

## What it shows

1. **With the self-energy the slab screens like a conductor at every η**: monotone profile, the sub-surface
   layer carrying the sign of its face (a screening tail), zero interior, and face charge 1.06–1.16× the
   Gauss limit, rising slowly as η → 0. The ON numbers match `NOTE_VALIDATE.md` (08-23) to every printed digit.
2. **Without it the operator is indefinite below η ≈ 3.5 eV.** λ_min falls linearly with η at unit slope,
   λ_min ≈ η − 3.47 eV, so the lowest eigenvalue of the off-diagonal lattice sum is ≈ −3.47 eV. With the
   self-energy, λ_min ≈ η + 12.78 eV, positive for every η ≥ 0. That is expected. With the self-energy on the
   diagonal the full Gaussian-charge Coulomb matrix is the energy of a real charge density, and so is
   positive semi-definite.
3. **CG reports convergence on the indefinite operators** (10–14 iterations, residual/b ≤ 8.5e-9) and
   returns unphysical charges: 9.3× Gauss at η = 2.5, and **inverted face charges** at η ≤ 1.0 (the field
   pulls positive charge to the face it should make negative). This is the same reading error as the
   paper's silent-substitution class: a converged residual on a neutral charge set does not certify the answer.
4. Even at η = 5.172, where the operator is still positive definite (λ_min 1.70 eV), the response without the
   self-energy is staggered and 2.2× Gauss. The soft lattice mode dominates the response.

## The width must be narrower than the pair kernel's (`run_defaultwidth.sh`, `lmin_defaultwidth.txt`)

`gself on` with no width uses the pair kernel's own Gaussian (standard deviation 2.38 Å for Au; J_ii = 3.416 eV):

| η (eV) | L1 | L2 | L3 | L4 | L5 | face/Gauss | λ_min (eV) |
|---|---|---|---|---|---|---|---|
| 5.172 | −1.5235 | +0.4233 | −0.0026 | −0.4152 | +1.5180 | 1.530 | 5.114 |
| 0.3 | −3.2708 | +3.4525 | −0.2574 | −3.0581 | +3.1337 | 3.221 | 0.2424 |

λ_min = η − 0.058 eV at both η: the Coulomb block, taken with the self-energy of its own Gaussians, is
positive semi-definite to within 0.06 eV, as it must be for a real charge density (the residual is
presumably mesh or self-term calibration error; not isolated). Its softest mode is nearly free, though, so the
response stays staggered and grows as η → 0. A self-energy narrower than the pair kernel's lifts that mode by
the difference of the two self-energies: 16.248 − 3.416 = 12.83 eV, against the measured η + 12.78.
The 5.172 row matches `ANALYSIS_samqeq_metal_screening.md` §12.1 to every printed digit.

**Statement for the paper:** with the self-energy of a Gaussian no wider than the pair kernel's, H is positive
definite for every η ≥ 0 (up to a ≤ 0.06 eV mesh residual). Only a strictly narrower width gives the lattice a
margin, and the margin is the self-energy difference.

## Caveats

- One slab thickness, one field, one width (0.5 Å, the Scalfi et al. value). The 6–16% excess over Gauss is not
  decomposed here (finite slab, effective boundary position). State it as measured; do not attribute it.
- The legacy probe fallback fires in every arm (`recip_self: probe tag 1 has no partner within the G2 window`),
  as it does in the committed golden. Not investigated here.
- The width is a free parameter of the conductor mode; the default (tied to the PQEq R_c) is recorded in
  `~/Research/electrochem/docs/ANALYSIS_samqeq_metal_screening.md` §12.3 as insufficient (1.5× Gauss).
