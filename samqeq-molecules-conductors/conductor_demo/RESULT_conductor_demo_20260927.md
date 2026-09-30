# RESULT — molecules and a conductor in one solve (2026-09-27)

Design: `DESIGN.md`. Binary: `~/codes/lammps/lammps-release/build/lmp` (md5 `55e1aaf5`, the release tree at capsule
patch `a7b41856`). Every run is a single-point solve, one rank, tolerance 1e-8. Reproduce:
`python3 build.py && python3 run_exp1.py && python3 run_exp2.py && python3 analyze_exp1.py && python3 analyze_exp2.py`;
figure `../figures/fig_conductor.py`. All 41 + 36 runs completed (no failures after the exp2 z-offset fix below).

## Experiment I — Au(111) + 81 TIP4P-FQ waters, the 2×2 (constraint × self-energy), η_Au scanned

Metal: the §2.4 slab (240 Au, 5 layers, Gaussian R_c 1.618 Å, χ 5.31 eV), self-energy width 0.5 Å. Water: shipped
Gaussian TIP4P-FQ port, film 22.1–30.0 Å (lowest O 3.0 Å above the top Au plane; closest water–Au 2.99 Å). No
field. η_Au = 0.05 eV stands for the conductor limit (the parser requires η > 0). No-metal reference film: mean
dipole 2.3486 D (about O, minimum image; the film's value, not the liquid's). **Corrected 2026-09-27**, see ledger.

| arm | η_Au | λ_min (eV) | CG | Q_water (e) | max \|Q_mol\| | max \|q\| | μ/μ₀ |
|---|---|---|---|---|---|---|---|
| global, no self-energy (QEq) | 5.172 | +1.530 | 70 | +1.491 | 0.279 | 1.36 | 1.037 |
| | 2.5 | **−1.020** | 593 | −3.847 | 0.428 | 1.65 | 0.973 |
| | 1.0 | −2.520 | 338 | +1.289 | 0.490 | 1.36 | 1.029 |
| | 0.3 | −3.219 | 368 | +1.507 | 0.273 | 1.36 | 1.034 |
| | 0.05 | −3.469 | 541 | +1.605 | 0.232 | 1.36 | 1.034 |
| global, self-energy | 5.172 → 0.05 | +1.531 (all) | 58–60 | +0.89 → +0.97 | 0.267–0.269 | 1.36 | 1.028–1.029 |
| per molecule, no self-energy | 5.172 | +1.654 | 41 | 0 | 0 | 1.29 | 1.001 |
| | 2.5 | −1.015 | 309 | 0 | 0 | **3.01** | 1.033 |
| | 1.0 → 0.05 | −2.515 → −3.465 | 154–190 | 0 | 0 | 1.29–1.35 | 0.989–1.016 |
| **per molecule, self-energy (samQEq)** | 5.172 → 0.05 | **+1.681 (all)** | 31–35 | **0** | **0** | 1.28 | **1.000** |

Slab layer charges (bottom → top): samQEq at every η a smooth image response, top +0.10…+0.13, sub-surface
−0.08…−0.11, converging as η → 0; the arms without the self-energy are staggered and, at η = 2.5 under the global
constraint, reach −15.4/+29.6/−5.4/−27.4/+22.5 e. Full table: `exp1.csv`.

Reading:
1. **Each ingredient repairs one limit.** Without the self-energy the operator's smallest eigenvalue follows the
   bare slab, λ_min ≈ η − 3.47 eV (the §2.4 value), whatever the constraint; the solve then takes 150–600
   iterations and returns charges that are not a minimum. Under the global constraint 0.9–1.6 e (3.8 e at η = 2.5)
   moves between the water and the metal and single molecules carry up to 0.49 e: the molecules are treated as
   part of the metal, **even with the self-energy on** (Q_water 0.89–0.97 e, the operator definite).
2. **Only the combination is well posed.** samQEq's λ_min is 1.681 eV at every η, set by the water, not the metal
   (the metal's modes sit above 13 eV with the self-energy, §2.4); CG needs 31–35 iterations at every η; every
   molecule is neutral; the water dipole is that of the metal-free film to 0.1%.

## Experiment II — the image charge: 2160-atom slab (3 × 3), fixed ±1 point charges at height z

Slab shifted +15 Å in z (the first set placed the lower charge outside the box at z = 12 Å, "Did not assign all
atoms correctly"; all runs repeated). Induced charge on the near half (top two layers; the middle layer excluded):

| arm | η_Au | induced charge, z = 2 → 12 Å | CG |
|---|---|---|---|
| self-energy | 5.172 | −1.215 → −1.234 (flat from 3 Å) | 18–21 |
| self-energy | 0.3 | −1.237 → −1.241 (flat) | 21–23 |
| no self-energy | 5.172 | −0.936 → −0.687 (falls with distance) | 46–50 |
| no self-energy | 0.3 | −2.53 → −3.19 (indefinite) | 289–317 |

With the self-energy the slab screens the whole unit charge at every distance and both hardnesses, as a conductor
must; the ~23 % excess over −1 is the same over-screening measured for the field in §2.4 (face charge 1.06–1.16×
Gauss) and is not decomposed here (the charges are nuclear-site bookkeeping of overlapping Gaussian densities).

**The energy is not a discriminating test here.** E(z) fitted to the single image law plus the uniform-sheet term
of the lateral image lattice, 2πk_e q²(z − z₀)/A per charge, leaves 30–60 meV rms structured residuals in *every*
arm, including the indefinite one (z₀ ≈ 0–0.3 Å with the self-energy, 1.2–1.4 Å without). Not used as a result;
a clean energy test needs the exact periodic image sum (or a much larger cell).

## Caveats
- One water configuration (a cut of the equilibrated TIP4P-FQ liquid; molecules overlapping across the new lateral
  boundaries removed), no dynamics. The result is a property of the operator, which is what it is meant to show.
- λ_min is the code's 200-step Lanczos estimate (an upper bound), as in §2.4.
- μ is taken about the O site; under the global constraint molecules are charged, so their μ is origin-dependent.

## Ledger
- **2026-09-27 · dipole analysis bug (found by the referee agent, confirmed).** `analyze_exp1.py` took H/M positions
  relative to O without the lateral minimum image; the dump is wrapped and 5 of 81 molecules straddle the boundary.
  The film dipole was reported as 5.2006 D; the correct value is **2.3486 D**. Corrected ratios μ/μ₀: samQEq **1.001**
  at every η (the paper's "0.1 %" stands); QEq 0.964–1.051; the μ/μ₀ column in the table above predates the fix —
  use `exp1.csv` (regenerated). Charges, λ_min, Q_water and layer charges are unaffected.
- **2026-09-27 · §7.1 misquoted the definite arm.** "Even where it is definite, 1.3–1.6 e … up to 0.5 e" drew on the
  indefinite points; conventional QEq is definite only at η = 5.172, where Q_water = 1.49 e and max|Q_mol| = 0.28 e.
  Text corrected.

## Referee item 3 additions (2026-09-27, night)

**Dense spectra** (item 2): `dense_exp1.csv` — see REVISION_referee_20260927.md §2.

**Electronegativity offset (P1)**, `run_offset.py` → `offset.csv`. H and M χ shifted together by Δχ:
- global, no self-energy (η 5.172): Q_water = +3.93, +2.71, +2.10, +1.49, +0.88, +0.27, −0.95 e for Δχ = −2 … +2 eV —
  linear, zero near Δχ ≈ +1.22 eV; max|Q_mol| 0.39 → 0.19 e, **never below 0.19 e** (charge moves between molecules at
  every offset, including where the net transfer vanishes).
- global, self-energy (η 5.172 / 0.05): same linear behaviour, zero near +1.3 eV; max|Q_mol| ≥ 0.20 e.
- per molecule (samQEq): Q_water = 0 and the charges are invariant to the offset to ≤ 2.2e-11 e.

**Screening by densities and the exact image law (P3)**, `analyze_exp2b.py` → `exp2b.csv`. Charge above the mid-plane
from Gaussian densities of the pair-kernel width (σ = 2.378 Å):
- self-energy: −1.13 (η 5.172) / −1.16 e (η 0.3), flat beyond 4 Å (sites gave −1.23/−1.24).
- no self-energy: −1.07 to −1.10 e at η 5.172 (definite: a legitimate TF metal, slightly distance-dependent);
  −1.9 e at η 0.3 (indefinite).
- The §2.4 field test on the same definition (`analyze_slabfield.py`): self-energy 1.13–1.16 × Gauss at every η;
  no self-energy 1.06 at η 5.172 but −0.63 … 2.24 × Gauss below J(0). **Both tests now give 13–16 % over-screening.**
- Exact lateral-lattice image energy, E = (k_e q²/2)(2π/A)[d − Σ_G e^{−Gd}/G], d = 2(z − z₀), fitted for z ≥ 4 Å:
  self-energy rms 8.1 / 10.3 meV, image plane z₀ = 0.47 / 0.72 Å above the top atoms; no self-energy 31 meV (η 5.172)
  and **10.7 meV at η 0.3 (indefinite)** — the energy is not a discriminating test, so it is reported but not used.

**Self-energy width (P4)**, `run_wscan.py` → `wscan.csv`, w = 0.3 … 2.378 Å (own width), η 5.172 and 0.3:
- field-test λ_min = Gram bracket − 0.11 eV at every w (the bound, confirmed across the scan).
- sub-surface layer carries the face's sign (a screening tail) for w ≤ 0.5 Å at both η; it flips (staggered mode) from
  w ≈ 1.0 Å (η 5.172) / 0.75 Å (η 0.3); at the own width and η 0.3 the layers reach ∓3.1/3.2 e and the face 0.81 × Gauss.
- density face ratio 1.06–1.18 for w ≤ 1.5 Å; water λ_min 1.68 eV and dipole within 0.2 % for w ≤ 1.5 Å; at w ≥ 2 Å
  (η 0.3) the metal's softest mode falls below the water's (0.85, 0.20 eV).
- Recommendation: w ≤ 0.5 Å (bracket ≥ 13 eV, monotone profile); 0.5 Å (Scalfi's value) is used throughout.

**Dynamics — the Gaussian TIP4P-FQ port is marginal on its own (found 2026-09-27).** At 1 fs, the samQEq arms engaged the
adaptive ridge within ~0.8 ps with M-site charges −1.6 to −1.8 e **3–7 Å from the metal**; the same film **without gold**
did the same (ridge engaged from step 813; max|q| 2.1–2.6 e), so it is the water model, not the conductor. The shipped
bulk example (as distributed, 0.2 fs) engages the ridge once in 1 ps. Consistent with the item-2 bound: the port's M site
has η 16.11 < J_MM(0) 17.53 eV. The 1 fs runs were stopped; dynamics rerun at 0.2 fs with a metal-free control
(`run_dyn.py`, runs/dyn/).

## Ledger (2026-09-28 night) — a force/solve substitution in MY decks, and a gap in the package's check

- **Retracted: "the Gaussian TIP4P-FQ port is marginal in dynamics" (written above, 09-27).** Every deck I wrote tonight
  (exp1, exp2, wscan, run_dyn, run_dyn_spc) gave `coul/shield/intra` coefficients per TYPE (`pair_coeff 1 1 …`,
  `2 2 …`) only. Under `pair_style hybrid/overlay` LAMMPS does not mix a sub-style across types when `* * lj/cut` and
  `* * coul/long` already cover the cross pairs, so every CROSS pair (O–H, H–M, Au–water, Au–ion) silently lost its
  shielding correction in the FORCES and ENERGY, while the charge SOLVE (the fix's own kernel) kept it. Single point on
  the SPC-FQ bulk box: charges identical (max|dq| 0), pe −208.41 vs −881.88 eV, max|dF| 6.63 eV/Å
  (`runs/sp_pairall` vs `runs/sp_pertype`). The dynamics "instabilities" were this inconsistency, not the water model;
  the metal-free film "control" had it too. The shipped examples list every pair explicitly (tip4pfq_liquid) or use
  `* *` (SPC-FQ), so they are unaffected.
- **What survives:** all charge-derived results (exp1 charges, dense λ_min, Q_water, offset scan, width scan, exp2
  induced/density charges) are unaffected — the solve does not read the pair coefficients. **Affected:** every ENERGY
  from these decks: the exp2 image-energy fit (Au–ion cross term missing) must be redone; the exp1 pe values were never
  used.
- **Package gap:** `check_shield_consistency()` compares kernel, range, λ and scope but not per-type-pair COVERAGE, so it
  does not refuse a hybrid/overlay deck in which some solved type pair has no `coul/shield/intra` term. To be closed in
  the code (and documented).
- **Also found while bisecting (package usability, referee P10):** the shipped native deck pins `recip_self 0.37933`;
  changing only its cutoffs (5/10 → 12/12 Å) makes PPPM retune g_ewald, the pinned self-term is then wrong, and the
  liquid goes near-critical (T 926 K, max|q| 3.8 e, 28 ridge engagements, `runs/native_check/v_cut12`); calibrated
  instead, the same deck is stable (`v_calib`, `v_mine`).
- Not yet understood: `rigid/npt/small` compression of the 30 Å box (z only) dies with "Out of range atoms - cannot compute
  PPPM" at ~9.5k steps, pinned or automatic mesh (`runs/spc_bulk180/in.npt_attempt.lammps`). Worked around by direct
  insertion at liquid density; open for the samQEq session.

## Referee item 3 closed (2026-09-28 morning) — dynamics and the image energy with every pair covered

**Dynamics** (`run_dyn.py` → `runs/dyn/`, harvest `python3 analyze_dyn.py` → `dyn.csv`). Frozen gold, rigid TIP4P-FQ,
rigid/nvt/small 298 K, 0.2 fs × 25 000 = 5 ps, tol 1e-8, all 10 coul/shield/intra pairs explicit:

| arm | max\|q\| (e) | ridge | CG fail | Q_water (e) | max\|Q_mol\| (e) | sd(q_Au) (e) |
|---|---|---|---|---|---|---|
| samQEq, η 5.172 | 1.325 | 0 | 0 | 0 | 0 | 0.0034 |
| samQEq, η 0.3 | 1.332 | 0 | 0 | 0 | 0 | 0.0048 |
| per molecule, no self-energy, η 0.3 | **1.882, rising** (1.34 → 1.79 → 1.88) | 0 | 0 | 0 | 0 | **0.63** |
| global QEq, η 5.172 | 1.362 | 0 | 0 | +3.73 mean | 0.28 | 0.051 |
| global QEq, η 0.3 | 3.975 | **≥ 8600 from step 1343 (0.27 ps)** | 4 | +4.63 mean, 9.49 max | 1.02 | 0.20 |
| metal-free film | 1.344 | 0 | 0 | 0 | 0 | — |

The per-molecule arm without the self-energy is the silent failure: no guard fires (max|q| below the threshold), but the
gold's staggered mode grows with time. LAMMPS stopped printing warnings after 100, so the ridge count is a floor.

**Image energy with the Au–ion cross pair** (`runs/exp2/`, old per-type decks in `runs/exp2_nocross_20260927`; old
CSVs `exp2_nocross.csv.bak`, `exp2b_nocross.csv.bak`). Charges unchanged (the solve never read the pair coefficients),
so Fig 3c is unchanged. Exact-law fit z ≥ 4 Å: self-energy rms 22.9 / 27.2 meV, z0 = 0.21 / 0.44 Å (was 8.1/10.3 meV,
0.47/0.72 Å); no self-energy 65.1 meV (η 5.172), 18.2 meV (η 0.3, indefinite). Below 4 Å the ion's Gaussian (0.5 Å)
overlaps the metal densities and the point-charge law does not apply (+0.3 to +1.7 eV at 2–3 Å). Conclusion unchanged:
the energy does not discriminate. §7.1 updated.
