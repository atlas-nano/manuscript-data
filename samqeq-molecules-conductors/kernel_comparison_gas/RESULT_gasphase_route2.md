# Gas-phase results recomputed on the projected-CG route (lr_ewald = 2)

2026-09-16. Author's decision of 2026-09-16: the paper's gas-phase results, produced on the `lr_ewald = 0` ACKS2
saddle, are to be recomputed on the plain projected-CG route (`lr_ewald = 2`), which is what the paper calls "the
atomic-charge solve". This directory holds every deck, log, result line and script; the manuscript was not edited.

Binary `~/codes/lammps/lammps-dev2/build/lmp_parallel`. It was **rebuilt by another session at 19:29:50** while the
alkane series was running (481188544 B, 16:05 → 481203488 B, md5 `675779c7828f37616cd54e79d49232a8`). Every result
line written after that carries a `bin=1789612190:481203488` stamp (`runs/*/res.*`). Bit test old vs new binary on
the same decks (`runs/bittest/`): C8 alkane, C12 L=200, C17 L=120, C10 bondsoft and the cbrt water monomer are
**identical to all printed digits (10 in μ, 8 in pe)**; the Slater-kernel water runs differ at 1e-6 (pe −4.539235 →
−4.539236, q_M by 2e-8), so the Slater path changed in the rebuild. Consequence: the whole water/kernel set was
regenerated on the 19:29 binary (the 16:05 set is archived in `runs/water_bin1605/`); the alkane set (cbrt kernel
path, shown unchanged) mixes the two binaries, 88 of 133 plain sweep points from 16:05 and the rest from 19:29.

## 0. Verdict

- **Route 2 works for all of these observables.** The field enters correctly (μ(0) = 5e-7 eÅ, antisymmetric to
  1e-6, α independent of the field step to six digits, energy ½αE² matched to 1e-4), the header fields 1-4 are
  inert (byte-identical under four headers), the mesh, g_ewald and self-term choices are converged, and the cubic
  periodic cell is removed analytically plus by extrapolation (below).
- **Table 3 changes in every numeric cell**, by ≤ 0.2% in the charges and dipoles and by 0.03-0.14 eV in E; the
  inversion stands. The "periodic monomer, full Ewald" pair of §7.1 (−0.8939 / +0.9099) was not an Ewald solve: it is
  reproduced exactly by the route-0 saddle at an 8 Å taper cutoff (`swb`), and it should go.
- **The atomic-charge alkane series changes a lot**: α(C8) 116.21 → **190.83 Å³**, α(C24) 1169.11 → **1636.71 Å³**,
  exponent n^1.888 → **n^1.833** (long half 2.121 → 2.022), local p(C8→C17) 2.085 → **1.898**. The route-0 numbers
  were low because the legacy gas path multiplies the kernel by the ReaxFF 7th-order taper from r = 0 (0.5 at half
  the 12 Å cutoff) and because the saddle's intra X-block (κ_bond = 1, r_ov = 2 Å) limits the response. An
  independent dense NumPy solve of the same functional (no cutoff, no taper, no images) agrees with route 2 to
  0.00-0.27% at every length (C8 190.8317 exactly).
- **Bondsoft does not change** (≤ 7e-7 relative at every n; κ→∞ rows identical to the published ones): the bondsoft
  keyword builds its own taper-cutoff saddle and forces MINRES on every route (`fix_qeq_sam_electrode.cpp:1101-1102`).
  The SI scan rows stand; only the plain reference they are compared against moves.
- **33 numbers in tables and text change**, plus the nine Figure 2 data points and three panel labels (Section 8).

## 1. Baselines: the published route-0 decks on the current binary

Run first, so the OLD side of every comparison is measured, not transcribed (`baseline_route0/`).

| Deck (as published) | Route banner | Result | Paper |
|---|---|---|---|
| `kernel_comparison/in.gas_slater` | saddle-BiCGStab 6 mv, resid/b 1.4e-1 | q_M −0.88870, q_H +0.44435, μ 1.86112 D, pe −4.396359 | Table 3 row 1 |
| `kernel_comparison/in.gas_cbrt` | saddle-BiCGStab 7 mv, resid/b 2.2e-2 | q_M +0.91109, q_H −0.45555, μ 1.90801 D, pe +0.859385 | Table 3 row 2 |
| `kernel_comparison/in.lr_slater` (pppm in the deck) | **saddle-BiCGStab** 6 mv | q_M −0.893898, μ 1.871999 | §7.1 "−0.8939" |
| `kernel_comparison/in.lr_cbrt` | saddle-BiCGStab | q_M +0.909866, μ 1.905439 | §7.1 "+0.9099" |
| `in.gas_slater` with `swb` 12 → 8 Å (no kspace) | saddle | **q_M −0.89390**, μ 1.87200 | = the "full Ewald" value |
| `in.gas_cbrt` with `swb` 12 → 8 Å | saddle | **q_M +0.90987**, μ 1.90544 | = the "full Ewald" value |
| `polarizability` C8, E = ±0.004 | saddle-BiCGStab 18 mv | μ_x(+0.004) = 0.0322799292 eÅ, α = 116.208 Å³ | 116.21 |

The `in.lr_*` decks carry the four-field header `0.30 50.0 2.0 1.6` (lr_ewald absent → 0), so the solve took the
real-space saddle with `swb = 8` and only the forces saw PPPM; the deviation from Table 3 is the shorter taper cutoff
(the taper polynomial at r = 0.96 Å is 0.9943 for swb = 8 against 0.9986 for swb = 12), not the electrostatic path.
Note also the logged relative residuals of the route-0 kernel solves (2e-2 and 1.4e-1 at a requested 1e-9).

## 2. Route-2 protocol

Deck templates are in `run_r2.py` (`ALK_DECK`, `WAT_DECK`); every run is `runs/<phase>/in.<tag>` with its log and a
`res.<tag>` line. Parameter files: `kernel/water4_r2.param` (header `0.30 50.0 2.0 1.6 0.0 2`, per-type lines as
published) and `alkane/alkane_r2.param` (`1.0 1.0 2.0 2.5 0.0 2`, identical to the repository's `alkane_lr.param`).

- Single molecule, centred (`change_box` on the centre of mass for the alkanes, the monomer at the origin) in a
  cubic periodic cell of side L; `boundary p p p`.
- `pair_style hybrid/overlay lj/cut 12.0 coul/long 12.0 coul/shield/intra 12.0 all [slater 2s 1 3 | cbrt]`, with the
  per-type shielding coefficients of the parameter file on the pair side (3.0803/1.7007/3.0803 for the monomer,
  0.2499 for C and H), so the pair energy mirrors the solve; `kspace_style pppm/samqeq 1.0e-6` with
  `kspace_modify gewald 0.35 mesh N N N order 5`.
  **The accuracy value cannot size the grid here**: the charges are zero at setup, so q2 = 0; stock PPPM then
  requires `gewald` (pppm.cpp:991 "Must use kspace_modify gewald for uncharged system") and its grid loop breaks at
  the coarsest step because the error estimate is zero (the log prints "estimated absolute RMS force accuracy = 0").
  pppm/samqeq does not override this. The mesh is therefore set explicitly and converged (below).
- `fix chg <group> qeq/sam 1 0.0 12.0 1.0e-10 <param>`; `fix_modify chg solver cg nofallback`
  (`fix_qeq_sam_electrode.cpp:1261-1278`); `fix_modify chg recip_self peratom` (`:777-803`, exact per-atom grid
  self-coefficient, no calibration). Route banner in every log: `lr solver = CG`, `recip_self = EXACT PER-ATOM`,
  `SOLVER-DIAG: CG k iters, resid/b ≤ 9e-11`.
- Route dispatch read: `fix_qeq_sam_lr.cpp:168` (`if (!lr_ewald)` → `FixACKS2Sam::pre_force`, the saddle), `:245`
  (`if (lr_bondsoft)` → saddle with the base X-block), `:295` (`qeq_solve()`, projected CG, `:2420`).

### 2a. External field on route 2 (code)

`FixQEqBaseSam::init()` refuses a **constant** field component along a periodic axis
(`fix_qeq_base_sam.cpp:512-528`, message "Must not have electric field component in direction of periodic
boundary ... use an atom-style potential variable"); this is exactly the error in the repository's own
`polarizability/log.lr_plain`, which is why the published series stayed on the open-boundary route 0. The guard is
bypassed only for `atom_pot = (varflag == ATOM && pstyle == ATOM)`, and `FixEfield::init()` sets `varflag = ATOM`
only from an atom-style field component (`fix_efield.cpp:235-236`), not from the potential alone. The deck therefore
uses

```
variable E   equal 0.004
variable ex  atom v_E
variable phi atom -1.0*v_E*x
fix     ef all efield v_ex 0.0 0.0 potential v_phi
```

`get_chi_field()` (`fix_qeq_base_sam.cpp:1370`) then takes `chi_field[i] = qe2f*efield[i][3] = φ_i = −E x_i`
(`:1428`; metal qe2f = 1), the same quantity the constant branch forms as `−(E·x_unwrapped)` (`:1416`), and
`qeq_solve()` adds it to the right-hand side as `qb[i] -= chi_field[i]` (`fix_qeq_sam_lr.cpp:2755` OpenMP, `:2776`
serial; `:229` calls `get_chi_field()` before the solve). No force or energy from `fix efield` is used; μ comes from
`compute dipole/samqeq` (unwrapped coordinates, image flags zero).

### 2b. Efield check (C8, L = 60 Å, mesh 120, `runs/efield/`)

| Check | Result |
|---|---|
| μ_x(E = 0) | 5.3e-7 eÅ (all lengths: ≤ 1.2e-6 except C17, 9.0e-5; the C17 geometry also gives 2.75e-4 in the published bondsoft data) |
| antisymmetry μ(+E) + μ(−E) | 1.06e-6 eÅ at E = 0.001, 0.002, 0.004 V/Å |
| α by central difference ±0.001 / ±0.002 / ±0.004 | 191.5630 / 191.5630 / 191.5630 Å³ |
| α by the seven-point least-squares slope | 191.5630 Å³, R² = 1 − 3e-10 |
| energy: pe(±E) − pe(0) (pe excludes the field term, so it is +½αE²) | +2.6610e-5 eV (E = 0.002) and +1.06430e-4 (0.004) against ½αE²/14.399645 = 2.6607e-5 and 1.06427e-4 |

### 2c. Header fields inert on route 2 (C8, E = +0.004, L = 60)

`1.0 1.0 2.0 2.5 0.0 2` (as published), `1e300 0.0 2.0 2.5 0.0 2` (gate closed, no intra X), `0.30 50.0 2.0 1.6 0.0 2`
(the kernel-table header), `1.0 1.0 7.0 9.0 0.0 2`: μ_x = 0.0532137966, pe = −0.85677954 in all four, byte-identical.

### 2d. Convergence

Mesh (order 5 unless stated), the quantity is α from ±0.004:

| System | h = 1.0 Å | 0.67 Å | 0.5 Å | 0.4 Å | 1.0 Å, order 7 |
|---|---|---|---|---|---|
| C24, L = 80 | 1662.8937 | 1662.8948 | 1662.8944 | 1662.8945 | 1662.8946 |
| C8, L = 60 | 191.56203 | 191.56302 | 191.56305 | | 191.56293 |
| monomer slater, L = 40 (q_M) | −0.889864 | | −0.889819 | −0.889819 (h = 0.33 also −0.889819) | |

Chosen: h = 2/3 Å for the alkanes (1.6e-7 relative), h = 1/2 Å for the monomer (its charges sit 0.15-0.96 Å apart
and need the finer grid).

g_ewald 0.30 / 0.35 / 0.40 Å⁻¹ at a 12 Å cutoff: C8 α = 191.5630 / 191.5630 / 191.5631; monomer q_M = −0.889785 /
−0.889786 / −0.889786.

Self term: `recip_self peratom` (used) against `recip_probes 16` and the default: C8 μ_x(+0.004) = 0.0532137966 /
0.0532123332 / 0.0532123332 (2.8e-5 relative, pe identical to 1e-8). For the monomer `recip_probes 16` is **refused**
("recip_self needs >= 4 group atoms for the multi-probe calibration (found 3)", `fix_qeq_sam_lr.cpp:675`; the fq
group is H, H, M) and the legacy default single-pair calibration shifts q_M by 4e-5 (−0.889743 against −0.889786 at
L = 60), so `peratom` is the only deterministic choice there. None of these choices moves a reported digit.

**Box size.** Under conducting (tinfoil) boundary conditions, which is what PPPM applies, a cubic lattice of
identical dipoles is not field-free: each molecule sits in the Lorentz cavity field 4πμ/(3L³) of its images (the
depolarising surface term that would cancel it is what tinfoil removes). The box polarizability is therefore
α_box = α/(1 − 4πα/(3L³)), i.e. 1/α = 1/α_box + 4π/(3L³), exact at dipole order; what remains after this analytic
correction is the coupling of the molecule's higher multipoles to its images, O(L⁻⁵). The term is small only when
4πα/(3L³) is: 5e-5 for the monomer (α ≈ 0.8 Å³) at L = 40 Å, but 0.7% for C24 (α ≈ 1640 Å³) at L = 100 Å. Measured,
C24 with the 12 Å cutoff, α from ±0.004 V/Å:

| L (Å) | α_box | Lorentz-corrected | note |
|---|---|---|---|
| 40 | 2075.90 | | image pairs inside the 12 Å cutoff (31 Å molecule); discarded |
| 60 | 1708.81 | 1654.00 | h = 0.5 |
| 80 | 1662.89 | 1640.57 | |
| 100 | 1649.27 | 1637.95 | |
| 120 | 1643.7281 | 1637.2046 | h = 2/3 (the series) |
| 160 | 1639.5699 | 1636.8254 | |
| 200 | 1638.1504 | 1636.7466 | |

The corrected residual falls 5-fold per step from 60 to 100 Å, the L⁻⁵ signature; a two-cell L⁻⁵ extrapolation
of the corrected values gives 1636.7074 (120,160), 1636.7081 (160,200) and 1636.7079 (120,200): consistent to
0.001 Å³, i.e. at the paper's 0.01 Å³. The raw L = 200 value would be 0.09% high. **Protocol for the series: L = 120
and 200 Å, Lorentz correction, L⁻⁵ extrapolation.** For C8 the corrected values are 190.85402 (L = 60), 190.83256
(120), 190.83181 (200), extrapolated 190.8317; the same number, 190.8317, comes from two image-free calculations
(Section 5, Section 6), which validates the 4π/3 coefficient.

**Real-space cutoff.** Route 2 applies the shielding correction J(r) − 1/r only inside `swb` (12 Å; the T3 taper is
off by default, `pair_coul_shield_intra.cpp:30-37`) and bare 1/r from PPPM beyond it. For molecules longer than 12 Å
this is a small bias relative to the untruncated kernel: C24 at L = 200, corrected, 1636.7466 (swb 12), 1636.3164
(24), 1636.1131 (36), against the dense untruncated 1636.0745 (`runs/swbcheck/`). The published protocol also used a
12 Å cutoff, so the series below keeps it; the dense column of Section 5 shows the size of the effect (≤ 0.27%, at
C12, where the chain length first exceeds the cutoff).

## 3. Table 3 (`tab:kernelcmp`) on route 2

Monomer in cubic cells L = 40-120 Å (mesh 2L), extrapolated in L⁻³ from the two largest cells; the (80,100) pair
reproduces the (100,120) extrapolation to 2e-6 in q_M and 1e-6 D (`runs/water/`, `water_route2.json`). CG converges
in 3 iterations to resid/b < 1e-16. All runs on the 19:29 binary.

| L (Å) | slater q_M | q_H | μ (D) | E (eV) | cbrt q_M | q_H | μ (D) | E (eV) |
|---|---|---|---|---|---|---|---|---|
| 40 | −0.889819 | 0.444910 | 1.863457 | −4.539236 | 0.908808 | −0.454404 | 1.903223 | 0.824843 |
| 60 | −0.889786 | 0.444893 | 1.863387 | −4.538943 | 0.908843 | −0.454422 | 1.903296 | 0.824855 |
| 80 | −0.889778 | 0.444889 | 1.863370 | −4.538872 | 0.908852 | −0.454426 | 1.903314 | 0.824858 |
| 100 | −0.889775 | 0.444887 | 1.863364 | −4.538847 | 0.908855 | −0.454427 | 1.903321 | 0.824859 |
| 120 | −0.889773 | 0.444887 | 1.863361 | −4.538836 | 0.908856 | −0.454428 | 1.903323 | 0.824859 |
| **L → ∞** | **−0.889770** | **0.444887** | **1.863357** | **−4.538821** | **0.908857** | **−0.454429** | **1.903326** | **0.824859** |

The two kernels move in opposite directions with L: the Lorentz feedback amplifies the dipole of the definite
(Slater) solve and damps that of the inverted one, as expected for an operator that is indefinite in the cube-root
reading.

**New Table 3** (paper's rounding; deviation = (μ − 1.85)/1.85):

| Kernel | q_M | q_H | μ (D) | Deviation | E (eV) |
|---|---|---|---|---|---|
| `slater` | −0.8898 (was −0.8887) | +0.4449 (+0.4444) | 1.8634 (1.8611) | +0.7% (+0.6%) | −4.5388 (−4.3964) |
| `cbrt` | +0.9089 (+0.9111) | −0.4544 (−0.4556) | 1.9033 (1.9080) | +2.9% (+3.1%) | +0.8249 (+0.8594) |

E on route 2 is pe = ecoul + elong + f_chg, e.g. slater at L = 120: −7.966230 − 3.337313 + 6.764708 = −4.538836 eV,
with the pair side shielding the same pairs with the same per-type coefficients as the solve (the kernel
consistency guard enforces it). The published E was pe of a deck whose pair style is `lj/cut/coul/cut`, i.e. the
**unshielded** Coulomb energy of shielded-solve charges plus the on-site term; the 0.14 eV (slater) and 0.03 eV
(cbrt) shifts in E are mostly that, not the route. The caption's claim that the inverted reading has a positive
total energy still holds (+0.825 eV). The §7.1 "periodic monomer" pair is explained in Section 1: it is the route-0
saddle at an 8 Å taper cutoff, and on route 2 the Ewald monomer is Table 3 itself.

## 4. Alkane series on route 2 (`runs/series_plain/`, `alpha_series_route2.csv`, `series_route2.json`)

Seven-point sweep E_x ∈ ±0.004 V/Å at L = 120 (mesh 180) and 200 Å (mesh 300), Lorentz correction, L⁻⁵ extrapolation;
tolerance 1e-10; CG 14-25 iterations. "dense" = `dense_qeq.py`, the same functional (χ, η on the diagonal, cube-root
kernel with γ_ij = √(γ_iγ_j), global neutrality, φ = −E x) solved densely with no cutoff, taper or images.

| n | α_box L=120 | α_box L=200 | corrected 120 | corrected 200 | **α (route 2)** | α/n | dense | route 0 (paper) | new/old |
|---|---|---|---|---|---|---|---|---|---|
| 4 | 61.6307 | 61.6233 | 61.6215 | 61.6213 | **61.6213** | 15.4053 | 61.6214 | 43.4115 | 1.419 |
| 6 | 119.2078 | 119.1805 | 119.1733 | 119.1731 | **119.1731** | 19.8622 | 119.1730 | 72.1582 | 1.652 |
| 8 | 190.9209 | 190.8509 | 190.8326 | 190.8318 | **190.8317** | 23.8540 | 190.8317 | 116.2111 | 1.642 |
| 10 | 285.6258 | 285.4684 | 285.4282 | 285.4257 | **285.4255** | 28.5425 | 285.3302 | 188.4246 | 1.515 |
| 12 | 403.3126 | 402.9968 | 402.9187 | 402.9117 | **402.9111** | 33.5759 | 401.8287 | 270.2019 | 1.491 |
| 14 | 542.7117 | 542.1353 | 541.9986 | 541.9815 | **541.9800** | 38.7129 | 541.1073 | 366.8347 | 1.477 |
| 17 | 799.3437 | 798.0766 | 797.7978 | 797.7433 | **797.7387** | 46.9258 | 796.9344 | 559.6294 | 1.425 |
| 20 | 1117.0635 | 1114.5491 | 1114.0468 | 1113.8991 | **1113.8866** | 55.6943 | 1113.0368 | 789.5027 | 1.411 |
| 24 | 1643.7281 | 1638.1504 | 1637.2046 | 1636.7466 | **1636.7079** | 68.1962 | 1636.0745 | 1169.1120 | 1.400 |

Fits (log-log least squares, as `run_alpha.py`/`extract_alpha.py`):

| quantity | route 2 (12 Å cutoff) | dense, no cutoff | paper (route 0) |
|---|---|---|---|
| exponent, full series n = 4-24 | **1.833** (1.8333) | 1.833 (1.8327) | 1.888 |
| exponent, long-chain half n ≥ 12 | **2.022** (2.0221) | 2.025 (2.0253) | 2.121 |
| α/n rise C4 → C24 | **4.43-fold** | 4.43 | 4.5-fold |
| α(C8) | **190.83** | 190.83 | 116.21 (7-pt), 116.23 (±0.002) |
| α(C17) | **797.74** | 796.93 | 559.63 |
| local p, C8 → C17 (±0.002 central difference, `scan_route2.json`) | **1.898** (1.8976; the 7-pt sweep gives the same) | 1.896 | 2.085 |
| max |μ_x(0)| | 9.0e-5 eÅ (C17); ≤ 1.2e-6 elsewhere | | "3e-4" |
| min R² of the field sweep | 1 − 1e-10 | | "≥ 0.9999999" |

The exponent is insensitive to the residual cell dependence (1.8334 / 1.8333 / 1.8333 from the corrected L = 120,
L = 200 and extrapolated columns) and to the 12 Å truncation (1.8333 against 1.8327 dense).

## 5. Bondsoft on route 2 (`runs/series_bondsoft/`, `runs/scan/`)

`fix_modify chg bondsoft 10.0 3.0` sets `lr_nrecip = 1` (Gaussian taper-cutoff H, no reciprocal term) and
`acks2_use_minres = 1` (`fix_qeq_sam_electrode.cpp:1098-1102`); `pre_force` then runs `FixACKS2Sam::pre_force` for
`lr_ewald ≠ 0` at `fix_qeq_sam_lr.cpp:245-247` and for `lr_ewald = 0` at `:168-180`. Confirmed in the logs
("bond-softness ACKS2 ON ... + Gaussian taper-cutoff H", "saddle-MINRES 64 matvecs") and by measurement: the
bondsoft C24 value is 15.044358 at L = 120 and 15.044373 at L = 200 (no cell dependence, as with no reciprocal term).

| n | route 2 | published route 0 | relative |
|---|---|---|---|
| 4 | 2.463712 | 2.463712 | −1.1e-7 |
| 6 | 3.679450 | 3.679449 | +8e-8 |
| 8 | 4.924004 | 4.924004 | +4e-8 |
| 10 | 6.181223 | 6.181223 | −3e-8 |
| 12 | 7.440149 | 7.440147 | +2e-7 |
| 14 | 8.701239 | 8.701238 | +1e-7 |
| 17 | 10.669406 | 10.669406 | −5e-8 |
| 20 | 12.498737 | 12.498737 | +1e-8 |
| 24 | 15.044358 | 15.044368 | −7e-7 |

Exponents 1.013 (full) and 1.016 (n ≥ 12), α/n flat to 1.8% from C4 to C24 (endpoints, the paper's figure; 2.3%
max/min), local p 1.026: **unchanged**. The SI scan rows were therefore not rerun except the bold row and the κ → ∞
limit at cutoff 3.0 (±0.002, C8 and C17):

| κ | α(C8) | α(C17) | p | published |
|---|---|---|---|---|
| 10 | 4.9240 | 10.6694 | 1.026 | 4.924 / 1.026 |
| 10⁴ | 125.8134 | 642.7897 | 2.164 | 125.81 / 2.164 |
| 10⁵ | 129.0490 | 682.4051 | 2.209 | 129.05 / 2.209 |
| 10⁶ | 129.3819 | 686.6407 | 2.214 | 129.38 / 2.214 |

The limit is unchanged, α(C8) → 129.4 Å³, but the comparison changes sign: it was "11% above" the route-0 atomic
charge solve (116.23) and is **32% below** the route-2 one (190.83), with p → 2.21 against 1.90. The stated reason
survives, since the bondsoft route solves its own taper-cutoff real-space operator, without the reciprocal-space
term, and its κ → ∞ limit is not the plain solve.

## 6. Why route 0 was low: two contributions, both measured (`baseline_route0/alkane/`)

The legacy gas path multiplies the kernel by the ReaxFF 7th-order taper Tap(r) from swa = 0 to swb
(`fix_qeq_base_sam.cpp:573-600` builds it, `calculate_H` `:962-978` applies it; the Slater override does the same,
`fix_qeq_sam.cpp:1510-1518`). With swb = 12 Å, Tap = 0.5 at r = 6 Å. Pushing swb out on route 0 (open box, κ_bond = 1):
C8 α = 116.208 (12 Å), 149.062 (24), 163.362 (48), 165.132 (96), 165.267 (192). The rest is the saddle's intra
X-block, whose κ_bond → ∞ limit is the plain solve: at swb = 192, κ_bond = 1, 10², 10⁴, 10⁶, 10⁸ give C8 165.27,
190.52, 190.81, 190.82, 190.82, and at swb = 2000 Å, κ = 10⁸: **190.8317**, the route-2 number to four decimals.
For C24 the same limit gives 1424.73, 13% below both route 2 (1636.71) and the dense solve (1636.07); the saddle's
κ → ∞ limit is not the plain solve for the long chains, which is a property of that formulation and is left to the
author (BiCGStab residuals 1e-6 to 1e-5 there).

## 7. Figure 2 candidate

`alpha_series_route2.csv` (new; columns `mode,n,npts,alpha,alpha_per_n,mu0,r2` as `figures/alpha_series.csv`, plus
the per-cell columns) and `fig_alpha_route2.py`, a copy of `figures/fig_polarizability.py` reading it, the
unchanged `figures/scan_grid.csv` and the new plain reference from `scan_route2.json`. Rendered:
`fig_alpha_route2.pdf` / `.png`. Panel (a) label n^1.833, panel (b) "×4.4, no saturation", panel (c) plain-solve
marker at α/n = 23.85, p = 1.898 (was 14.53, 2.085). The original files were not touched.

## 8. Exposure list, OLD → NEW, with proposed wording

Line numbers refer to `main.tex` at 17:20 and to labels in `si.tex` (its line numbers moved during this session).

**Main §7.1 (`sec:port`)**
1. l.832 "a gas-phase dipole of 1.8611 D" → **1.8634 D**.
   > Using the exact Slater kernel with the original Slater exponents of Rick, Stuart and Berne, in the intramolecular-shielding mode that matches the structure of their model, and adjusting nothing, we obtain a gas-phase dipole of 1.8634 D against their 1.85 D and a liquid dipole of 2.6404 D against their 2.637 D.
2. l.853-855 "The inversion persists under full Ewald: on the periodic monomer the exact-Slater kernel returns
   q_M = −0.8939 against +0.9099 for the cube-root reading, so it is a property of the kernel rather than of the
   short-range path." → the numbers are a real-space saddle solve at an 8 Å cutoff (Section 1) and Table 3 is now
   itself the Ewald solve. **Delete**, or replace by
   > Table 3 is a full-Ewald solve of the isolated monomer; the inversion is a property of the kernel and not of the electrostatic path.
3. l.857-861 J(0)/η = 2.15, 2.75, 1.60 and the eigenvalues −2.18 / +5.98 eV: parameter and rigid-molecule operator
   statements, **unchanged** (the image term perturbs the operator at the 1e-5 level). Flag only.
4. Table 3 (l.866-886): ten numeric cells, Section 3. Caption: "Both decks are distributed with the code, with this
   output as their reference" requires route-2 versions of `in.gas_slater`/`in.gas_cbrt` (and of the `tip4pfq_gas`
   regression case) if the table is to remain reproducible from the shipped decks. Proposed caption sentence:
   > Each run is one molecule in a cubic periodic cell solved with the reciprocal-space term inside the matrix-vector product and extrapolated to infinite dilution (Supplementary Material, Section S9); E is the reported potential energy, which includes the on-site term of Eq. (functional), the shielded real-space pair energy and the reciprocal-space energy.

**Main §7.5 (`sec:alpha`)**
5. l.1078-1079 "on the isolated molecule" → protocol changed:
   > obtained α∥ = dμ_x/dE_x by symmetric finite field along the chain axis on a single molecule in a cubic periodic cell, with the reciprocal-space term inside the matrix-vector product and the cell dependence removed as described in the Supplementary Material.
6. l.1081-1082 "vanishes to 3×10⁻⁴ eÅ" and "R² ≥ 0.9999999": both still hold (9e-5 and 2.8e-4 max; R² ≥ 1 − 1e-10);
   may stay or tighten to "1×10⁻⁴ eÅ" only if the bondsoft C17 point (2.8e-4) is also requoted. Flag.
7. l.1087-1089 "n^1.888 ... n^2.121 ... rising 4.5-fold" → **n^1.833, n^2.022, 4.4-fold**.
   > α∥ ~ n^1.833 over the series and n^2.022 over the long-chain half, with α∥/n rising 4.4-fold and no sign of saturation.
8. l.1093-1095 bondsoft n^1.013, n^1.016, 1.8%: **unchanged**. l.1095 "Two scaling exponents, 1.89 and 1.01" → **1.83 and 1.01**.
9. l.1103-1104 "0.62 Å³ per carbon", "1.27, then 1.86, then 2.16": bondsoft values, **unchanged**.
10. Figure 2 (`figures/fig_polarizability.pdf`, data `figures/alpha_series.csv` plain rows): Section 7.

**Main §6.1 (`sec:definite`)**: the rigid single-molecule eigenvalues 2.26/7.79 eV (native water) and the cbrt-reading
indefiniteness are operator statements, not route results: **unaffected**.

**Main §8.2 limitations** l.1153-1157 "about threefold low": bondsoft magnitude against experiment, **unchanged**.

**SI S9 (`sec:si-alpha`)**
11. Protocol paragraph "We treated each molecule as isolated, with open boundaries and no dispersion, and converged
    the charge solve to 10⁻¹⁰." →
    > We treated each molecule as a single molecule in a cubic periodic cell of side 120 and 200 Å, with pppm/samqeq (g_ewald = 0.35 Å⁻¹, mesh spacing 0.67 Å, order 5), a 12 Å real-space cutoff, the exact per-atom reciprocal self-term and no dispersion, and converged the projected conjugate-gradient solve to 10⁻¹⁰. The field enters as the atom-style potential φ = −E_x x. Under conducting boundary conditions the images of a molecule in a cubic cell exert the Lorentz field 4πμ/(3L³), so the cell polarizability is α_L = α/(1 − 4πα/(3L³)); we remove this term analytically and extrapolate the L⁻⁵ remainder from the two cells. For C24 the two corrections are 0.4% and 0.002%; a third cell of 160 Å reproduces the extrapolated value to 0.001 Å³, and a dense solve of the same functional without cutoff or images agrees to within 0.3% at every chain length.
12. Table S5 (`tab:si-alpha`), atomic-charge columns: 43.41/10.85 → **61.62/15.41**; 116.21/14.53 → **190.83/23.85**;
    270.20/22.52 → **402.91/33.58**; 559.63/32.92 → **797.74/46.93**; 1169.11/48.71 → **1636.71/68.20**; scaling
    n^1.888 → **n^1.833**. Bondsoft columns unchanged.
13. Scan text "against the atomic-charge-solve reference p = 2.085 with α(C8) = 116.23 Å³, the atomic-charge solve
    evaluated with the two-point central difference at ±0.002 V Å⁻¹ that the scan uses; the seven-point sweep of
    Table S5 gives 116.21 Å³" →
    > against the atomic-charge-solve reference p = 1.898 with α(C8) = 190.83 Å³; the two-point central difference at ±0.002 V Å⁻¹ that the scan uses and the seven-point sweep of Table S5 agree to four decimals.
14. "The limit is p → 2.21 with α(C8) → 129 Å³, 11% above the atomic-charge-solve reference, because the bond-softness
    route solves with its own tapered real-space operator" →
    > The limit is p → 2.21 with α(C8) → 129 Å³, 32% below the atomic-charge solve and with a higher exponent, because the bond-softness route solves with its own tapered real-space operator and without the reciprocal-space term (fix_modify <id> bondsoft switches the solve to the taper-cutoff form).
15. "the atomic-charge solve overshoots by about 7.5-fold at C8" → **about 12-fold** (190.83/15.5 = 12.3).
16. "threefold-low plateau" and Table S6 (`tab:si-scan`): bondsoft numbers, **unchanged**.
17. Table S4 (`tab:si-metal`, gold slab): route 2 already, **unaffected**.

Count: Table 3 10, §7.1 3 (1.8611 and the pair), §7.5 4 (1.888, 2.121, 4.5, 1.89), Table S5 11, SI text 5
(2.085, 116.23, 116.21, 11%, 7.5-fold): **33 numbers**, plus the nine Figure 2 data points and three panel labels.

## 9. Files

- `run_r2.py` driver (phases smoke, protocol, efield, water, series, boxcheck, kappa); `dense_qeq.py` reference
  solve; `fig_alpha_route2.py`; `alpha_series_route2.csv`, `series_route2.json`, `scan_route2.json`,
  `water_route2.json`, `boxcheck_route2.json`, `dense_series.txt`; `phase*.out` phase summaries.
- `runs/<phase>/{in,log,res}.<tag>`; `runs/water/charges_<kernel>_L<L>.dump` (12-digit charges);
  `runs/water_bin1605/` the 16:05-binary water set; `runs/bittest/` the binary comparison; `runs/swbcheck/`;
  `baseline_route0/` the published decks and the taper / κ_bond experiments.

## 10. Reproduce

```
W=~/manuscripts_todo/ms_alphasigma_overleaf/gasphase_route2_20260916; cd $W
# published route-0 decks on the current binary (Section 1)
( cd baseline_route0/kernel && for d in in.gas_cbrt in.gas_slater in.lr_slater in.lr_cbrt in.gas_slater_swb8 in.gas_cbrt_swb8; do
    mpirun -np 1 ~/codes/lammps/lammps-dev2/build/lmp_parallel -in $d -log log.$d < /dev/null > out.$d 2>&1; done; grep -h "^GAS_\|^LRKERNEL" out.* )
# route 2: checks, kernel table, series, box check, scan rows   (NP ranks per run, WORKERS concurrent runs)
NP=2 WORKERS=2 python3 -u run_r2.py protocol      # box/mesh/g_ewald/self-term/header checks (C24, C8)
NP=2 WORKERS=2 python3 -u run_r2.py efield        # mu(0), antisymmetry, step independence, energy check
NP=2 WORKERS=2 python3 -u run_r2.py water         # Table 3, L = 40..120, -> water_route2.json
NP=2 WORKERS=2 python3 -u run_r2.py series        # 9 lengths x 7 fields x 2 cells + bondsoft -> alpha_series_route2.csv
NP=2 WORKERS=2 python3 -u run_r2.py boxcheck      # third cell (C24 L=160), C8 cells, bondsoft cell independence
NP=2 WORKERS=2 python3 -u run_r2.py kappa         # plain +-0.002 reference, bondsoft kappa = 10, 1e4, 1e5, 1e6
python3 dense_qeq.py 4 6 8 10 12 14 17 20 24      # image-free, cutoff-free reference
python3 fig_alpha_route2.py                       # fig_alpha_route2.pdf/png
# a cached run is skipped; delete runs/<phase>/res.<tag> to repeat it. Each res line ends with the binary stamp.
```
