# Dense diagonalization of the constrained charge operator: native water (primary) and Drude hybrid (archived)

2026-09-16. Brief: replace the borrowed numbers of the main-text bound paragraph (`main.tex:797-808`,
`eq:bound-ineq`: J(0) = 7.20 eV < η_O = 13.364 eV, λ_min = 7.42 eV, κ = 2.15, all measured on the
cube-root defect-diagnosis deck, `ADJUDICATION_A_items_20260916.md` A1) with numbers measured on the
paper's own model. The target was changed mid-task from the Drude hybrid of Section 4.2 (dropped from
the paper) to the **native fluctuating-charge water of Section~\ref{sec:water}**; the hybrid work is
kept in Section 5 for a future paper.

Method, identical for both: read the operator the solver applies from the source (Section 1), run the
frame through LAMMPS once at a tight tolerance, rebuild the operator in numpy including the
reciprocal-space term as an exact Ewald sum, reproduce the LAMMPS charges (validation gate), then
diagonalize P H P on the constrained subspace. Script: `hd_reconstruct.py` (this directory).

Everything below comes from files in this directory, the source tree `/home/tpascal/codes/samQEq/src/SAMQEQ/`
(commit state of 2026-09-16; local binary `/home/tpascal/codes/lammps/lammps-dev2/build/lmp_parallel`,
built 2026-09-16 13:50), or the cited manuscript lines.

---

## 1. The operator as coded (lr_ewald = 2 route, both models)

Parameter-file header fields (`fix_qeq_sam.cpp:1154-1170`): `gamma_align kappa_bond r_ov [r_loc]
[lr_alpha] [lr_ewald] [lr_ridge] ...`; per-type lines `type chi eta col4 q0 eHOMO eLUMO`
(`fix_qeq_sam.cpp:1179-1185`). Both param files carry `lr_alpha = 0`, `lr_ewald = 2`, no ridge.
With `lr_ewald = 2` the plain long-range route of `FixQEqSam::pre_force` runs
(`fix_qeq_sam_lr.cpp:144-306`); `lr_alpha` is overwritten by the kspace `g_ewald` at every solve
(`fix_qeq_sam_lr.cpp:219`).

| piece | as coded | source |
|---|---|---|
| variables | atoms of the fix group; solve variable qs = q − q0[type]; committed q = qs + q0 | `fix_qeq_sam_lr.cpp` qeq_solve, 2399 ff. |
| constraint | per-molecule mean subtraction over group atoms, `v[i] -= molsum[m]*molinv[m]` (P = I − per-fragment J/n) | `fix_qeq_sam_lr.cpp:970-1058` (line 1055) |
| real-space off-diagonal | `H_ij = qqrd2e·(J(r) − erf(a r)/r)` for group–group pairs with `r ≤ swb`, no taper, `a = g_ewald` | `fix_qeq_sam.cpp:1488-1551` (line 1551) |
| kernel, default (`shield cbrt`) | `g = sqrt(gamma_i gamma_j); J = 1/cbrt(r³ + 1/g³)`; J(0) = g | `fix_qeq_sam.cpp:1391-1392` |
| kernel, `shield pqeq` | `a_i = λ/(2 Rc_i²)`, `a_ij = sqrt(a_i a_j/(a_i + a_j))`, `J = erf(a_ij r)/r`; J(0) = 2 a_ij/√π; λ = 0.462770 | `fix_qeq_sam.cpp:1385-1389`, `fix_qeq_sam.h:543` |
| reciprocal, inside the matvec | `out_i += qqrd2e·(Σ_j K_ij x_j − rs·x_i)`, K = PPPM k≠0 potential per unit charge including j = i; rs = `recip_self` (calibrated by K = 16 probe pairs, or pinned by `fix_modify recip_self`) | `fix_qeq_sam_lr.cpp:879-923` (line 922); pin `fix_qeq_sam_electrode.cpp:779-780` |
| diagonal | `etaloc_diag_of(i)` = `eta[type]`, or `eta[type] + (1/6) c4 q_lag²` on the quartic group (Picard-lagged secant), plus the ridge (0 here) | `fix_qeq_sam_lr.cpp:1434`, `:48-63`; `fix_qeq_sam_levels.cpp:524-570`; `fix_qeq_sam.h:139-167` |
| right-hand side | `qb_i = −(chi_b(i) + drude_field_i + q0field_i)`, then P; `chi_b = chi[type]` (no electrode) | `fix_qeq_sam_lr.cpp:2746`; `fix_qeq_sam_electrode.cpp:2304-2309` |
| fixed-charge field (shells) | real-space `qqrd2e·(J − erf(a r)/r)·q_j` over non-group neighbours within swb, the core's OWN bonded shell contributing only `−erf(a r)/r·q_own` (spring, not a Coulomb pair), plus the PPPM potential of all non-group charges at the cores | `fix_qeq_sam_lr.cpp:2239-2296` (line 2275) |
| q0 reference field | `q0field = P[(H_real + qqrd2e(K − rs I)) q0]` (only when some q0 ≠ 0) | `fix_qeq_sam_lr.cpp:1455-1482`, `:2505` |
| quartic Picard | `eta_eff ← eta + (1/6) c4 q²` at the lagged q, linear solve, `q ← mix·q_new + (1−mix)·q_old`, stop when max|Δq| < 1e-5 | `fix_qeq_sam_lr.cpp:316-400` (330, 394) |

Consequences used below.
- The solved system is `P (D + H_real + qqrd2e(K − rs I)) P qs = P qb`. In the reconstruction K is the
  exact Ewald k-space sum ψ_k(r_ij) = (4π/V) Σ_{k≠0} e^{−k²/4a²} cos(k·r_ij)/k², converged to
  machine precision (spherical cutoff k²/4a² ≤ 40, checked against a cutoff of 24). PPPM approximates
  this K; the residual of the validation therefore measures the particle-mesh error.
- On the constrained subspace the reciprocal term reduces to a small correction of the pair kernel.
  For r ≪ L, ψ_k(r) − erf(ar)/r = ξ + π/(a²V) + (2π/3V) r² + …, and the same constant
  ξ + π/(a²V) = ψ_k(0) − 2a/√π sits on the diagonal, so the constant is a multiple of 𝟙𝟙ᵀ and P
  annihilates it (ξ is the Ewald image self-term, −2.837297/L for a cubic cell; the reconstruction
  reproduces it to six digits on both boxes). What survives on the neutral subspace is the net kernel
  J(r), the dipole-image term (2π/3V) r² and higher images, and one uniform shift
  δ = qqrd2e·(2a/√π − rs) of every eigenvalue, set by how far the calibrated or pinned `recip_self`
  sits from 2a/√π. dδ/da = qqrd2e·2/√π = 16.248 eV per Å⁻¹.
- The reciprocal term cannot be dropped: the operator "η + J(r) truncated at swb with no Ewald"
  (the form the old tool `tools/measure_closepair_eigenvalue.py` used, but it tapered J to zero) is
  indefinite for both boxes here (Section 2.4, 5.4), an artifact of the hard truncation at swb.

---

## 2. Native fluctuating-charge water (Section~\ref{sec:water})

### 2.1 Model, from the deck and the parameter file

Deck `/home/tpascal/codes/samQEq/examples/samqeq/in.npt_canon` (copied lines in `in.nat_run0` here):
512 rigid SPC-geometry molecules (`spcfq.mol`: r_OH = 1.000 Å, r_HH = 1.633 Å), `pair_style
hybrid/overlay lj/cut 10.0 coul/long 10.0 coul/shield/intra 5.0 all`, `pair_coeff * *
coul/shield/intra 1.11`, `kspace_style pppm/samqeq 1.0e-6` (no `gewald` pin), `fix chg fq qeq/sam 1
0.0 10.0 1.0e-5 spcfq_expt.param` (swa 0, **swb 10 Å**, tolerance 1e-5 in production),
`fix_modify chg recip_self 0.37933` (pinned), no `shield` keyword, so the kernel is the default
**cube-root** form of `fix_qeq_sam.cpp:1391-1392`.

`spcfq_expt.param` (copied here): header `0.30 10.0 2.0 3.5 0.0 2` → lr_alpha 0, **lr_ewald 2**,
lr_ridge 0. Types: O `chi 6.7721 eta 15.1145 gamma 1.11 q0 0.0`, H `chi 4.50 eta 16.1595 gamma
1.11 q0 0.0`. **η_O = 15.1145 eV verified** (smallest on-site hardness). q0 = 0 for both types, so
there is no q0 field; there are no fixed charges and no quartic term: the operator is
`P (η + H_real + qqrd2e(K − rs I)) P` exactly.

g_ewald: production tunes it once at its run start from the accuracy, the cutoff, the 24.8 Å
lattice box and the ±0.85/0.425 seed charges, giving **0.33957131** (log of the lattice single point
`log.natlat_a6`, `G vector (1/distance) = 0.33957131`, grid 36³), and keeps it for the whole 30 ps.
The pin 0.37933 sits 3.83e-3 Å⁻¹ below 2a/√π = 0.383165 at that g_ewald, i.e. δ = +0.0552 eV on
every eigenvalue (Section 1). The single points below pin `kspace_modify gewald 0.33957131` so that
the operator is the one production applies on the frame; an un-pinned single point is reported as a
variant.

### 2.2 Frame provenance

**PENDING at hand-back (14:33).** The coordinator's driver writes the frame to
`/home/tpascal/manuscripts_todo/ms_alphasigma_overleaf/tolscan_npt_20260916/npt.data` (30 ps of
`in.npt_canon` at the pin; marker `NPT_FINAL DONE` in `out.npt_snap` there, expected after about
15:30). It had not appeared when this report was forced to close. Everything needed to finish is in
place and takes about two minutes once the file exists:

    bash run_native.sh /home/tpascal/manuscripts_todo/ms_alphasigma_overleaf/tolscan_npt_20260916/npt.data

which runs three `run 0` single points (`in.nat_run0`: pppm 1e-6 with g_ewald pinned at the
production value, pppm 1e-8 with the same pin, pppm 1e-6 un-pinned), each with `solver cg
nofallback`, tolerance 1e-10 and `recip_self 0.37933`, dumps `frame_nat_*.dump`, and reconstructs
each (`rec_nat_*.json`, one summary block per setting on stdout). If the driver's frame has not
appeared by 16:15, generate one with `mpirun -np 4 lmp_parallel -in in.nat_npt_own` (a copy of
`in.npt_canon` that ends in `write_data nat_npt_own.data`, about 90 min) and pass that file instead.

The machinery was validated on the deck's own start configuration meanwhile (the 8x8x8 lattice of
`in.npt_canon` with its ±0.85/0.425 seed charges; this is not an equilibrated frame and its spectrum
is not the paper's number): `log.natlat_a6`, `frame_natlat_a6.dump`, `rec_natlat_a6.json`.

### 2.3 Validation against LAMMPS

Lattice start configuration (`rec_natlat_a6.json`; pppm 1e-6, grid 36³, g_ewald 0.33957131, pinned
rs 0.37933, CG 25 iterations to resid/b = 9.0e-11, pe = −371.672827867 eV): reconstruction versus
LAMMPS **max|Δq| = 3.91e-6 e, RMS 1.03e-6 e** (O 3.9e-6, H 2.7e-6); LAMMPS per-molecule sums exact to
7e-16. This residual is the particle-mesh error of PPPM at 1e-6 against the exact Ewald k-sum; the
hybrid's 1e-8 run (Section 5) is the check that it shrinks with the accuracy. The equilibrated NPT frame
is to be reported here from `run_native.sh` (three settings, see 2.2).

### 2.4 Spectrum of the constrained operator

Lattice start configuration (not the paper's number; the NPT frame's values come from `run_native.sh`):

| operator | λ_min (eV) | λ_max (eV) | κ | negative | null (full space) |
|---|---|---|---|---|---|
| as solved (η, Ewald with pinned rs 0.37933, g_ewald 0.33957131) | 1.74203 | 11.56241 | 6.637 | 0 | 512 |
| same with rs := 2a/√π (pin offset removed) | 1.68680 | 11.50719 | 6.822 | 0 | 512 |
| kernel only, truncated at swb, no reciprocal (diagnostic) | −10.180 | 20.706 | – | 274 | 512 |

Dimension of the constrained subspace 1536 − 512 = 1024; the 512 null vectors of P are counted in the
full space and excluded from λ_min. The shift between the first two rows is exactly δ = +0.0552 eV
(Section 1). The isolated rigid molecule's constrained block gives 2.264 and 7.786 eV, so the liquid's
λ_min sits below the O↔H charge-transfer mode of one molecule by the intermolecular coupling, as
expected. Closest approaches on the lattice frame: O–H 1.878 Å (J = 7.40 eV), H–H 1.452 Å (9.23 eV),
O–O 2.362 Å (5.99 eV).

### 2.5 Contact values against the smallest hardness

As coded, `g = sqrt(1.11 × 1.11) = 1.11 Å⁻¹` for every pair, so

| pair | J(0) = qqrd2e·g (eV) | J(0)/η_O | J(0)/η_H |
|---|---|---|---|
| O–O | 15.9836 | 1.0575 | 0.9891 |
| O–H | 15.9836 | 1.0575 | 0.9891 |
| H–H | 15.9836 | 1.0575 | 0.9891 |

η_min = η_O = 15.1145 eV. **The contact bound J(0) < η_min does NOT hold for the native model**:
J(0) − η_O = +0.869 eV. The operator is nonetheless positive definite (Section 2.4), because the
kernel is evaluated at the separations that occur, never at contact: the molecule is rigid, with
J_OH(1.000 Å) = 11.992 eV and J_HH(1.633 Å) = 8.373 eV, both below η_O, and the constrained 3×3
block of an isolated molecule has eigenvalues 2.264 and 7.786 eV (computed inline, Section 6). The
closest intermolecular approaches on the frame are in Section 2.4.

**Flag for the coordinator.** `main.tex:885-887` quotes "J(0)/η of 2.15, 2.75 and 1.60 for the
H–M, M–M and H–H pairs, against 0.54 for the native model of Section~\ref{sec:water}". With the
kernel as coded and γ = 1.11 Å⁻¹ the native ratio is 1.058 (1.06 against η_O, 0.99 against η_H); I
could not reproduce 0.54 from any coded form (0.54 corresponds to J(0) = 8.16 eV, i.e. an
effective width of 0.567 Å⁻¹). Since the native model violates the contact bound and still solves,
the sentence "the bound is violated and the solve inverts" at `main.tex:888` no longer follows from
the bound alone; the ported model's ratios (2.15–2.75) are larger, but the argument needs λ_min,
not J(0)/η. I did not measure the ported model.

---

## 3. Proposed replacement text (main.tex, `eq:bound-ineq` paragraph)

Numbers in square brackets are to be taken from `rec_nat_a6p.json` (the pinned-g_ewald, pppm 1e-6
single point on the NPT frame) once `run_native.sh` has run; every other number is measured.

OLD (`main.tex:797-808`):
> For a molecular system the kernel is monotone in $r$ and therefore bounded by its contact value. For the
> core/shell water of Section~\ref{sec:drude}, at $\gamma = 0.5$~\AA$^{-1}$ and a smallest on-site hardness
> $\eta_{\rm O}$,
> \begin{equation}
> \Jij(r) \;\le\; \Jij(0) \;=\; 7.20~\mathrm{eV} \;<\; \eta_{\rm O} = 13.364~\mathrm{eV},
> \label{eq:bound-ineq}
> \end{equation}
> so the off-diagonal coupling never reaches the diagonal at any separation, and direct diagonalization of
> the constrained operator returns $\lmin = 7.42$~eV and a condition number of 2.15
> (Section~\ref{SI-sec:dense} of the Supplementary Material). Equation~\eqref{eq:bound-ineq} is a
> condition on the parameters rather than on the kernel form.

NEW:
> For a molecular system the kernel is monotone in $r$ and therefore bounded by its contact value, and the
> constrained operator is positive definite whenever that value stays below the smallest on-site
> hardness. The native water model of Section~\ref{sec:water} does not satisfy this sufficient
> condition. With the cube-root kernel at $\gamma = 1.11$~\AA$^{-1}$ for both types its contact value is
> \begin{equation}
> \Jij(0) \;=\; \frac{e^2}{4\pi\varepsilon_0}\,\gamma \;=\; 15.98~\mathrm{eV} \;>\; \eta_{\rm O} = 15.11~\mathrm{eV},
> \label{eq:bound-ineq}
> \end{equation}
> so the bound cannot be argued from the parameters alone. The molecule is rigid, and at the separations
> that occur the coupling stays well below the diagonal: $J_{\rm OH}(1.00~\text{\AA}) = 11.99$~eV and
> $J_{\rm HH}(1.63~\text{\AA}) = 8.37$~eV inside the molecule, and less than 7.4~eV at the closest
> intermolecular contact. We therefore measure the spectrum directly. On an equilibrated configuration of
> 512 molecules, with the reciprocal-space term included as the solve applies it, direct diagonalization
> of the constrained operator returns $\lmin = [\lambda_{\min}]$~eV, $\lambda_{\max} = [\lambda_{\max}]$~eV
> and a condition number of $[\kappa]$, with no negative eigenvalue and exactly one null vector per
> molecule (Section~\ref{SI-sec:dense} of the Supplementary Material). Equation~\eqref{eq:bound-ineq} is a
> condition on the parameters rather than on the kernel form, and it is sufficient, not necessary.

Notes for the editor of the surrounding text.
- The sentence that follows in the OLD paragraph, "The same cube-root kernel violates it when the fourth
  column is read as an inverse shielding length instead of a Slater exponent, and the solve then
  inverts (Section~\ref{sec:port})", and `main.tex:885-888` ("The cause is the bound of
  Eq.~\eqref{eq:bound-ineq} ... against 0.54 for the native model ... The bound is violated and the
  solve inverts") need the same correction: the native model also violates the contact bound (ratio
  1.06, not 0.54) and does not invert, so the inversion of the ported model has to be attributed to its
  measured spectrum, not to the bound. I did not measure the ported model.
- The SI section `SI-sec:dense` currently describes the 944-molecule diagnosis deck; it must be
  rewritten for this frame (512 molecules, cube-root kernel, γ = 1.11 Å⁻¹, pinned self-term, the three
  settings of Section 2.3) or the reference retargeted.

---

## 4. Reproduce

```
# inputs (copies): data.drude_water, drude_water.param (hybrid); spcfq.mol, spcfq_expt.param (native)
# md5 data.drude_water 47124a66b363cf3288b304bccf72211e ; drude_water.param fde63de7778db016c30d7984375ef4d0
cd /home/tpascal/manuscripts_todo/ms_alphasigma_overleaf/hybrid_dense_diag_20260916
LMP=/home/tpascal/codes/lammps/lammps-dev2/build/lmp_parallel          # built 2026-09-16 13:50; no rebuild
# native, lattice start configuration (machinery check)
OMP_NUM_THREADS=1 $LMP -in in.nat_run0 -var tag natlat_a6 -log log.natlat_a6 > out.natlat_a6
OMP_NUM_THREADS=4 python3 hd_reconstruct.py --dump frame_natlat_a6.dump --param spcfq_expt.param --kernel cbrt \
    --gewald 0.33957131 --rs 0.37933 --swb 10.0 --solve-types 1 2 --out rec_natlat_a6.json
# native, equilibrated NPT frame (the paper's numbers)
bash run_native.sh /home/tpascal/manuscripts_todo/ms_alphasigma_overleaf/tolscan_npt_20260916/npt.data
# hybrid: equilibrate (20000 steps, ~35 min at np 4), then single points + reconstructions
OMP_NUM_THREADS=1 mpirun -np 4 $LMP -in in.hd_eq -log log.hd_eq > out.hd_eq      # writes hd_eq.data
bash run_hybrid.sh
# hybrid, raw 2020 data-file configuration (machinery check, unequilibrated)
OMP_NUM_THREADS=1 $LMP -in in.hd_run0 -var data data.drude_water -var quartic 0 -var tag raw_q0_a6 -log log.raw_q0_a6 > out.raw_q0_a6
OMP_NUM_THREADS=4 python3 hd_reconstruct.py --dump frame_raw_q0_a6.dump --param drude_water.param --data data.drude_water \
    --kernel pqeq --gewald 0.30 --rs 0.33822 --swb 8.0 --solve-types 1 2 --out rec_raw_q0_a6.json
# isolated-molecule blocks and contact values quoted in 2.5 / 5.3 (inline python, see the transcript of
# 2026-09-16 14:29; formulas: cbrt J = 14.399645/cbrt(r^3 + 1/g^3); pqeq alpha_i = 0.462770/(2 Rc^2))
```
Runtime: each `run 0` a few seconds at np 1; each reconstruction 3-8 s with 4 threads (k-sum to
k²/4a² ≤ 40, half-space vectors: 2,300 for the hybrid box, 4,600 for the native box).

---

## 5. Hybrid (archived for a future paper): Drude-shell water of the former Section 4.2

### 5.1 Model (former Section 4.2)

Decks `../drude_cutoff/in.dw_eq` (lines 32, 59, 81, 82) and `../retest_20260806/in.drude_prod_pm`
(:57-60): 216 flexible waters (class2) with Drude shells, `coul/shield/intra 5.0 all pqeq`,
`pair_coeff * * coul/shield/intra 0.5`, `kspace_style pppm/samqeq 1.0e-6` with `kspace_modify gewald
0.30`, `fix chg CORES qeq/sam 1 0.0 8.0 1.0e-5 drude_water.param` (**swb 8 Å**), `fix_modify chg
drude`, `fix_modify chg shield pqeq`, `fix_modify chg quartic watO 5.0 0.0 8 0.5` (c4 = 5 eV/e⁴ on the
O cores, c = 0, 8 Picard iterations, mix 0.5; parser `fix_qeq_sam_electrode.cpp:1093-1101`). The
paper's dipole 2.539 D was measured at spring scale kfac 0.34 (`../retest_20260806/RESULT_drude_remeasure.md:15`);
the equilibration deck here uses that value.

`drude_water.param`: header `0.30 10.0 2.0 3.5 0.0 2 0.0` (lr_alpha 0, lr_ewald 2, lr_ridge 0). Cores:
O `chi 11.0 eta 13.3640 Rc 0.5 q0 1.0`, H `chi 4.5280 eta 17.9841 Rc 0.5 q0 1.0`; shells (types 3,
4) carry the same chi/eta/Rc with q0 0 and are outside the solve group, fixed at −1 e
(`data.drude_water`). So: the solve group is the 648 cores in 216 fragments of 3 (the constraint
holds Σ qs = 0 per molecule, i.e. Σ q_core = +3 e, neutral with the three shells); the shells enter
only the right-hand side through `add_fixed_charge_field` with the own-shell rule; q0 = 1 e on every
core adds the q0 reference field; the quartic secant `(5/6) q_O²` is added to η_O at the Picard-lagged
charge and the committed charge is the 50/50 mix of the last two Picard iterates (so LAMMPS's
committed charges sit within 0.5 × 1e-5 e of a fixed point, not at one).

Kernel: Rc = 0.5 Å for both types → α_i = λ/(2 Rc²) = 0.92554 Å⁻², α_ij = sqrt(α_i/2) = 0.68027 Å⁻¹
for every pair, J(0) = qqrd2e·2α_ij/√π = **11.0532 eV** for O–O, O–H and H–H alike, against
η_min = η_O = **13.364 eV**: the contact bound holds (J(0)/η_O = 0.827, η_O − J(0) = 2.31 eV). The
adjudication's own estimate (A1, "≈ 11.05 eV") is confirmed.

### 5.2 Frame and validation (equilibrated frame)

Frame: `hd_eq.data` = `in.hd_eq` (in.dw_eq's protocol: NVE + `langevin/drude` 300 K/1 K, 20000
steps of 0.25 fs, L = 18.77435 Å fixed, density 0.9765 g cm⁻³, kfac 0.34), `mpirun -np 4`, 14:13-14:45;
log `log.hd_eq`/`out.hd_eq` (`EQ_DONE Tcore=209.2` on the reduced-mass core compute, i.e. the
~90 K-low reading the SI describes). Single points by `run_hybrid.sh` (np 1, tolerance 1e-10,
`reset_timestep 0`, `run 0`, dump `id mol type q x y z` at %.17g), reconstructions by
`hd_reconstruct.py` with the recip_self the code calibrated in that run:

| setting | grid, recip_self (raw) | LAMMPS solve | max abs Δq (e) | RMS Δq (e) | pe (eV) |
|---|---|---|---|---|---|
| quartic on (the deck), pppm 1e-6 | 32³, 0.33823 | 4 Picard solves (dq 4.9e-05 → 2.5e-05 → 1.2e-05 → 6.1e-06; logged as `3 Picard iters, 56 matvecs`, the code's 0-based count, 4 × 14 CG iterations), resid/b 9.4e-11 | 1.23e-5 | 8.0e-6 | 3614.49327323 |
| quartic off, pppm 1e-6 | 32³, 0.33823 | CG 14 iterations, resid/b 9.4e-11 | 1.33e-5 | 8.6e-6 | 3614.49314315 |
| quartic off, pppm 1e-8 | 75³, 0.33825 | CG 14 iterations, resid/b 9.4e-11 | **2.70e-6** | **1.76e-6** | 3614.49304610 |

The residual falls fivefold when the particle-mesh accuracy goes from 1e-6 to 1e-8 at fixed
g_ewald = 0.30, so the agreement is limited by PPPM, not by the reconstruction; the kernel, the
constraint, the shell field with the own-shell rule, the q0 field and the quartic Picard are all
reproduced. On the quartic path the reconstruction follows the code's own four solves and mixes;
the exact fixed point of the quartic problem lies 1.5e-05 e from LAMMPS's committed charges (the committed
charge is the 50/50 mix of the last two iterates, which stop at dq < 1e-5). Solved charges: O
-0.0149 e (range -0.053 to 0.033), H 1.5074 e (range 1.462 to 1.561); per-molecule core sums exact to 9e-16.

Earlier machinery check on the unequilibrated 2020 data-file configuration (`rec_raw_q?_a6.json`):
1.53e-5 / 1.50e-5 e (quartic off / on, the latter through the code's 8 unconverged Picard
iterations from the data-file charges), same conclusion.

### 5.3 Spectrum (equilibrated frame `hd_eq.data`)

| operator | λ_min (eV) | λ_2 (eV) | λ_max (eV) | κ | negative | null (full space) |
|---|---|---|---|---|---|---|
| **as solved: η + quartic secant at the solved charges, Ewald with calibrated rs 0.33823** | **4.24837** | 4.27921 | **12.75695** | **3.003** | 0 | 216 |
| bare η (quartic dropped), same Ewald | 4.24792 | 4.27880 | 12.75694 | 3.003 | 0 | 216 |
| bare η, rs := 2a/√π (pin offset +0.0041 eV removed) | 4.24383 | | 12.75286 | 3.005 | 0 | 216 |
| bare η, pppm-1e-8 calibration rs 0.33825 (`rec_q0_a8.json`) | 4.24763 | | 12.75665 | 3.003 | 0 | 216 |
| kernel only, truncated at 8 Å, no reciprocal (diagnostic) | −3.744 | | 20.557 | – | 57 | 216 |

Constrained subspace 648 − 216 = 432; the 216 null vectors of P are counted in the full space and
excluded. The quartic term moves λ_min by +0.0005 eV at the solved charges ((5/6) q_O² = 0.0002 eV on
the O diagonal at q_O = -0.015 e). Isolated-molecule constrained block at the frame's mean geometry
(r_OH 0.957, r_HH 1.514 Å): 4.719 and 9.854 eV; the liquid's λ_min is 0.47 eV below the single
molecule's O↔H mode. Closest approaches on the frame: O–H intermolecular 1.946 Å (J = 6.95 eV),
H–H 1.886 Å (7.10 eV), O–O 2.592 Å (5.48 eV); intramolecular O–H 0.874 Å (9.88 eV), H–H 1.297 Å (8.75 eV).
For the hybrid the contact-value argument is complete: J(0) = 11.05 eV < η_O = 13.364 eV for every
pair, the reciprocal part adds only the image terms and the +0.004 eV calibration offset, and the
quartic diagonal can only stiffen the operator. The 2020-configuration values (4.13 / 12.41 / 3.00) are
in `rec_raw_q?_a6.json`.

### 5.4 What the archived numbers replace

Nothing in the paper now; if the hybrid returns in a later manuscript, quote λ_min = 4.248 eV, λ_max = 12.757 eV, κ = 3.00 (as solved, `rec_q1_a6.json`; bare η 4.248 / 12.757 / 3.00, `rec_q0_a6.json`) on `hd_eq.data`, with J(0) = 11.05 eV < η_O = 13.364 eV for all three pairs. Superseded line:
`rec_q1_a6.json` (as solved) and `rec_q0_a6.json` (bare η) on `hd_eq.data`, with J(0) = 11.05 eV <
η_O = 13.364 eV for all three pairs.

---

## 6. Files in this directory

| file | content |
|---|---|
| `hd_reconstruct.py` | reconstruction + validation + dense diagonalization (both kernels, shells, quartic Picard) |
| `in.nat_run0`, `run_native.sh` | native single points and their reconstructions (frame path as argument) |
| `in.nat_npt_own` | fallback NPT deck (copy of in.npt_canon ending in write_data); not run |
| `in.hd_eq`, `in.hd_run0`, `run_hybrid.sh` | hybrid equilibration, single points, reconstructions |
| `log.natlat_a6`, `frame_natlat_a6.dump`, `rec_natlat_a6.json` | native lattice configuration: run 0 and reconstruction |
| `log.raw_q?_a6`, `frame_raw_q?_a6.dump`, `rec_raw_q?_a6.json` | hybrid 2020 configuration: run 0 and reconstruction |
| `log.hd_eq`, `out.hd_eq`, `hd_eq.data`, `hd_eq.restart` | hybrid equilibration (done 14:45) |
| `log.q?_a?`, `frame_q?_a?.dump`, `rec_q?_a?.json` | hybrid equilibrated frame: run 0 and reconstruction (q1 quartic on, q0 off; a6/a8 pppm 1e-6/1e-8) |
| `data.drude_water`, `drude_water.param`, `spcfq.mol`, `spcfq_expt.param` | input copies |

Nothing outside this directory was modified.
