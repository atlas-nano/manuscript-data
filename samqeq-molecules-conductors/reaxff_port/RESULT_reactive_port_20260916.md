# RESULT — zero-refit reactive port: `fix qeq/sam` reproduces `fix qeq/reaxff` (2026-09-16)

For α·σ highlight H4 ("reactive and nonreactive force fields"). Plan: `../PLAN_restructure_highlights_20260916.md` §9.

**Result: on two stock LAMMPS ReaxFF systems, samQEq reproduces the ReaxFF charge step to solver tolerance.
The atomic charges agree to ≤1.2e-9 e, the ReaxFF energy to ≤8e-9 kcal/mol, and the forces to 7e-8 kcal/mol/Å.
A 20-step NVT trajectory tracks the reference to 5.5e-8 kcal/mol. No parameter is refit. The only
reconciliation is a units constant, read from the code (below).**

## Versions

- Binary `~/codes/lammps/lammps-dev2/build/lmp_parallel` (2026-09-12), LAMMPS (30 Mar 2026). Its
  `src/SAMQEQ` is byte-identical to samQEq `deb4a08` for `fix_qeq_sam.cpp`, `fix_qeq_sam_lr.cpp`,
  `fix_qeq_sam_electrode.cpp`, `fix_qeq_base_sam.{cpp,h}` and `fix_acks2_sam.cpp`. REAXFF is stock
  (`fix_qeq_reaxff.cpp` last touched `fa4e558889`, 2026-02-18).
- Systems: `examples/reaxff/water` (3000 atoms, QEq water of Achtyl et al. 2015, `qeq_ff.water`) and
  `examples/reaxff` RDX (21 atoms, `ffield.reax`, C/H/O/N). Inputs are copied here unchanged.

## How the two operators map (read from the code, not tuned)

| | `fix qeq/reaxff` | `fix qeq/sam` (units real) |
|---|---|---|
| χ | ffield `chiEEM`, eV | param column 2, kcal/mol |
| η | **2 ×** ffield `etaEEM` (`reaxff_ffield.cpp:200`), eV | param column 3, kcal/mol |
| kernel | Taper(r)·14.4/∛(r³+(γᵢγⱼ)^(-3/2)) (`fix_qeq_reaxff.cpp:734-749`) | Taper(r)·qqrd2e/∛(r³+1/g³), g=√(γᵢγⱼ) (`fix_qeq_sam.cpp:1540-1545`, `cutoff on`) |
| Coulomb constant | **14.4 eV·Å exactly** (`fix_qeq_reaxff.cpp:52`) | `force->qqrd2e` = 332.06371 kcal/mol·Å |
| taper | 7th order, swa=0, swb=10 | same polynomial (`Tap[]`, same fix arguments) |
| constraint | s/t two-solve, global neutrality | projected CG, one global fragment (no molecule IDs) |

Charges depend only on the ratios χ:η:J. Two samQEq parameter sets were built:
- **`opmatch`**: χ and η in eV × (332.06371/14.4 = 23.0599799), which makes the operator identical to ReaxFF's.
- **`phys`**: × 23.060549, samQEq's own eV→kcal/mol factor (`fix_qeq_base_sam.cpp:108`). Its J is 14.399645/14.4,
  i.e. 2.5e-5 relative weaker in effect: ReaxFF rounds the Coulomb constant, and this is the physically
  correct conversion.

**Only the samQEq `lr_ewald=2` + `fix_modify cutoff on` route is plain QEq.** On `lr_ewald=0` the solve goes
through the inherited ACKS2 saddle with samQEq's band-alignment X-block (`FixQEqSam::compute_X`, `calc_w`):
without molecule IDs every pair counts as inter-fragment and gets a Lorentzian-gated weight, which is the
charge-transfer mode, not QEq. That route was not used.

## Numbers

Per-atom charges (`python3 compare.py <ref> <sam>...`; `compare_water.txt`, `compare_rdx.txt`, `compare_ctrl.txt`):

| system | deck | max\|Δq\| (e) | rms Δq (e) | max rel. | CG matvecs |
|---|---|---|---|---|---|
| water, 3000 | opmatch | 1.19e-09 | 1.68e-10 | 1.6e-09 | 69 |
| water, 3000 | phys | 6.66e-05 | 3.89e-05 | 8.3e-05 | 69 |
| RDX, 21 | opmatch | 5.46e-10 | 2.41e-10 | 1.5e-08 | 19 |
| RDX, 21 | phys | 9.25e-05 | 4.72e-05 | 8.3e-04 | 19 |
| water, **negative control** (η without ReaxFF's factor 2) | opmatch | **3.12** | 1.12 | 8.8 | — |

Reference charge ranges: water [−0.8109, +0.4238] e, RDX [−0.3112, +0.5265] e. Total charge ≤5e-12 e in every run.

Energies (kcal/mol; ReaxFF sub-style energy `c_reax` in the samQEq decks, where pe also carries coul/long + kspace):

| system | reference pe | opmatch `c_reax` | phys `c_reax` |
|---|---|---|---|
| water | −257042.978119953 | −257042.978119945 | −257042.946057478 |
| RDX | −1884.30808269515 | −1884.30808269516 | −1884.30795039755 |

**Drop-in deck** (`water_dropin/`): `fix_modify chg energy no` (pair reaxff already carries the polarization
energy), `pair_modify pair coul/long compute no`, `kspace_modify compute no`. Then the **total** pe is
−257042.978119945 against −257042.978119953, **forces agree to 7.27e-8 kcal/mol/Å** (max |F| 81.9), and the
charges are as above.

**Dynamics** (`md_ref/`, `md_sam/`): the example's 300 K NVT, same velocity seed, 20 steps at 0.5 fs:

| step | ref T (K) | sam T (K) | ref pe | sam pe |
|---|---|---|---|---|
| 0 | 300.0000000000 | 300.0000000000 | −257042.9781199530 | −257042.9781199448 |
| 10 | 301.2919028691 | 301.2919028698 | −257054.9326198313 | −257054.9326198930 |
| 20 | 297.9061081021 | 297.9061080941 | −257024.6186365511 | −257024.6186366064 |

## Commands

```
cd reactive_port
(cd water_ref && mpirun -np 1 $LMP -in in.ref -log log.lammps -screen none)
(cd water_sam_opmatch && mpirun -np 1 $LMP -in in.sam ...)      # likewise water_sam_phys, rdx_*, water_ctrl_noeta2,
                                                                # water_ref_forces, water_dropin, md_ref (in.md), md_sam (in.md)
python3 compare.py water_ref/charges.dump water_sam_opmatch/charges.dump water_sam_phys/charges.dump
python3 compare.py rdx_ref/charges.dump rdx_sam_opmatch/charges.dump rdx_sam_phys/charges.dump
```
The parameter files are regenerated by `python3 make_params.py` (into `regen/`, which reproduces the four used files value for value); the negative-control file halves η by hand.
`sam_*.param` / `rdx_sam_*.param` / `sam_ctrl_noeta2.param` are the files actually used; `rdx_eem.json` holds
the extracted C/H/O/N EEM values.

## Caveats (state these if the paper uses the result)

1. **Units constant.** Agreement at solver tolerance needs χ and η scaled by 332.06371/14.4, because ReaxFF's
   QEq rounds the Coulomb constant to 14.4 eV·Å. With samQEq's physical conversion the charges differ by
   ≤9e-5 e and the energy by 0.03 kcal/mol (water). That is a property of ReaxFF, not a defect of either code.
2. **The composition needs scaffolding.** The plain-QEq samQEq route (`lr_ewald=2`, `cutoff on`) reads
   g_ewald from a kspace style, so the deck carries `coul/long` + `pppm/samqeq` (switched off with
   `compute no`), and needs `kspace_modify gewald` when the data file is uncharged (RDX: "Must use
   kspace_modify gewald for uncharged system"). The on-site energy must be excluded (`energy no`), or pe
   double-counts ReaxFF's own polarization term. A user who forgets any of the three gets wrong pe/forces
   with correct charges.
3. **No guard covers this composition.** `check_shield_consistency()` compares the fix against
   `coul/shield/intra` only and returns early when that style is absent (`fix_qeq_sam.cpp:898-899`). Nothing
   checks that the fix's kernel (e.g. `shield pqeq`) matches pair reaxff's cube-root Coulomb, so a mismatched
   deck would run silently. This is the paper's second silent-substitution mechanism, in a new composition.
4. Single configuration per system, one rank, global neutrality only. ReaxFF's other charge fixes
   (`acks2/reaxff`, `qtpie/reaxff`) were not tested. The `lr_ewald=0` route is the band-alignment mode and
   does **not** reproduce QEq.
5. Not a performance comparison. The loop times (11.4 s vs 10.0 s for 20 steps) come from one run each at
   tolerance 1e-10.
