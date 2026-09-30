# RESULT: clean NVE for the charge solve, and the grid self-term sensitivity (referee items 4 and P10)

2026-09-28. Binaries: `lmp_snapshot_nve` (a copy of lammps-release build `757d00a8`, the item-2 build; charge code
identical to the release for these paths). Local workstation, one rank per run.

## Item 4: NVE, `run_nve.py` → `runs/`, harvest `python3 analyze_nve.py` → `nve.csv`

The system is the 64-molecule SPC-FQ benchmark, rigid molecules, at 1.00 g/cc (sc 3.1 lattice box).
- Protocol: 100 ps NVT at 298 K, then 200 ps NVE.
- "Matched" deck: coul/shield/intra cutoff 10 Å = solve cutoff, `pair_modify shift yes` on lj/cut.
- Drift: slope of etotal, in K/ns over 381 rigid-body dof, ± the standard error of 10 segment slopes.
- sd(E): scatter of etotal about its linear trend.

| arm | matvec/step | drift (K/ns) | sd(E) (meV) | note |
|---|---|---|---|---|
| matched, tol 1e-5 | 12.94 | +7.7 ± 12.0 | 19.1 | |
| matched, tol 1e-7 | 17.99 | +11.0 ± 13.8 | 21.8 | **table row "full solve"** |
| matched, tol 1e-9 | 22.95 | +18.3 ± 21.1 | 22.4 | |
| matched, tol 1e-9, dt 0.5 fs | 22.87 | +41.7 ± 19.4 | 22.1 | |
| matched, warm start, tol 1e-7 | 13.59 | +19.6 ± 18.6 | 22.0 | **183 ps only**: the process died at NVE step 183 300 when the session was interrupted |
| matched, aspc 2, tol 1e-7 | 2.00 | −1181.8 ± 28.9 | 27.6 | <T> 244 K; μ 2.89 → 3.10 D; accepted resid/b median 4.5e-4 (acceptance rtol 1e-2) |
| matched, xl 5e-5 0.5 | — | (freezes) | 1520 | T 320 → 69 K in 10 ps, → 1.8 K at 110 ps; μ → 3.42 D |
| **mismatched** (shield 5 Å, no shift), tol 1e-5 | 12.97 | **+80.8 ± 28.2** | **70.7** | reproduces the old paper row (+82) |

Timing (`timing/`, isolated, 2000 NVT + 5000 NVE steps, ms/step): full 15.02 / 14.28 / 15.90; warm 16.53 / 12.55 / 14.84
(noise ±1.5 ms, so the quarter fewer matvecs is not resolved); aspc 9.21; xl 8.88.

Conclusions: with the cutoffs matched, the full solve conserves energy within the resolution of 200 ps, independent of
tolerance. The old table's +82 K/ns was the cutoff mismatch. The predictor-corrector and XL do not sustain NVE at the
documented settings. Paper table `tab:propagators` and its text are rewritten; the README aspc/xl claims are corrected.

## P10: grid self-term vs dipole, `run_selfterm.py` → `runs_selfterm/`, harvest `python3 run_selfterm.py --analyze`

Method: one 70 ps trajectory of the benchmark as shipped (shield 5 Å, pin 0.37933). The 51 frames (20–70 ps, every 1 ps)
were re-solved (`read_dump` + `reset_timestep 0` + `run 0`) at fixed `kspace_modify gewald 0.33957131 mesh 18 18 18`.

| self-term | value (raw) | <μ> (D) | vs pin |
|---|---|---|---|
| pinned | 0.37933 | 2.5592 | — |
| probe calibration (K=16) | 0.38207 | 2.6225 | +2.47 % (frame spread 0.0015 D) |
| exact per-atom | (per atom) | 2.6488 | +3.50 % (0.0021 D) |

The converged-mesh limit is 2 g/√π = 0.38317. The pin sits 1.0 % below it and the probe 0.3 % below.

**Ledger, found while measuring:**
1. **PPPM re-tunes g_ewald at every `run` from the current charges.** Here 0.3396 with the initial charges and 0.3452
   with the solved ones, and the mesh goes from 18³ to 20³. A pinned recip_self is therefore stale in every run after
   the first unless `kspace_modify gewald/mesh` is fixed. The benchmark trajectory's own production run carried a pin
   2.6 % off its split, which is the P10 hazard inside the shipped deck. The new warning (`fix_qeq_sam_lr.cpp`, checked
   once per g_ewald) catches it.
2. My first two re-solve attempts were invalid: `rerun` does not trigger the charge solve, and `run 0` after `read_dump`
   skips it because `setup_pre_force` re-solves only from step 0. Both gave identical dipoles in all arms. Fixed with
   `reset_timestep 0`.
