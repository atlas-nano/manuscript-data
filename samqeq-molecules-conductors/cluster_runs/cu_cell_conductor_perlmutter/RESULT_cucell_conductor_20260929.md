# RESULT: the 13,216-atom Cu/ZnCl2(aq) cell in the conductor mode (referee M5), 2026-09-29

**Why:** in the published benchmark deck every Cu atom (and every Zn2+) has its own molecule ID, so under per-molecule
constraints each is a single-atom fragment pinned at q0 (Cu 0, per `cucell.param` "INERT q0=0"; Zn +2). The electrode
was electrostatically inert and the paper's "solves the charges of ... the copper" was false.

**Conductor deck** (`in.cucell_conductor`, `data.cell_L0_conductor.data`): the lower and upper electrodes are molecule
90001 and 90002 (800 Cu each); `fix_modify chg gself on width 0.5 types 7`; the eigenvalue ridge is removed. Everything
else is as published (ASPC 2, quartic wall on O, SHAKE, 0.25 fs, 500 steps, metal frozen). Inert = the published deck.
The only other change to both decks is the kernel keyword pqeq → gaussian, which the current binary requires (same kernel).
Perlmutter debug job **59102811**, `lmp_cpu` md5 dbf62bae (dev build, 09-28; the deck needs pair rexpon, which is not in the release).
The first two submissions (59102032, 59102336) failed at parse: a missing `-var GRID none`, then the pqeq rename.

| arm | np | ts/s | Pair % | Modify (solve) % | ASPC acc/rej |
|---|---|---|---|---|---|
| conductor | 16 | 8.738 | 75.20 | 18.08 | 395/0 |
| conductor | 32 | 12.612 | 72.16 | 19.92 | 395/0 |
| conductor | 64 | 16.429 | 65.05 | 26.00 | 395/0 |
| conductor | 128 | 18.984 | 61.06 | 28.89 | 390/1 |
| inert | 16 | 8.591 | 74.49 | 18.90 | 395/0 |
| inert | 32 | 12.319 | 71.04 | 20.99 | 395/0 |
| inert | 64 | 16.780 | 66.78 | 24.30 | 395/0 |
| inert | 128 | 19.502 | 62.30 | 27.53 | 395/0 |

- The inert baseline reproduces the published figures (paper 8.0/12.2/16.9; published logs 8.69/12.32/16.75).
- The conductor arm is within 4 % of it.
- Conductor efficiency: 72 % and 65 % per doubling; 47 % from 16 to 64.
- Step-0 full solve: 24 CG iterations (conductor) against 8 (inert).
- Electrodes at step 500 (np64, `metal_q.dump`): each neutral (sums −2.1e-8, +4.0e-7 e); max|q| 0.090 / 0.072 e; sd 0.010 e.
- Lower electrode layers, from the back face to the electrolyte face (z 26 → 31): −0.158, −0.002, −0.031, **+0.191 e** (face,
  lateral sd 0.020 e). Upper electrode (z 127 → 132): **−0.054 e** (face, lateral sd 0.022), −0.109, +0.004, +0.158.
- The warnings are common to both arms: the quartic Picard iteration is not converged in 8 passes on O; the ASPC capped
  corrector stays above tolerance; the ionfield notice. The ridge Lanczos warning appears only in the inert arm.
- The paper §9 is rewritten with these numbers; the Limitations item is updated.

## Mobile electrodes (Tod: "rerun with a Cu EAM/FS or MEAM potential, so that the electrode need not be fixed")

Deck `in.cucell_mobile`, which differs from `in.cucell_conductor` in four ways:
- Cu–Cu is **Mendelev Cu2 (2012) EAM/FS** (`~/ff/EAM/Mendelev_Cu2_2012.eam.fs.real`; Mendelev & King, Phil. Mag. 93, 1268 (2013),
  doi 10.1080/14786435.2012.747012, Crossref-verified).
- The copper is integrated and thermostatted with the electrolyte (the setforce was removed; metal added to `mobile`).
- Zn2+ leaves the metal potential. Zn2+–Cu uses an illustrative LJ equal to Li+–Cu (0.10 kcal/mol, 1.80 Å).
- MEAM was not used: `~/ff/EAM` has Cu and Zn only as single-element files, with no Cu–Zn cross terms.
- Mendelev a0 = 3.639 Å (relaxed locally) against the electrode's built a = 3.614 Å: 0.7 % lateral compression in the fixed box.

Perlmutter interactive job **59104218** (nid200012), same `lmp_cpu` md5 dbf62bae; 5-step smoke passed.

| np | ts/s | Pair % | solve (Modify) % | ASPC acc/rej |
|---|---|---|---|---|
| 16 | 8.515 | 72.80 | 20.76 | 390/1 |
| 32 | 12.488 | 71.92 | 19.73 | 395/0 |
| 64 | 16.378 | 64.84 | 26.34 | 390/1 |
| 128 | 20.127 | 64.45 | 24.85 | 395/0 |

Efficiency 73 % / 66 % per doubling, 48 % 16→64; within 3 % of the inert fixed-copper baseline (−0.9, +1.4, −2.4 %).

5 ps stability run (`runs_mobile/stab_np64`, 20000 steps, 17.8 ts/s):
- T_water about 300 K; T_Cu 290–307 K.
- Cu MSD grows to 0.136 Å², vibrational: the electrodes stay crystalline.
- Electrode sums |q| < 6e-16 e throughout.
- max|q_Cu| rises from 0.10 to about 0.2–0.25 e.
- ASPC acc/rej 19960/7; no ridge; no near-critical warning.
- Final profile: the induced charge is in the electrolyte-facing layers (lower-electrode face z ≈ 31–32: −0.076 / +0.158 e;
  upper face z ≈ 126–127: +0.240 / −0.267 e, lateral sd 0.02–0.03 e); the back layers are within 0.01 e.
- The only warnings are pre-existing (quartic Picard on O at steps 0–4, ASPC capped corrector, ionfield notice).

The paper's §9 is rewritten for the mobile-electrode run (backup main.tex.bak_premobile; refs.bib + Mendelev2013).
The Limitations and the §8 cell description are updated.
