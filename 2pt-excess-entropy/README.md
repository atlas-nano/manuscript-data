# 2PT excess entropy of water — data and analysis scripts

Data, simulation inputs, and figure-generation scripts for:

**"Structural vs. Dynamical Estimators of the Excess Entropy of Water Across Its Phase
Diagram"** (Do, Johnson, Pascal). Preprint under the earlier title "Beyond Two-Body
Correlations: The Many-body, Excess Entropy of Water Across the Entire Phase Diagram":
DOI [10.26434/chemrxiv.15003979/v1](https://doi.org/10.26434/chemrxiv.15003979/v1).

Compares the two-phase thermodynamics (2PT) excess entropy against the two-body
approximation $S_2$ for four water models (SPC/E, TIP4P/2005, TIP4P/2005f, TIP4P/ICE)
across a 220–360 K, 0.85–1.3 g/cm³ grid.

**Archived on Zenodo:** https://doi.org/10.5281/zenodo.21824565
(all versions: https://doi.org/10.5281/zenodo.21824564)

The 2PT entropies here were produced with the reference implementation archived at
[`codes/2pt-legacy`](https://github.com/atlas-nano/codes/tree/main/2pt-legacy)
(v1.4, https://doi.org/10.5281/zenodo.7731073). Its successor,
[py-xPT](https://github.com/atlas-nano/codes/tree/main/py-xPT), implements the same
partition.

## Layout

```
2pt-excess-entropy/
  figures/
    reproduce_figures.py       regenerates main-text Figs. 5 and 6 and SI Fig. S7
    build_iapws_reference.py   rebuilds data/reference/iapws_sex_kb_corrected.csv
    fig_5_plotter_original.ipynb  the notebook used for the originally submitted Fig. 5,
                              one model and one isotherm per run, against the uncorrected
                              reference; kept as provenance, superseded by the script above
  data/
    entropy/
      2body_translational/   S_2^tr, from the O–O radial distribution function
      2body_orientational/   S_2^or, F7 factorization (Lazaridis & Karplus)
      2body_total/           S_2 = S_2^tr + S_2^or
        adjusted/            + per-model correction factor (manuscript Table 3)
      excess/                2PT excess entropy S_ex = S_liq − S_id
        adjusted/            + per-model correction factor (manuscript Table 3)
      ideal_gas/             S_id = 3·wsr_ideal + wsp_ideal (from the 2PT code output)
        rigid_rotor/         rotational ideal-gas component
        translational/       translational ideal-gas component
      liquid/                S_liq = Sq / (n · R), from the 2PT code's Sq output
      classical_weighting/   quantum vs classical weighting of the solid-like DoS,
                              TIP4P/2005 220 K isotherm (response Table R2)
    reference/
      iapws_sex_kb.csv           IAPWS-95 excess entropy as originally evaluated
      iapws_sex_kb_corrected.csv the same grid with a validity status per point (below)
    diffusivity/             2PT self-diffusivities (cm^2/s), all four models, full grid
    trajectory_counts/       windows retained per state point, by model, for the 2PT
                              excess entropy, the two-body entropy, and the diffusivity
    power_spectra/           vibrational density-of-states (.pwr) files, five state
                              points per model spanning the density and temperature range
  lammps/                    LAMMPS input decks used for the production trajectories
  2pt_input/                 2PT-code input files (rigid_2pt.in, flexible_2pt.in)
```

Each entropy CSV has columns `Density, Temperature, thermoValue` and **one row per
independent 200 ps trajectory window**, not one row per state point. Mean, standard
deviation and window count are recovered by grouping on `(Density, Temperature)`.

No raw MD trajectories are included (too large). The LAMMPS and 2PT inputs are sufficient
to regenerate them. Every figure in the paper and its Supporting Information can now be
regenerated from what is deposited here.

This directory is the authoritative copy. Earlier working copies of the same data used a
different directory layout and, in places, different column conventions; where they differ,
this one governs.

## Quickstart

```bash
pip install numpy pandas matplotlib CoolProp
cd figures
python reproduce_figures.py --ref corrected --outdir corrected   # as published in the revision
python reproduce_figures.py --ref original  --outdir control     # reproduces the original artwork
python build_iapws_reference.py                                  # rebuilds the corrected reference
```

## Figure → script + data

| Manuscript item | Script | Input data |
|---|---|---|
| Fig. 5 (S_ex, S₂ and the IAPWS-95 reference vs density) | `figures/reproduce_figures.py` | `data/entropy/{excess,2body_total}/`, `data/reference/` |
| Fig. 6 (residual heat maps) | `figures/reproduce_figures.py` | `data/entropy/*/adjusted/`, `data/reference/` |
| SI Fig. S7 (IAPWS-95 excess entropy) | `figures/reproduce_figures.py` | `data/reference/` |
| Response Table R2 (classical vs quantum weighting) | — | `data/entropy/classical_weighting/` |
| SI Fig. S8 (self-diffusivity) | — | `data/diffusivity/` |

`--ref original` exists so the script can be validated against the originally published
artwork before the corrected reference is applied. It reproduces the original worst-case
2PT residual of 1.781 k_B exactly.

## The IAPWS-95 reference, and why two versions are deposited

The reference values used in the submitted version were produced by an implementation that
silently **flashes to a two-phase mixture** when the requested density falls inside the
vapour–liquid envelope, returning the saturated-liquid entropy rather than the entropy at
the requested density. The signature is unambiguous: below saturation the values reduce to
a constant plus ln ρ, the ideal-gas term alone.

`iapws_sex_kb_corrected.csv` carries a `status` column classifying all 70 state points:

| status | n | meaning |
|---|---|---|
| `valid` | 26 | in-range single-phase liquid; the original value stands |
| `re-referenced` | 9 | inside the envelope but mechanically stable at negative pressure, the standard metastable stretched-liquid reference. Recomputed there; shifts by up to 0.385 k_B |
| `drop-spinodal` | 12 | inside the envelope **and** beyond the liquid spinodal. No homogeneous-liquid reference exists; excluded from every comparison |
| `extrap-T`, `extrap-P` | 14, 9 | outside the certified T/P limits but single-phase and stable. Retained and flagged: the limits bound the accuracy IAPWS-95 guarantees, not the existence of a solution |

Worth knowing before reusing these numbers: the largest |S₂ − IAPWS| deviation on the grid
falls from 4.058 k_B to 1.995 k_B once the unphysical points are removed. The 2PT side is
unchanged at 1.781 k_B, because its largest residual lies where the reference was always
valid.

## Calculation notes

- **Total two-body entropy**: $S_2 = S_2^{tr} + S_2^{or}$, with $S_2^{or}$ from F7
  factorization (Lazaridis & Karplus).
- **2PT excess entropy**: $S_{ex} = S_{liq} - S_{id}$, where $S_{liq} = S_q / (n \cdot R)$
  from the 2PT code's `Sq` output.
- **Ideal-gas entropy**: $S_{id} = 3\,\mathrm{wsr\_ideal} + \mathrm{wsp\_ideal}$, from the
  2PT code's printed values, i.e. three rotational degrees of freedom plus the
  translational term. The factor of three is not a typo: checked against the deposited
  components in `ideal_gas/rigid_rotor/` and `ideal_gas/translational/`, this reproduces the
  totals in `ideal_gas/` to within 0.002 k_B at all 70 state points, whereas
  $\mathrm{wsr} + \mathrm{wsp}$ is wrong by up to 3.8 k_B.
- **Adjusted values**: the `adjusted/` directories add the per-model correction factor of
  manuscript Table 3, which references the IAPWS-95 excess entropy at the single state
  point ρ = 1.0 g/cm³, T = 300 K.
- **Classical weighting**: `classical_weighting/` holds the liquid and excess entropies
  under quantum and classical weighting of the solid-like density of states. All four files
  carry their columns in the order the header states; see the directory's own `README.txt`
  for a note on the normalisation applied for this deposit.

## Licence

[CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/), as for the whole repository; the
full text is in [`LICENSE`](../LICENSE) at the repository root and is included alongside
this directory in the Zenodo archive.
