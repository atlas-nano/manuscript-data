# Weakly coordinating anions and solvents — figure data, processing scripts, inputs and force fields

Data and scripts that reproduce the figures of A. Liu, G. Raghavendran, E. S. Yang, K. J. Anderton, V. V. Singh,
A. Larson, K. Ryu, K. J. Griffith, J. S. Figueroa, Y. S. Meng and T. A. Pascal, *Weakly coordinating anions and
solvents break the conductivity–transference–solvation tradeoff*, Science (submitted). Archived at Zenodo: [10.5281/zenodo.23741204](https://doi.org/10.5281/zenodo.23741204) (version 1.0.2). The deposit also carries the
scripts that turn simulation output into the plotted quantities, one example input for each type of calculation,
and the force fields exactly as the simulations used them.

## Quickstart

```bash
python3 reproduce.py               # writes every figure to figures/out/
python3 reproduce.py --only fig1   # one figure
```

Requires Python 3 with numpy, pandas, matplotlib and Pillow (openpyxl to open `data/data_S1.xlsx`, PyYAML for `scripts/crystal_thermo.py`) (the figures were produced with Python 3.12, numpy 2.4,
pandas, matplotlib 3.10 and Pillow 12). Text files larger than 200 kB are stored gzipped (`*.gz`); every script reads
either form.

Optional environment variables:

| variable | effect |
|---|---|
| `FIGURE_FONTS` | a directory holding `arial.ttf`, `arialbd.ttf`, `ariali.ttf`, `arialbi.ttf`. Figures 1, 2, 4, S5, S9, S14, S15 and S17 are set in Arial; without it matplotlib falls back to Liberation Sans or DejaVu Sans and the glyphs differ |
| `FIGURE_PANEL_IMAGES` | a directory holding the image-only panels listed below. Without it those panels are left blank and every data panel is drawn unchanged |
| `FIGURE_OUTDIR` | where the figures are written (default `figures/out/`) |

## Layout

```
README.md            this file
reproduce.py         runs every figure script in order
figures/             one script per figure, plus shared style and reader modules
data/experimental/   measured data behind Figs. 1, 2 and 4 and figs. S5–S7 and S14–S16 (Data S1 sheets), eNMR spectra;
                     raw/ holds the instrument values behind Fig. 4 and figs. S7, S14, S16 and S18; figS17/ holds
                     the rate-test data behind fig. S17
data/simulation/     simulation results the figures read (npz, json, csv, yaml); plotted/ holds the Data S1 sheets of the
                     simulation panels; xpt/ holds the reduced 2PT/3PT output read by scripts/xpt_solvation.py;
                     vdw_structures/ holds the anion and Li⁺–solvent shell structures behind fig. S1 and tables S6 and S7;
                     esp/ holds the anion wavefunction files behind fig. S2; rdf/ holds the g(r) arrays behind fig. S20
data/tables/         tables S1–S15 (Data S1) and the tables written by the processing scripts
data/data_S1.xlsx    Data S1 as one workbook (README sheet = figure, panel, sheet, content, rows)
scripts/             processing scripts: simulation output -> plotted quantities
inputs/              one example input per calculation type (LAMMPS, PLUMED, Q-Chem, VASP, packmol)
forcefield/          the force fields as run, one merged parameter file and pair-coefficient include per system
```

## Figures

| figure | script | data |
|---|---|---|
| Fig. 1 B–D | `figures/fig1.py` | `data/experimental/Fig1B.csv`, `Fig1C.csv`, `Fig1D.csv` |
| Fig. 2 B–D | `figures/fig2.py` | `data/experimental/Fig2B.csv`, `Fig2C.csv`, `Fig2D.csv`, `Fig1C.csv`; table S2 values in `figures/figstyle.py` (= `data/tables/tableS2.csv`) |
| Fig. 3 A–E | `figures/fig3_panels_ae.py` (with `fig3_layout.py`, `fig_sources.py`) | binding energies (table S3) in `figures/fig_sources.py`; `data/simulation/snapshot_stats.json`; `data/simulation/opes_pmf/{litpfa,lipf6,libf4}_fm.npz` |
| Fig. 3 F | `figures/fig3f.py` | `data/simulation/xpt_dgmix_293K.csv`, `xpt_dgmix_anion_clamp_293K.csv`, `pairing_cycle.json`, `data/tables/crystal_dissociation.csv`, `ideal_mixing.csv` |
| Fig. 3 (plate) | `figures/fig3.py` | composes the two plates above |
| Fig. 4 A, B, D–F | `figures/fig4.py` | `data/experimental/Fig4A.csv` … `Fig4F.csv`; instrument values in `data/experimental/raw/Fig4*.csv` |
| fig. S1, tables S6 and S7 | `scripts/vdw_volumes.py` (volumes and k; the structure and envelope images are renders) | `data/simulation/vdw_structures/*.xyz` |
| fig. S3 | `figures/figS3.py` | as Fig. 3F |
| fig. S5 | `figures/figS5.py` | `data/experimental/Fig2B.csv`; table S2 values in `figures/figstyle.py` |
| fig. S6 | `figures/figS6.py` | `data/experimental/enmr/enmr_{1H,7Li,19F}.png` (instrument phase-shift plots) |
| fig. S8 | `figures/figS8.py` (with `rdf_cn.py`) | `data/simulation/conc_ladder/litpfa_fm_{n1,n2,n4,n8}/rdf.txt`, `litpfa_fm_n8/ions.lammpstrj` |
| fig. S9 B–D | `figures/figS9.py` | `data/simulation/opes2d_surfaces.npz` |
| fig. S11 | `figures/figS11.py` | `data/simulation/opes_pmf/*.npz` (14 salt–solvent pairs) |
| fig. S14 | `figures/figS14.py` | `data/experimental/figS14AB.csv`, `figS14C.csv`; fit uncertainties in `data/experimental/raw/figS14C_fit_uncertainty.csv` |
| fig. S15 | `figures/figS15.py` | `data/experimental/figS15A.csv`, `figS15B.csv`, `figS15C.csv` |
| fig. S17 | `figures/figS17.py` | `data/experimental/figS17/figS17_rate_cycling.csv` |
| fig. S19 | `figures/figS19.py` | `data/simulation/phonons/<salt>/mesh.yaml`, `thermal_properties.yaml` |

Figs. S7 and S16 are given as data only: `data/experimental/figS7.csv` and `figS16.csv` hold the plotted coordinates
read from the published vector graphics, and `data/experimental/raw/figS7_transference.csv` and
`figS16_cyclic_voltammetry.csv.gz` hold the measured values behind them. The fitted XPS components behind fig. S15 are in
`data/experimental/figS15/`.

Figs. S2 and S20 are also given as data only. `data/simulation/esp/{BF4,PF6,FSI,TFSI,TPFA}.fchk.gz` are the Gaussian 16
formatted checkpoint files of the five isolated anions (optimized geometry, wavefunction and density) behind the
electrostatic-potential maps of fig. S2: BF₄⁻, PF₆⁻, FSI⁻ and TFSI⁻ at B3LYP/6-31+G(d,p), TPFA⁻ at HF/6-31G(d).
`data/simulation/rdf/figS20_rdf.csv` holds the g(r) curves of fig. S20, from the dilute LiTPFA simulations with the GAFF
force field that supplied the Li⁺ shells (see Van der Waals volumes below) (columns `r_A`, `g_ACN`, `g_DME`, `g_Me2O`,
`g_DEE`, `g_FM`; r from 0.05 to 10.05 Å in 0.10 Å steps); the first local minimum after the main peak of each column is
the distance marked in the figure.

`data/experimental/raw/` holds the measured values as exported by the experimental groups, with neutral column names:

| file | figure | columns |
|---|---|---|
| `Fig4A_cell_voltage.csv.gz` | Fig. 4A | electrolyte, time (h), cell voltage (V), Li‖Li symmetric cells |
| `Fig4B_overpotential.csv` | Fig. 4B | electrolyte, current density (mA cm⁻²), average overpotential (mV) |
| `Fig4D_rate_cycling.csv` | Fig. 4D | electrolyte, cycle, specific discharge capacity (mAh g⁻¹), Coulombic efficiency (%), current density (mA cm⁻²) |
| `Fig4E_voltage_profiles.csv` | Fig. 4E | temperature (°C), charge and discharge rate, branch, point, specific capacity (mAh g⁻¹), voltage (V) |
| `Fig4F_voltage_profiles.csv` | Fig. 4F | as Fig. 4E, at −40 °C and four charge rates |
| `figS7_transference.csv` | fig. S7 | solvent, temperature (°C), Li⁺ transference number |
| `figS14C_fit_uncertainty.csv` | fig. S14C | electrolyte, temperature (°C), fitted element, value (Ω), fit uncertainty (Ω) |
| `figS16_cyclic_voltammetry.csv.gz` | fig. S16 | electrolyte, cycle, point, voltage (V vs Li/Li⁺), current (µA), current density (µA cm⁻²) |
| `figS18_extended_cycling.csv` | fig. S18 | electrolyte, temperature (°C), charge and discharge rate, cycle, specific discharge capacity (mAh g⁻¹), Coulombic efficiency (%) |

The plotted values in the Data S1 sheets agree with these files to the resolution of the published vector graphics.
Cycles whose Coulombic efficiency lies outside the plotted axis range in Fig. 4D are listed with their measured values.

Image-only panels, not included: Fig. 1A (schematic), Fig. 2A (apparatus schematic), Fig. 3B and 3C (molecular
renders; Fig. 3C shows configurations from the runs behind Fig. 3D), Fig. 4C (cryo-FIB/SEM micrographs), fig. S9A
(simulation-cell render), figs. S4, S10 and S12 (photographs), fig. S13 (equivalent circuit).

`data/experimental/figS17/figS17_rate_cycling.csv` holds the rate test of fig. S17: electrolyte, cycle, charge rate (C),
specific discharge capacity (mAh g⁻¹), Coulombic efficiency (%), discharge current density (mA cm⁻²) and `ce_shown`. The
areal loadings in the figure are 10 × the mean C/10 discharge current density of cycles 1–3, with 1C = 200 mA g⁻¹.
Coulombic efficiency is drawn only where the discharge capacity is at least 5 % of the cell's C/10 capacity
(`ce_shown` = 1); every value is in the file.

Not included and available from the authors on request: the eNMR phase-shift points behind fig. S6 and the PFG-NMR
spectra behind table S2.

## Processing scripts

The raw trajectories, COLVAR files and 2PT/3PT window output are not included (several gigabytes); the scripts take
their locations as arguments. The scripts marked * run from the deposit as it stands and rewrite files it ships.

| script | produces | from |
|---|---|---|
| `scripts/analyse_opes_pmf.py` | `data/simulation/opes_pmf/<salt>_<solvent>.npz` (Fig. 3E, fig. S11) | OPES COLVAR files of the bias groups |
| `scripts/opes_assoc.py` | `data/tables/opes_assoc.csv` (table S4) | OPES COLVAR files |
| `scripts/analyse_opes2d.py` | `opes2d_surfaces.npz` (fig. S9), ΔG_C | interface OPES COLVAR and KERNELS files |
| `scripts/species_count.py` | `data/simulation/species_count/species_count.json` (the species counts of `data/simulation/snapshot_stats.json`, Fig. 3D) | ion trajectories of the bulk runs |
| `scripts/desolv_pmf.py` | ΔG1 of Li⁺ desolvation | desolvation OPES COLVAR files |
| `scripts/xpt_solvation.py` * | `data/simulation/xpt_dgmix_293K.csv`, `xpt_dgmix_anion_clamp_293K.csv` (Fig. 3F, fig. S3) | `data/simulation/xpt/` (2PT/3PT window entropies and energies, group/group energies, neat and vacuum references) |
| `scripts/crystal_thermo.py` * | `data/tables/crystal_dissociation.csv`, `ideal_mixing.csv`, `crystal_thermo.csv`, `gas_ion_thermo.csv` (tables S8, S10, S12–S14) | `data/simulation/phonons/`, `data/simulation/qchem_freq/` (needs PyYAML) |
| `scripts/pairing_cycle.py` * | `data/simulation/pairing_cycle.json` (fig. S3) | the xPT tables, the thermochemistry tables and `data/simulation/opes_pmf/` |
| `scripts/smd_harvest.py` | `data/tables/smd_solvation.csv` (table S15) | Q-Chem SMD outputs |
| `scripts/cluster_be.py` | Li⁺–solvent binding energies (table S3, Fig. 3A) | Q-Chem cluster outputs |
| `scripts/table_s5.py` * | `data/tables/table_s5.csv` (table S5, shipped as `tableS5.csv`) | constants in the script |
| `scripts/li_phi_table.py` | `forcefield/interface/li_phi.table` | `Li_v2.eam.fs` |
| `scripts/vdw_volumes.py` * | `data/tables/vdw_volumes.csv`, `vdw_grid_convergence.csv`, `vdw_k.csv` (fig. S1, tables S6 and S7; the k of Figs. 1B and 2B) | `data/simulation/vdw_structures/` |


## Calculation notes

- Ion-pair potentials of mean force (Fig. 3E, fig. S11, table S4): frames of every walker are reweighted by
  exp[(V_bias − c(t))/k_BT] after discarding the first 30 %; W(r) = −k_BT ln P(r) + 2k_BT ln r on 0.1 Å bins, zeroed on
  10–12 Å; the band is the standard deviation over independent bias groups. K_a = 4π∫ r² exp(−W/k_BT) dr up to
  4.25 Å and ΔG_assoc = −k_BT ln(K_a c°). The table S4 uncertainties add in quadrature the spread over the bias groups,
  the spread over the four quarters of each trajectory and the spread over the choice of zeroing window; the LiBF₄
  values carry the spread over the bias groups alone.
- Solvation and mixing free energies (Fig. 3F, fig. S3): one Li⁺–anion pair per cell in 313 to 641 solvent molecules (641 in fluoromethane);
  2PT/3PT entropies of the ions, the ion solvation shells and the bulk solvent, referenced to a neat-solvent cell at
  the same volume and to the isolated anion in vacuum; ΔH from group/group ion–solvent energies plus the solvent
  reorganization energy. ΔG_mix = ΔG_solv + ΔG_diss − TΔS_ideal at 293.15 K. The second table removes the
  finite-subset offset of the shell entropies with equal-number random subsets of the neat solvent. Values plotted are
  medians over the cells; the Li⁺–anion interaction energy of each cell is added to ΔG_solv and ΔG_mix (Fig. 3F).
- Dissociation free energies (tables S8–S14): crystals from VASP relaxations and phonopy harmonic phonons at
  293.15 K, gas ions as rigid-rotor harmonic ideal gases at 1 bar.
- Interface free energies (fig. S9): two-dimensional OPES on the Li⁺ height above the Li(100) face and the Li⁺–F(FM)
  coordination number, reweighted after a 30 % discard per walker; region levels C = [1.4, 3], L1 = (3, 6] and
  R = [11, 14] Å.
- Van der Waals volumes (fig. S1, tables S6 and S7): the union of atom-centred spheres (radii in Å: H 1.20, Li 1.82,
  B 1.92, C 1.70, N 1.55, O 1.52, F 1.47, Al 1.84, P 1.80, S 1.80) on a cubic voxel grid, a voxel counting when its
  centre lies inside a sphere, with a 2h margin; V = N h³ at h = 0.035 Å, with h = 0.12, 0.08 and 0.05 Å and a
  100-point-per-axis grid (5 Å padding, `volume_grid100_A3`) for comparison. r = (3V/4π)^(1/3) and k = r₋/r₊ from the
  unrounded volumes. One structure per species; the Li⁺ shells come from dilute LiTPFA MD (one formula unit per 150
  solvent molecules) in each solvent, run with the GAFF force field.
  `python3 scripts/vdw_volumes.py --openbabel <path to obabel>` also prints the bonded fragments of each structure
  (Li⁺ bonds removed), a check that every solvent molecule is complete; the volumes do not use it.
- Ion aggregation (Fig. 3D): Li⁺ and an anion are in contact when any Li–anion atom pair is closer than 3.0 Å; the
  species count is the number of connected ionic clusters among the 16 ions, over the second half of each run.

## Inputs and force fields

`inputs/README.md` lists every example input with its launch line. `forcefield/merged/` holds the force fields
exactly as the simulations read them (comment text removed), one parameter file and one pair-coefficient include per
system. `forcefield/README.md` describes the file formats, states which file each figure uses, and maps every type
and pair row to the single-component parameter file it was taken from (`forcefield/component_map.csv.gz`; the
components are not distributed separately) or lists it as a term that exists only in the merged files. `forcefield/interface_controls/` holds, in the same form, the
two interface cells of section M.3 that no figure draws (the LiPF₆ contact-pair control and the acetonitrile set). The example decks have not been re-run with the released
samQEq 1.0.1 build.

The Li embedded-atom potential of the interface cell, `Li_v2.eam.fs` (A. Nichol and G. J. Ackland, Phys. Rev. B 93,
184101 (2016)), is not included. The NIST Interatomic Potentials Repository distributes it
(https://www.ctcms.nist.gov/potentials/entry/2016--Nichol-A-Ackland-G-J--Li/), and the file distributed there is byte-identical to the file the simulations read:

MD5 checksum of `Li_v2.eam.fs`: `fa7cc2b5c0bc41e6ee0b48b09e41464b`

Place it at `forcefield/interface/Li_v2.eam.fs` before running `inputs/lammps/interface_opes2d/`. The deck also reads the table `forcefield/interface/li_phi.table`, which is not included. It holds the pair term of the same potential, and `python3 scripts/li_phi_table.py forcefield/interface/Li_v2.eam.fs > forcefield/interface/li_phi.table` regenerates it from that file.

## Versions

- **1.0.2**: Fig. 1B plots no LiBF4/DEE conductivity. LiBF4 is insoluble in DEE (table S2), so `data/experimental/Fig1B.csv`
  and `data/data_S1.xlsx` carry no row for it, and Fig. 1 is re-rendered from them. No other value changes. Author order:
  A. Liu first, G. Raghavendran second.
- **1.0.1**: dielectric constants follow table S1 (Figs. 1B, 1C and 2C).
- **1.0.0**: first release.

## Licence

CC BY 4.0 (see `LICENSE` at the repository root).
