# Force field

Fluctuating-charge (samQEq) force field for the Li-salt electrolytes of this study, in LAMMPS `units metal`.
Every file in `merged/` is the file the production simulations read, with comment text removed; no number,
type index or line order differs from the run copy.

## Files

| path | content |
|---|---|
| `merged/<salt>_<solvent>.param` | `fix qeq/sam` parameter file for one electrolyte (salt = `litpfa`, `lipf6`, `libf4`; solvent = `fm`, `acn`, `me2o`, `dme`, `dee`) |
| `merged/paircoeff_<salt>_<solvent>.inc` | the matching `pair_coeff` lines (`include` file) |
| `merged/paircoeff_litpfa_fm_xpt.inc`, `merged/paircoeff_lipf6_fm_xpt.inc` | the pair set of the neat-fluoromethane reference cells of the 2PT/3PT analysis (see "Which file each figure uses") |
| `merged/interface_cell.param`, `merged/paircoeff_interface_cell_geometric.inc` | Li(100) slab with LiTPFA in fluoromethane |
| `interface_controls/interface_cell_lipf6.param`, `interface_controls/paircoeff_interface_cell_lipf6.inc` | Li(100) slab with LiPF₆ in fluoromethane (the LiPF₆ contact-pair control of section M.3) |
| `interface_controls/interface_cell_acn.param`, `interface_controls/paircoeff_interface_cell_acn.inc` | Li(100) slab with LiTPFA in acetonitrile (the acetonitrile interface set of section M.3) |
| `interface/Li_v2.eam.fs` | not included: the Li embedded-atom potential of A. Nichol and G. J. Ackland, Phys. Rev. B 93, 184101 (2016), used unmodified for the Li(100) electrode (see "Li embedded-atom potential" below) |
| `interface/li_phi.table` | not included: the pair term phi(r) of `Li_v2.eam.fs` alone, tabulated for the Li⁺–Li⁰ interaction (`pair_style table`, keyword `LI_PHI`); `scripts/li_phi_table.py` writes it from `Li_v2.eam.fs` (see "Li embedded-atom potential") |
| `component_map.csv` | one line per type row and per pair row of every merged file, naming the single-component parameter file and row it was taken from |

Files ending in `.gz` are gzip-compressed; decompress them before use (`gunzip *.gz`).

## Units and formats

- `.param`: line 1 is the global solver line `gamma_align kappa_bond r_ov r_loc lr_alpha lr_ewald`; each further line is
  `type chi eta Rc q0 eHOMO eLUMO`, with electronegativity chi and hardness eta in eV, Gaussian screening radius Rc in Å,
  reference charge q0 in e, and two unused placeholder columns (format of samQEq 1.0.1).
- `pair_coeff i j coul/shield/intra Rc_ij`: Gaussian-kernel screening radius in Å. The diagonal equals the Rc column of the
  `.param` file and the off-diagonal is the geometric mean sqrt(Rc_i Rc_j).
- `pair_coeff i j morse D alpha r_e`: well depth D in eV, alpha in Å⁻¹, r_e in Å, cutoff 12.0 Å from `pair_style`.
- `pair_coeff i j buck A rho C`: A in eV, rho in Å, C in eV Å⁶ (C = 0 everywhere; these rows are repulsive overlays).
- Atom types are numbered solvent first, then anion, then Li⁺ (and Li⁰ in the interface cell):

| system | solvent types | anion types | Li types (LiTPFA / LiPF₆ / LiBF₄) |
|---|---|---|---|
| `*_fm` | 1–5 | LiTPFA 6–10, LiPF₆ 6–7, LiBF₄ 6–7 | 11 / 8 / 8 |
| `*_acn` | 1–6 | LiTPFA 7–11, LiPF₆ 7–8, LiBF₄ 7–8 | 12 / 9 / 9 |
| `*_me2o` | 1–9 | LiTPFA 10–14, LiPF₆ 10–11, LiBF₄ 10–11 | 15 / 12 / 12 |
| `*_dme` | 1–16 | LiTPFA 17–21, LiPF₆ 17–18, LiBF₄ 17–18 | 22 / 19 / 19 |
| `*_dee` | 1–15 | LiTPFA 16–20, LiPF₆ 16–17, LiBF₄ 16–17 | 21 / 18 / 18 |
| `interface_cell` | 1–5 (fluoromethane) | 6–10 (TPFA⁻) | 11 Li⁺, 12 Li⁰ (electrode) |
| `interface_cell_lipf6` | 1–5 (fluoromethane) | 6–7 (PF₆⁻) | 8 Li⁺, 9 Li⁰ (electrode) |
| `interface_cell_acn` | 1–6 (acetonitrile) | 7–11 (TPFA⁻) | 12 Li⁺, 13 Li⁰ (electrode) |

## Which file each figure uses

| merged pair (`.param` + `.inc`) | used for |
|---|---|
| `litpfa_fm.param` + `paircoeff_litpfa_fm.inc` | aggregation runs (Fig. 3C, D), ion-pairing OPES (Fig. 3E, fig. S11, table S4), concentration ladder (fig. S8), Li⁺ desolvation free energy (ΔG1), 2PT/3PT solvation cells (Fig. 3F, fig. S3) |
| `litpfa_fm.param` + `paircoeff_litpfa_fm_xpt.inc` | neat-fluoromethane 2PT/3PT references of the LiTPFA cells (Fig. 3F, fig. S3), and the geometric-combining check of section F |
| `lipf6_fm.param` + `paircoeff_lipf6_fm.inc` | aggregation runs (Fig. 3C, D), ion-pairing OPES (Fig. 3E, fig. S11, table S4), 2PT/3PT solvation cells (Fig. 3F, fig. S3) |
| `lipf6_fm.param` + `paircoeff_lipf6_fm_xpt.inc` | neat-fluoromethane 2PT/3PT references of the LiPF₆ cells (Fig. 3F, fig. S3), and the geometric-combining check of section F |
| `libf4_fm` | ion-pairing OPES (fig. S11, table S4) and 2PT/3PT cells (Fig. 3F, fig. S3), one file for both |
| `litpfa_acn` | ion-pairing OPES, 2PT/3PT cells and Li⁺ desolvation in acetonitrile, one file for all three |
| `{litpfa,lipf6,libf4}_{acn,me2o,dme}`, `litpfa_dee`, `lipf6_dee` | ion-pairing OPES (fig. S11, table S4) and 2PT/3PT cells (Fig. 3F, fig. S3), one file for both |
| `libf4_dee` | 2PT/3PT cells (Fig. 3F, fig. S3) |
| `interface_cell` + `Li_v2.eam.fs` + `li_phi.table` | Li(100)/electrolyte OPES (fig. S9, ΔG_C, section M) |
| `interface_controls/*_lipf6` + `Li_v2.eam.fs` + `li_phi.table` | LiPF₆ contact-pair control (ΔG_L1, section M.3) |
| `interface_controls/*_acn` + `Li_v2.eam.fs` + `li_phi.table` | acetonitrile interface set (section M.3) |

In fluoromethane, the two LiTPFA pair files and the two LiPF₆ pair files differ only in their anion–solvent cross terms:

- `paircoeff_*_fm.inc` (Fig. 3C–F, figs. S3, S8 and S11, table S4, ΔG1 and the interface cell) uses the anchored
  PF₆⁻–fluoromethane cross terms (with the P–FM rows at D = 0). For
  TPFA⁻ it uses F–FM rows obtained by applying the same correction factors to TPFA⁻'s geometric baseline, which gives
  F(TPFA⁻)–F(FM) D = 0.00514 eV.
- `paircoeff_*_fm_xpt.inc` uses geometric means of the anion and fluoromethane diagonals for every
  anion–fluoromethane Morse row. For example, F(TPFA⁻)–F(FM) has D = 0.00339 eV and P–F(FM) has D = 0.00658 eV.
  It is read by the neat-fluoromethane reference cells of Fig. 3F and fig. S3, which hold no anion, so the
  anion–fluoromethane rows, the only rows in which the two files differ, do not act in them; and by the
  geometric-combining check of section F (the same sixteen LiTPFA and LiPF₆ solution cells run with geometric
  anion–fluoromethane depths).

The interface cell's pair file is named `_geometric` after its Li⁰ cross terms; its anion–fluoromethane rows are the
anchored rows of `paircoeff_litpfa_fm.inc`.

Both pair files in each pair use the same `.param` file. The solution-cell example in `inputs/lammps/xpt_cell/solution/`
reads `paircoeff_litpfa_fm.inc`, and the neat reference in `inputs/lammps/xpt_cell/neat/` reads `paircoeff_litpfa_fm_xpt.inc`.

## Interface control files

`interface_controls/` holds the files of the two interface datasets of section M.3 that are not drawn in a figure: the
LiPF₆ contact-pair control in fluoromethane (ΔG_L1) and the acetonitrile interface set. Each file is the file the
production simulations read, with comment text removed; from the LiPF₆ pair file two `variable` lines that hold no
parameter are also removed. No number, type index or line order differs from the run copy.

- `paircoeff_interface_cell_lipf6.inc`: its first 72 lines are the lines of `merged/paircoeff_lipf6_fm.inc`, in the
  same order; the remaining 20 lines are the Li⁰ (type 9) rows.
- `paircoeff_interface_cell_acn.inc`: its first 158 lines are the lines of `merged/paircoeff_litpfa_acn.inc`, in the
  same order; the remaining 30 lines are the Li⁰ (type 13) rows.

These four files are not covered by `component_map.csv`. Like `interface_cell`, both cells take their Li⁰–Li⁰
interaction from `Li_v2.eam.fs` and their Li⁺–Li⁰ short-range interaction from `li_phi.table` (see below).

## Provenance of the rows

The files in `merged/` are complete: the example inputs in `inputs/` read them, and no other parameter file is needed.
Each merged file was assembled from single-component parameter files (one per solvent, per anion and for Li⁺), which are
not distributed separately. `component_map.csv.gz` records, for every type row and pair row, the component file and row
it was taken from (column `component_file`).

Type rows: every type row is identical to a row of the component file of its species (solvent, anion or Li⁺).

The Li⁰ electrode row (type 12 of `interface_cell.param`) repeats the Li⁺ chi, eta and Rc with q0 = 0.

Pair rows: each pair row falls into one of three classes, and `component_map.csv` gives the class and the source file for every row.

1. **Identical to a component line.** The source is either the component's own pair file or a cross-term file, which
   holds the Li⁺–F(FM) and Li⁺–N(ACN) rows, the Li⁺–TPFA⁻ Morse and Born–Mayer rows, or the PF₆⁻–fluoromethane rows.
2. **Geometric mean of component diagonals.** Rc for `coul/shield/intra`, and D, alpha and r_e for `morse`, matching to the printed precision.
3. **Merged files only.** These rows exist only in the merged files:

| merged file | rows | term |
|---|---|---|
| `paircoeff_litpfa_fm.inc`, `paircoeff_interface_cell_geometric.inc` | `1 9`, `2 9`, `3 9`, `4 9`, `5 9` morse | F(TPFA⁻)–F, –C, –H(FM): D 0.005144 / 0.006019 / 0.000948 eV, alpha 1.84125 / 1.67602 / 1.98309 Å⁻¹, r_e 3.27486 / 3.59772 / 3.04063 Å |
| every `paircoeff_lipf6_*.inc` | Li⁺–P and Li⁺–F(PF₆⁻) morse | D 0.003888 / 0.00629 eV, alpha 2.00509 / 2.19596 Å⁻¹, r_e 3.2 / 2.734 Å |
| `paircoeff_interface_cell_geometric.inc` | `1 12` morse | Li⁰–F(FM) |
| `paircoeff_interface_cell_geometric.inc` | `3 12`, `4 12`, `5 12` morse and buck | Li⁰–H(FM) Morse plus a Born–Mayer wall A = 0.1 exp(1.65/0.2) eV, rho = 0.2 Å |

The interface cell takes its Li⁰–Li⁰ interaction and the Li⁰ embedding from `Li_v2.eam.fs` (`pair_style eam/fs`, mapped to
type 12 only), and its Li⁺–Li⁰ short-range interaction from `li_phi.table`.

## Li embedded-atom potential

`Li_v2.eam.fs` is not part of this deposit. It is the LAMMPS `eam/fs` file of the Li potential of A. Nichol and
G. J. Ackland, Phys. Rev. B 93, 184101 (2016), and the interface cell uses it unmodified. The NIST Interatomic
Potentials Repository distributes it (entry `2016--Nichol-A-Ackland-G-J--Li`, LAMMPS implementation
`2016--Nichol-A--Li--LAMMPS--ipr2`): https://www.ctcms.nist.gov/potentials/entry/2016--Nichol-A-Ackland-G-J--Li/
(file: https://www.ctcms.nist.gov/potentials/Download/2016--Nichol-A-Ackland-G-J--Li/3/Li_v2.eam.fs). The file distributed there is byte-identical to the file the simulations read
(786275 bytes):

MD5 checksum of `Li_v2.eam.fs`: `fa7cc2b5c0bc41e6ee0b48b09e41464b`

Place it at `forcefield/interface/Li_v2.eam.fs`, where `inputs/lammps/interface_opes2d/in.opes2d` reads it.
`li_phi.table` holds only the pair term of this potential, tabulated. It is not part of this deposit. With `Li_v2.eam.fs` in place, run
`python3 scripts/li_phi_table.py forcefield/interface/Li_v2.eam.fs > forcefield/interface/li_phi.table` from the deposit root. The script reads the
r*phi(r) column of the file, interpolates it with a cubic spline (not-a-knot ends) onto 8000 points from 0.5 to 4.563 Å, and writes phi(r) = (r*phi)/r and
the force -dphi/dr in LAMMPS metal units. It needs numpy and scipy.
