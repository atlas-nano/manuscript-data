# Example inputs

There is one example input for each type of calculation in the paper. The LAMMPS decks use samQEq 1.0.1 syntax
(`fix qeq/sam`, `coul/shield/intra … all gaussian`, `fix_modify q shield gaussian`, `pppm/samqeq`) and read their
force field from `../forcefield/` through relative paths, so they run from their own directory.

Files ending in `.gz` are gzip-compressed. Run `gunzip` on them before you launch a deck. The decks name the
uncompressed files.

In the launch lines, `lmp` is a LAMMPS build that includes samQEq 1.0.1, PLUMED (with OPES) and `fix xpt`. The rank
counts are the ones used in production.

## LAMMPS

| directory | calculation | launch |
|---|---|---|
| `lammps/bulk_fqmd/` | Bulk fluctuating-charge MD of LiTPFA in fluoromethane, 8 ion pairs (Fig. 3C, D). Capped-displacement relaxation at 5 K, heating to 293 K, NPT at 90 bar, then NVT production with an ion-only trajectory. `in.litpfa.cont` continues from the restart file of the first segment, renamed `litpfa.restart`. | `srun -n 16 lmp -in in.litpfa` |
| `lammps/opes_ionpair/` | Ion-pair association free energy (Fig. 3E, fig. S11, table S4). The bias is OPES_METAD on the Li⁺–anion centre-of-mass distance (`plumed.dat`). Four walkers share one bias through `WALKERS_MPI`, the LAMMPS `world` variables give each walker its own seed, and production used 8 such independent groups. | `srun -n 64 lmp -partition 4x16 -in in.litpfa_fm` |
| `lammps/xpt_cell/solution/` | 2PT/3PT solvation cell (Fig. 3F, fig. S3): one LiTPFA ion pair in 641 fluoromethane molecules, with `fix xpt` on the Li⁺ and anion solvation shells, the bulk solvent and the ions. | `srun -n 32 lmp -in in.litpfa_fm` |
| `lammps/xpt_cell/neat/` | Neat-solvent 2PT/3PT reference at the volume of the solution cell. | `srun -n 32 lmp -in in.fm_neat` |
| `lammps/xpt_cell/vacuum/` | Gas-phase reference of the TPFA⁻ anion: 30 independently thermalized windows with explicit equipartition, each analysed by `fix xpt` (2PT). `rotaxes.inc` holds the normal deviates that set the initial angular momentum of each window. | `srun -n 24 lmp -in in.tpfa_vac` |
| `lammps/opes_desolv/` | Li⁺ desolvation free energy ΔG1: OPES_METAD on the Li⁺–fluorine coordination number `cn` (`plumed.dat`), 4 walkers per group. | `srun -n 128 lmp -partition 4x32 -in in.desolv_fm` |
| `lammps/conc_ladder/` | Concentration ladder (fig. S8), lowest rung: one ion pair in 641 fluoromethane molecules, NPT then NVT with Li⁺ radial distribution functions. | `srun -n 32 lmp -in in.litpfa_fm` |
| `lammps/interface_opes2d/` | Li(100)/electrolyte interface OPES (fig. S9, ΔG_C): 2D bias on the Li⁺ height above the electrode face and on `cn`. Sixteen walkers share one bias, and the `world` variables at the top of the deck give each walker its Li⁺ index, face and start frame. The start frame is read from `traj.lammpstrj` (the walkers marked `END` start from `equil_final.data`). `plumed.w1.s1.dat` is the input for walker 1, segment 1. The other walkers' files differ only in the biased Li⁺ atom index and in which face distance (`dt` or `db`) defines `h`. The deck reads the Li embedded-atom potential from `../../../forcefield/interface/Li_v2.eam.fs`, which is not included: obtain it as `forcefield/README.md` describes and place it there. The deck also reads `../../../forcefield/interface/li_phi.table`, which is not included either: regenerate it from the same file with `python3 scripts/li_phi_table.py forcefield/interface/Li_v2.eam.fs > forcefield/interface/li_phi.table` from the deposit root. | `srun -n 512 lmp -partition 16x32 -in in.opes2d -var CONT 0 -var SEG 1` |
| `lammps/neat_solvent/density/` | Neat-solvent liquid density (SI section K.2): 256 acetonitrile molecules, NPT at 293.15 K and 1 bar with the 12.0 Å shield cutoff. The production length is set on the command line. | `srun -n 16 lmp -in in.density -var NRUN 10000000` |

The 2PT/3PT decks use `fix_modify q pinfrag pin`. The charge-conservation fragments come from the integer per-atom
property `i_pin`, which holds the molecule ids of the data file, so the anion keeps its sub-fragments. The molecule-id
column is then reset by `set group anionall mol …`, so that `fix xpt` treats the anion as one whole molecule. The two
roles of the molecule id are thus kept apart, and the charge solve is the same as without the reset.

## Packmol

`packmol/pack.inp` packs 641 fluoromethane molecules, one TPFA⁻ and one Li⁺ into a 38 Å cubic cell (run as
`packmol < pack.inp`). The molecule templates are `_solv.xyz`, `_an.xyz` and `_li.xyz`. Packmol uses only their
coordinates, so the single element label carried by every atom of the solvent template has no effect. Atom types,
charges and topology are assigned afterwards when the LAMMPS data file is written. Every electrolyte box in the paper
was packed with Packmol.

## Q-Chem

| file | calculation |
|---|---|
| `qchem/cluster_be/li_fm_cluster.in` | Li⁺(CH₃F)ₙ cluster energies at ωB97X-D/def2-TZVPD (table S3): fluoromethane, Li⁺, and the n = 1–6 clusters at their optimized geometries |
| `qchem/cluster_be/li_fm_ligands_at_cluster.in` | the n ligands alone at their in-cluster coordinates (E(ligandsₙ@cluster) of Eq. S1) |
| `qchem/smd/tpfa_me2o_smd.in` | SMD solvation of TPFA⁻ in dimethyl ether at ωB97X-D/def2-TZVPD (table S15): gas-phase job, then the SMD job with the user-defined solvent descriptors in `$smx` (ε 6.2, n 1.30, α 0, β 0.45, γ 17.9915, φ 0, ψ 0) |
| `qchem/anion_freq/bf4_freq.in` | harmonic frequencies of BF₄⁻ at ωB97M-rV/def2-TZVP (tables S12–S14) |
| `qchem/chelpg_field/fm_chelpg.in` | B3LYP/def2-TZVP optimization with Becke–Johnson-damped Grimme dispersion (`DFT_D d3_bj`), and ChElPG charges of fluoromethane (charge targets of SI section K.2) |
| `qchem/chelpg_field/bf4_field_xm.in` | one finite-field job of the anion polarizability (SI section K.3): BF₄⁻ in a uniform field of −0.001 au along x at ωB97M-rV/def2-TZVPD. The other five jobs differ only in the field component and its sign. `SCF_GUESS READ` takes the guess from the zero-field job run in the same scratch directory. |

## VASP and phonopy

| directory | calculation |
|---|---|
| `vasp/relax/` | full cell and ion relaxation of LiTPFA (SI section A), PBE with zero-damped Grimme dispersion (`IVDW = 11`), 520 eV, Γ-centred 2×2×1 mesh |
| `vasp/phonon/` | single-point forces on the 348 displaced 232-atom cells (SI sections B and G), at Γ, `EDIFF = 1E-8`. `phonopy_command.txt` gives the phonopy commands that generate the displacements, collect the forces and evaluate the thermal properties on an 8×8×8 mesh. |

POTCAR files are not included. The PAW potentials are, from the VASP PBE set: `PAW_PBE Al 04Jan2001`, `PAW_PBE C 08Apr2002`,
`PAW_PBE F 08Apr2002`, `PAW_PBE Li_sv 10Sep2004` and `PAW_PBE O 08Apr2002`, in the order of the POSCAR species line
(Al C F Li O).
