# Tailoring Energy Landscapes for Vehicular Transport in Single-Ion Conducting Organo Ionic Solids
## Molecular Dynamics Simulations

**Archived on Zenodo:** https://doi.org/10.5281/zenodo.21448069

This repository contains classical molecular dynamics and metadynamics simulations for studying sodium ion (Na+) conductivity in glyme-based electrolytes, focusing on the NA-Orion conductor system.

## Repository Structure
## Simulation Systems

### Classical MD
- **glyme_tetramer_8_1**, **glyme_tetramer_6_1**, **glyme_tetramer_4_1**, **glyme_tetramer_2_1**: glyme:tetramer=8:1, 6:1, 4:1, 2:1, Na+ is the counterion
- **Temperatures**: 25°C (room temperature), 60°C, 100°C
- **Run types**: Short equilibration runs + extended production runs (long_run)


### Metadynamics
- **Systems**: 
  - glyme_tetramer_2_1_fix_momentum (glyme:tetramer= 2:1, represented as [1O:Na]G1 in the paper)
  - glyme_tetramer_8_1_fix_momentum (glyme:tetramer= 8:1, represented as [4O:Na]G1 in the paper)
- **Collective Variables**: Biased free energy calculations via COLVARS


## Input Files

Each simulation directory contains:
- `in.lammps` — LAMMPS input script for simulation
- `data.lammps` — Initial atomic coordinates and connectivity
- `*.slurm` — SLURM job submission scripts
- `in_restart.lammps` — Restart input for continuation runs (long_run)

### Force Field Files
- `na_conductor.frcmod` — AMBER force field bonded parameters

- `resp_full.mol2` - na_conductor structure file updated with the DFT partial charges
- `submit.sh` - The command to generate na_conductor.frcmod force field parameters
- `na_conductor.mol2` — NA-Orion conductor structure mol2 file(charge not updated)
- `na_conductor.bgf` — NA-Orion conductor structure bgf file (charge not updated)
- `gaff2.dat` - GAFF froce field parameters to use with na_conductor.frcmod, use both gaff2.dat and na_conductor.frcmod for the simulations

## Running Simulations

### Classical MD
```bash
# Navigate to simulation directory
cd classical/298K/glyme_tetramer_8_1/

# Run LAMMPS
sbatch lammps.lammps.slurm

# For long runs, use restart input
cd long_run/
sbatch lammps_restart.lammps.slurm


### MetaDynamics MD
# Run classical MD for smaller systems
cd metadynamics/glyme_tetramer_8_1_fix_momentum/classical/
sbatch lammps.lammps.slurm
# Navigate to metadynamics directory
cd metadynamics/glyme_tetramer_8_1_fix_momentum/meta/
sh mk-folders-cp.sh
sh submit.sh

## Analysis Scripts
All analysis scripts are in classical/analysis and metaD/analysis.
1. classical/analysis
- `clusters_classicfication.ipynb` - classify clusters into Structural/CIP/SSIP
- `ele_machine.ipynb` - extract clusters from trajectories, run before classifying clusters 
- `solvation_cli_na_only.py` - run solvation structure analysis and residence time analysis
- `transport_unfix_mom.ipynb` - run transport analysis including self-diffusivities and conductivites.