# 3PT — data and analysis scripts

Data, post-processing scripts, and run templates for:

**"An anharmonic liquid-entropy functional from the Mori–Zwanzig memory kernel"**
(Three-Phase Explicit Anharmonic Thermodynamics, 3PT).

**Archived on Zenodo:** https://doi.org/10.5281/zenodo.21447740

The 2PT/3PT entropy code itself (py-xPT) lives in a separate repository:
<https://github.com/atlas-nano/codes/tree/main/py-xPT>. The `.ini` control files in
`pyxpt_control/` are run with that code.

## Layout

```
3pt-construction/
  figures/          analysis + figure-generation scripts (self-contained: numpy/scipy/matplotlib)
  data/             post-processed inputs the scripts read
    lj_fullgrid_cage.csv      71-point LJ (ρ*,T*) grid, all entropy functionals (= 8 fs)
    lj_grid_8fs.csv,_64fs.csv timestep-resolved LJ grids (cage-efficiency figure)
    fig_d_data.npz            d=2 / metals p* calibration
    cage_sens/{8fs,64fs,metals,water}/cin_*.npz   per-state kernel+DoS arrays
    lj_4d/*_result.npz        4D-LJ isotherm/isobar analysis results
    pwr/*.pwr                 gas|cage|solid DoS decompositions (py-xPT output)
  pyxpt_control/    py-xPT .ini control files (LJ grid, metals, water; lj_4d/ engine .ini)
  lammps/           LAMMPS data files, EAM potentials, input templates; 4D MD engine source
  reproduce.py      regenerate every figure/table
```

## Quickstart

```bash
pip install numpy scipy matplotlib
python reproduce.py            # regenerate all figures + Table S1 into figures/
# or individually, e.g.:
python figures/make_fig_lj_cage.py
```

## Figure / Table → script + data

| Manuscript item | Script | Input data |
|---|---|---|
| Fig. 1 (signed bias), Fig. 2 (per-state ΔS*) | `figures/make_fig_lj_cage.py` | `data/lj_fullgrid_cage.csv` |
| Fig. 3 (liquid-metal accuracy) | `figures/make_fig_metals_accuracy.py` | in-script (Sun 2017 TI + this work) |
| Fig. 4 (gas\|cage\|solid DoS) | `figures/make_fig_cage_dos.py` | `data/pwr/*.pwr` |
| Fig. 5 (cage-efficiency universal law) | `figures/plot_universal_with_water.py` | `data/cage_sens/*`, `data/lj_grid_{8fs,64fs}.csv` |
| Fig. 6 (p=1/d cross-dimensional) | `figures/make_fig_d.py` | `data/fig_d_data.npz`, `data/lj_fullgrid_cage.csv` |
| Fig. 7 (state-function / path independence) | `figures/make_fig_pathindep.py` | `data/lj_fullgrid_cage.csv` (+ `mbwr_lj.py`) |
| Table I, Table S1 (full LJ grid) | `figures/gen_table_s1.py` | `data/lj_fullgrid_cage.csv` |
| SI Fig. (d=4 prefactor test) | `figures/fig_d4_prefactor.py` | `data/lj_4d/*_result.npz` |
| SI Fig. (cage f0/nu_c sensitivity) | `figures/sweep_gate_filter.py` | `data/cage_sens/8fs/cin_*.npz` |

Helper modules: `figures/kernel_tools.py` (standalone Volterra-kernel inversion +
cage-entropy routines, vendored from py-xPT so the figures need no external code),
`figures/cage_decompose_compare.py`, `figures/strengthen_nonmarkovian.py`,
`figures/mbwr_lj.py` (Johnson–Zollweg–Gubbins LJ EOS).

## Notes

- **Trajectories are not included.** The raw MD trajectories (~200 GB) are excluded
  for size; the shipped `data/` are the post-processed inputs (per-state kernel/DoS
  arrays, DoS decompositions, and the tabulated entropies) sufficient to regenerate
  every figure. Trajectories are available from the authors on request.
- The **d=4 decisive-point** marker (ρ*=1.10, T*=1.0) in the SI prefactor figure
  requires the 4D production trajectory and is omitted here; the figure regenerates
  from the isotherm/isobar result files.
- **LAMMPS files are templates** (one representative input per system class, plus the
  data files / EAM potentials and the custom 4D MD engine source under `lammps/lj_4d/`).
- Water standard molar entropies (manuscript Tables III–IV) are assembled from the
  FEP/Frenkel–Ladd routes described in the manuscript SI; the cage contributions come
  from the `data/pwr/` decompositions.
