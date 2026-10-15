# An anharmonic liquid-entropy functional from the Mori–Zwanzig memory kernel — figure data

Data and scripts that reproduce the figures of R. Buarque, M. Gascon and T. A. Pascal, *An anharmonic liquid-entropy
functional from the Mori–Zwanzig memory kernel*, The Journal of Chemical Physics (submitted). The entropies were
computed with py-xPT, https://github.com/atlas-nano/codes/tree/main/py-xPT.

```bash
python3 reproduce.py    # writes every figure to figures/out/
```

or one figure at a time from `figures/`:

```bash
cd figures
python3 fig_lj_bias.py              # Figures 1 and 2, from data/lj_grid.csv
python3 fig_metals_accuracy.py      # Figure 3, from data/metals.csv
python3 fig_cage_dos_partition.py   # Figure 4, from data/cage_dos_partition.npz
python3 fig_cage_efficiency.py      # Figure 5, from data/cage_efficiency.json
python3 fig_path_independence.py    # Figure 6, from data/lj_grid.csv (with mbwr_lj.py)
python3 fig_SI_cage_sensitivity.py  # Figure S1, from data/cage_sensitivity.csv
```

Requires Python 3 with numpy (2.0 or later) and matplotlib.

| file | contents |
|---|---|
| `data/lj_grid.csv` | 43 Lennard-Jones liquid states (0.70 ≤ ρ* ≤ 1.02): reduced density and temperature, MBWR entropy, 3PT entropy, and the signed deviation of the rigorous-HS, Lin-2003, Desjarlais, R2PT and 3PT entropies from MBWR, in k_B per atom |
| `data/metals.csv` | seven Sutton–Chen liquid metals: temperature, thermodynamic-integration entropy (Sun *et al.*, 2017), and the rigorous-HS, Lin-2003, Desjarlais, R2PT (Sun *et al.*, 2017) and 3PT entropies, in k_B per atom |
| `data/cage_dos_partition.npz` | for each panel of Figure 4: frequency (cm⁻¹), total, gas, cage and solid densities of states (DoF per particle per cm⁻¹), and the gas, cage and solid integrals over the full spectrum |
| `data/cage_efficiency.json` | cage efficiency and γ/Ω₀ for the Lennard-Jones liquids, the metals, the time-step ladders and water, with the fitted rate constant, its R², the mean entropy deficit and the resolution threshold Δt·Ω₀ |
| `data/cage_sensitivity.csv` | for 40 Lennard-Jones states: density, temperature, fluidicity, the cage entropy at the default settings, and the cage entropy (k_B per atom) at each fluidicity-gate value f₀ and each filter-cutoff scale |

License: CC BY 4.0.
