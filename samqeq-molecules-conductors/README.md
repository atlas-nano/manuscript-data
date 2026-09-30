# samQEq: unified charge equilibration from molecules to conductors — figure data

Data and scripts that reproduce the figures of N. Solan, D. Sun and T. A. Pascal, *samQEq: unified charge
equilibration from molecules to conductors*, Computer Physics Communications (submitted). The code is samQEq 1.0.0,
https://github.com/atlas-nano/codes/tree/main/samqeq, archived at https://doi.org/10.5281/zenodo.23043176.

```bash
cd figures
python3 fig_architecture.py    # Figure 1 (schematic)
python3 fig_callstructure.py   # Figure 2 (schematic)
python3 fig_conductor.py       # Figure 3, from data/dense_exp1.csv, offset.csv, exp2b.csv, wscan.csv
python3 fig_gpuscale.py        # Figure 4, from data/scaling.csv
python3 fig_tolscan.py         # Figure S1 (data listed in the script)
```

Requires Python 3 with numpy and matplotlib.

| file | contents |
|---|---|
| `data/dense_exp1.csv` | smallest eigenvalue of the constrained operator, Au(111) slab with a water film, by constraint, self-energy and gold hardness |
| `data/offset.csv` | water-to-metal charge and the largest single-molecule charge against an offset of the water electronegativities |
| `data/exp2b.csv` | charge above the slab mid-plane opposite a unit point charge at height z |
| `data/wscan.csv` | face and sub-surface charges of the slab in a uniform field against the self-energy width |
| `data/scaling.csv` | GPU throughput and time per task against node count |

License: CC BY 4.0.
