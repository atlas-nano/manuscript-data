# samQEq: charge equilibration for molecules and conductors — data

Input decks, analysis scripts and outputs behind N. Solan, D. Sun and T. A. Pascal, *samQEq: charge equilibration
for molecules and conductors*, Computer Physics Communications (submitted). The code is samQEq 1.0.0,
https://github.com/atlas-nano/codes/tree/main/samqeq, archived at https://doi.org/10.5281/zenodo.23043176.

## Where each result lives

| paper item | directory | harvest |
|---|---|---|
| §2.4 bare Au(111) slab in 0.5 V/Å; Fig. 3d (self-energy width) | `conductor_slab_field/` (base deck), `conductor_demo/runs/wscan/` | `conductor_demo/analyze_slabfield.py`, `analyze_wscan.py` |
| §7.1 molecules and a conductor in one solve; Fig. 3a–c; dense λmin | `conductor_demo/` (`runs/exp1`, `exp2`, `offset`) | `analyze_exp1.py`, `dense_exp1.py`, `analyze_exp2b.py`; `RESULT_conductor_demo_20260927.md` |
| §7.1 dynamics of the film (5 ps, six arms) | `conductor_demo/runs/dyn/` | `conductor_demo/analyze_dyn.py` → `dyn.csv` |
| §7.2 ported TIP4P-FQ liquid, 2.635 ± 0.002 D | `port_dipole/` | `python3 run_port.py --analyze` |
| §7.2 all-pairs-shielded liquid, 2.340 ± 0.001 D | `cluster_runs/allpairs_liquid_dipole_perlmutter/` | `anA.py` |
| Table 3 (kernel comparison, gas monomer) | `kernel_comparison_gas/` | `run_r2.py water`; `water_route2.json` |
| §7.3 ReaxFF charge step | `reaxff_port/` | `compare.py`; `RESULT_reactive_port_20260916.md` |
| §8 grid self-term sensitivity (51 frames) | `nve_clean/runs_selfterm/` | `nve_clean/run_selfterm.py --analyze` |
| Table 5 (charge dynamics, NVE) | `nve_clean/runs/`, `nve_clean/timing/` | `nve_clean/analyze_nve.py` → `nve.csv`; `RESULT_nve_clean_20260928.md` |
| §9 13,216-atom Cu/electrolyte cell scaling (Perlmutter CPU) | `cluster_runs/cu_cell_scaling_perlmutter/` | logs per rank count |
| Fig. 4 (GPU scaling, Perlmutter) | `gpu_scaling_perlmutter/` | `pm_harvest.py` → `scaling.csv` |
| Table 4 (device benchmark, Expanse V100 vs CPU) | `cluster_runs/gpu_benchmark_expanse/` | job logs |
| SI §S3 dense spectrum of the native water operator | `native_dense_diag/` | `hd_reconstruct.py` |
| SI Fig. S1 (charge-solve tolerance scan) | `tolerance_scan/` | `scan_npt.txt` (values inlined in `figures/fig_tolscan.py`) |
| Figures 1–4, S1 | `figures/` | `python3 fig_*.py` (run from `figures/`) |

The regression cases, the ReaxFF example and the TIP4P-FQ decks that ship with the code are in the code
repository and are not repeated here.

## Notes

- The `RESULT_*.md` files are the working records of each campaign, kept as written, including corrections of
  earlier measurements. Where a record and the paper differ, the paper is final.
- Decks and logs carry the absolute paths of the machines they ran on (a workstation, NERSC Perlmutter and
  SDSC Expanse); adjust them to rerun. LAMMPS restart files were removed (binary and machine-specific).
- The operator dumps of `conductor_demo/runs/exp1/*/opdump.txt.gz` are gzipped; `dense_exp1.py` reads them
  directly and reproduces `dense_exp1.csv` exactly.
- `conductor_demo/runs/exp2_nocross_20260927/` holds the image-energy runs made before the Au–ion shielding pair
  was added; the charges are unchanged, the energies are superseded by `runs/exp2/`.

## License

CC BY 4.0.
