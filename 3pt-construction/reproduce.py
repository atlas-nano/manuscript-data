#!/usr/bin/env python3
"""Regenerate every figure and table for the 3PT manuscript.

Runs each script in figures/ in turn; outputs are written into figures/.
Requires: numpy, scipy, matplotlib.
"""
import subprocess
import sys
from pathlib import Path

FIG = Path(__file__).resolve().parent / "figures"

SCRIPTS = [
    "make_fig_lj_cage.py",        # Fig. 1, Fig. 2
    "make_fig_d.py",              # Fig. 6
    "make_fig_pathindep.py",      # Fig. 7
    "make_fig_metals_accuracy.py",  # Fig. 3
    "make_fig_cage_dos.py",       # Fig. 4
    "fig_d4_prefactor.py",        # SI: d=4 prefactor test
    "sweep_gate_filter.py",       # SI: cage f0/nu_c sensitivity (fig_SI_cage_sensitivity)
    "plot_universal_with_water.py",  # Fig. 5 (universal cage-efficiency law)
    "gen_table_s1.py",            # Table I / Table S1
]

if __name__ == "__main__":
    fail = 0
    for s in SCRIPTS:
        print(f"==> {s}")
        r = subprocess.run([sys.executable, str(FIG / s)], cwd=str(FIG))
        if r.returncode != 0:
            print(f"    FAILED: {s}")
            fail += 1
    print(f"\nDone ({len(SCRIPTS) - fail}/{len(SCRIPTS)} succeeded).")
    sys.exit(1 if fail else 0)
