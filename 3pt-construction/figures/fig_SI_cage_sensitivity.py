"""Figure S1: sensitivity of the Lennard-Jones cage entropy to the fluidicity gate and the spectral-filter cutoff. Reads ../data/cage_sensitivity.csv."""
import csv
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, os.pardir, "data", "cage_sensitivity.csv")
OUT = os.path.join(HERE, "out")
os.makedirs(OUT, exist_ok=True)

rd = csv.reader(open(CSV))
hdr = next(rd)
rows = np.array([[float(v) for v in r] for r in rd])
f0c = [i for i, h in enumerate(hdr) if h.startswith("cage_f0_")]
ncc = [i for i, h in enumerate(hdr) if h.startswith("cage_nuc_")]
f0_grid = np.array([float(hdr[i][len("cage_f0_"):]) for i in f0c])
nuc_grid = np.array([float(hdr[i][len("cage_nuc_"):]) for i in ncc])
states = [(r[0], r[1], r[3], r[f0c], r[ncc], r[2]) for r in rows]

cols = cm.viridis(np.linspace(0, 1, len(states)))
fig, ax = plt.subplots(1, 2, figsize=(12, 4.8))
a = ax[0]
for (rho, T, base, f0r, ncr, f), c in zip(states, cols):
    a.plot(f0_grid, f0r / base, '-o', color=c, ms=3, label=f"rho*={rho:.2f} (f={f:.2f})")
    a.axvline(f, color=c, ls=':', lw=0.6, alpha=0.5)
a.axvline(0.01, color='k', ls='--', lw=1, label='f0=0.01 (default)')
a.set_xscale('log')
a.set_xlabel("gate parameter $f_0$")
a.set_ylabel(r"$\Delta S_{cage}(f_0)/\Delta S_{cage}(0.01)$")
a.set_title("(a) gate sensitivity across density\n(dotted = each state's fluidicity f; plateau for $f_0\\ll f$)")
a.legend(fontsize=6, ncol=2)
a.grid(alpha=.3)
a.set_ylim(0, 1.15)
b = ax[1]
for (rho, T, base, f0r, ncr, f), c in zip(states, cols):
    b.plot(nuc_grid, ncr / base, '-o', color=c, ms=3, label=f"rho*={rho:.2f}")
b.axvline(1.0, color='k', ls='--', lw=1, label='$\\nu_c=\\Omega_0/2\\pi c$ (parameter-free)')
b.set_xscale('log')
b.set_xlabel(r"filter cutoff scale  $\nu_c/(\Omega_0/2\pi c)$")
b.set_ylabel(r"$\Delta S_{cage}(\nu_c)/\Delta S_{cage}(1)$")
b.set_title("(b) filter-cutoff sensitivity across density\n(ratio to the nominal cutoff; not flat)")
b.legend(fontsize=6, ncol=2)
b.grid(alpha=.3)
fig.suptitle("Cage $\\Delta S$ sensitivity to the fluidicity gate ($f_0$) and the spectral-filter "
             "cutoff ($\\nu_c$) across density (LJ, T*$\\approx$1.4)", fontsize=10)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(os.path.join(OUT, "fig_SI_cage_sensitivity.pdf"), dpi=140)
