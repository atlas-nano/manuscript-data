#!/usr/bin/env python3
"""SI cage-sensitivity sweep (fig_SI_cage_sensitivity): cage entropy delta-S_cage
versus the fluidicity gate f0 and versus the spectral-filter cutoff scale nuc_scale,
across density, using the py-xPT cage_memory_entropy on the saved per-state inputs.
Self-contained: the cage function is vendored in kernel_tools.py (numpy/scipy only)."""
from pathlib import Path
import re, glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from kernel_tools import cage_memory_entropy

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"

files = sorted(glob.glob(str(DATA / "cage_sens" / "8fs" / "cin_*_grp1.npz")))

def rhoT(name):
    m = re.search(r"cin_r(\d+)_T(\d+)", name)
    r, T = m.group(1), m.group(2)
    return float(r[0] + "." + r[1:]), float(T[0] + "." + T[1:])

f0_grid  = np.array([1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1])
nuc_grid = np.array([0.25, 0.354, 0.5, 0.71, 1.0, 1.41, 2.0, 2.83, 4.0])

states = []
for fpath in files:
    d = np.load(fpath); rho, T = rhoT(fpath)
    args = dict(dt=float(d['dt']), C_scalar=d['C_scalar'], nu_cm=d['nu_cm'],
                dos_total=d['dos_total'], dos_gas=d['dos_gas'], T_K=float(d['T_K']),
                mass_amu=float(d['mass_amu']), vol_per_atom_A3=float(d['vol_per_atom']),
                dimension=int(d['dimension']), prefactor=1.0 / int(d['dimension']))
    base = cage_memory_entropy(**args, gate_f0=0.01, nuc_scale=1.0)
    if base is None:
        print(fpath, "-> None"); continue
    f0row = [cage_memory_entropy(**args, gate_f0=f0, nuc_scale=1.0) for f0 in f0_grid]
    ncrow = [cage_memory_entropy(**args, gate_f0=0.01, nuc_scale=nc) for nc in nuc_grid]
    f = float(np.trapezoid(d['dos_gas'], dx=(d['nu_cm'][1] - d['nu_cm'][0])) / int(d['dimension']))
    states.append((rho, T, base, np.array(f0row, float), np.array(ncrow, float), f))
    print(f"rho*={rho:.3f} T*={T:.2f}  dS_base={base:+.4f}  f={f:.3f}")

states.sort()
cols = cm.viridis(np.linspace(0, 1, len(states)))
fig, ax = plt.subplots(1, 2, figsize=(12, 4.8))
a = ax[0]
for (rho, T, base, f0r, ncr, f), c in zip(states, cols):
    a.plot(f0_grid, f0r / base, '-o', color=c, ms=3, label=f"rho*={rho:.2f} (f={f:.2f})")
    a.axvline(f, color=c, ls=':', lw=0.6, alpha=0.5)
a.axvline(0.01, color='k', ls='--', lw=1, label='f0=0.01 (default)')
a.set_xscale('log'); a.set_xlabel("gate parameter $f_0$")
a.set_ylabel(r"$\Delta S_{cage}(f_0)/\Delta S_{cage}(0.01)$")
a.set_title("(a) gate sensitivity across density\n(dotted = each state's fluidicity f; plateau for $f_0\\ll f$)")
a.legend(fontsize=6, ncol=2); a.grid(alpha=.3); a.set_ylim(0, 1.15)
b = ax[1]
for (rho, T, base, f0r, ncr, f), c in zip(states, cols):
    b.plot(nuc_grid, ncr / base, '-o', color=c, ms=3, label=f"rho*={rho:.2f}")
b.axvline(1.0, color='k', ls='--', lw=1, label='$\\nu_c=\\Omega_0/2\\pi c$ (parameter-free)')
b.set_xscale('log'); b.set_xlabel(r"filter cutoff scale  $\nu_c/(\Omega_0/2\pi c)$")
b.set_ylabel(r"$\Delta S_{cage}(\nu_c)/\Delta S_{cage}(1)$")
b.set_title("(b) filter-cutoff sensitivity across density\n(flat near 1 => entropy insensitive to the cutoff)")
b.legend(fontsize=6, ncol=2); b.grid(alpha=.3)
fig.suptitle("Cage $\\Delta S$ sensitivity to the fluidicity gate ($f_0$) and the spectral-filter "
             "cutoff ($\\nu_c$) across density (LJ, T*$\\approx$1.4)", fontsize=10)
fig.tight_layout(rect=[0, 0, 1, 0.94])
for ext in ("png", "pdf"):
    fig.savefig(HERE / f"fig_SI_cage_sensitivity.{ext}", dpi=140)

f0m = (f0_grid <= 0.03); ncm = (nuc_grid >= 0.5) & (nuc_grid <= 2.0)
print("\n=== insensitivity (max % dev of dS over the relevant range) ===")
for rho, T, base, f0r, ncr, f in states:
    df0 = 100 * np.nanmax(np.abs(f0r[f0m] / base - 1))
    dnc = 100 * np.nanmax(np.abs(ncr[ncm] / base - 1))
    print(f" rho*={rho:.3f}: f0<=0.03 dev {df0:5.1f}%   nuc 0.5-2x dev {dnc:5.1f}%")
print("saved fig_SI_cage_sensitivity.png/.pdf")
