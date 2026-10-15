"""Figure 6: 3PT entropy differences against MBWR thermodynamic integration over all Lennard-Jones liquid pairs. Reads ../data/lj_grid.csv and uses mbwr_lj.py."""
import csv
import itertools
import math
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import mbwr_lj as mb

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, os.pardir, "data", "lj_grid.csv")
OUT = os.path.join(HERE, "out")
os.makedirs(OUT, exist_ok=True)

liq = []
for r in csv.DictReader(open(CSV)):
    rho, T = float(r["rho_star"]), float(r["T_star"])
    if 0.70 <= rho <= 1.02:
        liq.append((rho, T, float(r["3PT"])))

dS3, dSti = [], []
for (rA, TA, t3A), (rB, TB, t3B) in itertools.combinations(liq, 2):
    dS3.append(t3B - t3A)
    dSti.append(mb.s_total_reduced(rB, TB) - mb.s_total_reduced(rA, TA))
dS3, dSti = np.array(dS3), np.array(dSti)
resid = dS3 - dSti
rms = math.sqrt((resid ** 2).mean())

plt.rcParams.update({"font.size": 10, "axes.linewidth": 0.8, "figure.dpi": 150, "savefig.bbox": "tight"})
fig, (axA, axB) = plt.subplots(1, 2, figsize=(8.0, 3.6))
lim = [min(dSti.min(), dS3.min()) - 0.1, max(dSti.max(), dS3.max()) + 0.1]
axA.plot(lim, lim, "-", color="0.5", lw=1.0, zorder=1)
axA.plot(dSti, dS3, "o", ms=3, color="#55a868", alpha=0.5, zorder=2)
axA.set_xlim(lim)
axA.set_ylim(lim)
axA.set_aspect("equal")
axA.set_xlabel(r"$\Delta S_{\rm TI}$ (MBWR path integral, $k_B$/atom)")
axA.set_ylabel(r"$\Delta S_{\rm 3PT}=S^*_{\rm 3PT}(B)-S^*_{\rm 3PT}(A)$")
axA.set_title(f"(a) all {len(dS3)} liquid pairs", fontsize=9.5)
axA.grid(ls=":", lw=0.5, alpha=0.6)
axB.hist(resid, bins=30, color="#55a868", alpha=0.85)
axB.axvline(0, color="k", lw=0.9)
axB.set_xlabel(r"$\Delta S_{\rm 3PT}-\Delta S_{\rm TI}$  ($k_B$/atom)")
axB.set_ylabel("pairs")
axB.set_title(f"(b) residual: mean {resid.mean():+.3f}, RMS {rms:.3f}", fontsize=9.5)
axB.grid(axis="y", ls=":", lw=0.5, alpha=0.6)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig_path_independence.pdf"))
