#!/usr/bin/env python
"""Metals accuracy figure for the 3PT/cage paper (the Table 2 visual).

Seven Sutton-Chen liquid metals vs the Sun et al. (2017) thermodynamic-integration
(TI) reference. All five entropy functionals are shown:
  rig-HS, Lin-2003 (+lnZ), Desjarlais, R2PT, 3PT.
rig-HS / Lin-2003 / Desjarlais / 3PT are computed by py-xPT from the same
metals_sc trajectories (this work); R2PT(delta=1.5) is Sun's tuned literature
benchmark (Sun 2017, Table II). The reference is Sun's TI.

Left  : per-metal signed dS = S_method - S_TI (k_B/atom), on a BROKEN y-axis so the
        large rig-HS deficit (~-0.3 k_B) does not squash the +-0.15 k_B spread of the
        corrected methods.
Right : RMS deviation vs TI by method (all five on one axis).
"""
import math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent

# k_B/atom.  TI and R2PT(1.5) from Sun et al. 2017, Table II.
# rig-HS, Lin-2003 (+lnZ), Desjarlais, 3PT: py-xPT, this work
# (metals_sc/*_{hs_rig,lin2003,desjarlais,3pt}.thermo; S_q[J/mol/K]/R).
METALS = ["Al", "Ag", "Ni", "Pd", "Rh", "Au", "Ir"]
TI    = np.array([ 9.25 , 10.24 , 10.31 , 11.08 , 11.10 , 11.71 , 12.25 ])
RIG   = np.array([ 9.089, 9.897 , 10.023, 10.725, 10.757, 11.365, 11.848])  # py-xPT rigorous-HS
LINZ  = np.array([ 9.547, 10.261, 10.478, 11.164, 11.169, 11.804, 12.250])  # py-xPT Lin-2003 (+lnZ)
DESJ  = np.array([ 9.576, 10.254, 10.489, 11.170, 11.163, 11.818, 12.242])  # py-xPT Desjarlais
R2PT  = np.array([ 9.40 , 10.25 , 10.38 , 11.05 , 11.09 , 11.73 , 12.19 ])  # Sun 2017 R2PT(1.5), tuned
S3PT  = np.array([ 9.314, 10.148, 10.294, 11.051, 11.078, 11.674, 12.160])  # py-xPT 3PT

# (label, data, color, marker) — colors match Figs 1/2
SERIES = [("rig-HS (uncorrected)",       RIG,  "#8172b3", "v"),
          ("Lin-2003 ($+\\ln Z$)",       LINZ, "#c44e52", "o"),
          ("Desjarlais (MF gas)",        DESJ, "#937860", "^"),
          ("R2PT ($\\delta{=}1.5$, Sun)", R2PT, "#dd8452", "s"),
          ("3PT (parameter-free)",       S3PT, "#55a868", "D")]

plt.rcParams.update({"font.size": 10, "axes.linewidth": 0.8,
                     "figure.dpi": 150, "savefig.bbox": "tight"})


def rms(x): return math.sqrt(np.mean(x**2))


def main():
    fig = plt.figure(figsize=(8.6, 4.1))
    gs = fig.add_gridspec(2, 2, width_ratios=[2.3, 1.0], height_ratios=[2.6, 1.0],
                          hspace=0.07, wspace=0.30)
    axT = fig.add_subplot(gs[0, 0])   # upper part of broken left panel
    axB = fig.add_subplot(gs[1, 0], sharex=axT)   # lower part (rig-HS deficit)
    axR = fig.add_subplot(gs[:, 1])   # RMS panel spans both rows
    x = np.arange(len(METALS))

    # --- left (broken y): per-metal signed deviation vs TI ---
    YT = (-0.20, 0.37)    # upper segment: corrected-method spread + band + rig-HS(Al) outlier
    YB = (-0.45, -0.25)   # lower segment: the rig-HS deficit cluster (6 metals)
    for ax in (axT, axB):
        ax.axhline(0, color="k", lw=0.9, zorder=2)
        for lab, S, col, mk in SERIES:
            ax.plot(x, S - TI, mk, color=col, ms=6, label=lab, zorder=3)
        ax.grid(axis="y", ls=":", lw=0.5, alpha=0.6)
    axT.axhspan(-0.05, 0.05, color="0.85", alpha=0.6, zorder=1, label=r"$\pm0.05\,k_B$")
    axT.set_ylim(*YT); axB.set_ylim(*YB)

    # hide the shared spine + add diagonal break marks
    axT.spines["bottom"].set_visible(False)
    axB.spines["top"].set_visible(False)
    axT.tick_params(labelbottom=False, bottom=False)
    d = 0.012
    kw = dict(transform=axT.transAxes, color="k", clip_on=False, lw=0.9)
    axT.plot((-d, +d), (-d, +d), **kw); axT.plot((1 - d, 1 + d), (-d, +d), **kw)
    kw.update(transform=axB.transAxes)
    axB.plot((-d, +d), (1 - 3*d, 1 + 3*d), **kw); axB.plot((1 - d, 1 + d), (1 - 3*d, 1 + 3*d), **kw)

    axB.set_xticks(x); axB.set_xticklabels(METALS)
    axB.set_xlabel("Sutton–Chen metal (increasing $T$, $S_{\\rm TI}$)")
    axT.set_title("Deviation from thermodynamic integration (Sun 2017)", fontsize=9.5)
    fig.supylabel("")  # avoid duplicate; set a single y-label spanning both
    axT.set_ylabel(r"signed $\Delta S = S_{\rm method}-S_{\rm TI}$  ($k_B$/atom)")
    axT.yaxis.set_label_coords(-0.10, 0.30)
    axT.legend(fontsize=7.5, loc="upper right", framealpha=0.9, ncol=2,
               columnspacing=1.0, handletextpad=0.3)

    # --- right: RMS by method (all five) ---
    for i, (lab, S, col, mk) in enumerate(SERIES):
        r = rms(S - TI)
        axR.bar(i, r, color=col, width=0.7, zorder=3)
        axR.text(i, r + 0.006, f"{r:.2f}", ha="center", va="bottom",
                 fontsize=8.5, fontweight="bold")
    axR.set_xticks(range(len(SERIES)))
    axR.set_xticklabels(["rig", "+lnZ", "Desj", "R2PT", "3PT"], fontsize=8, rotation=0)
    axR.set_ylabel(r"RMS deviation vs TI  ($k_B$/atom)")
    axR.set_title("RMS", fontsize=9.5)
    axR.set_ylim(0, 0.40)
    axR.grid(axis="y", ls=":", lw=0.5, alpha=0.6)

    for ext in ("png", "pdf"):
        fig.savefig(HERE / f"fig_metals_accuracy.{ext}")
    plt.close(fig)
    print("wrote fig_metals_accuracy.png/.pdf")
    for lab, S, col, mk in SERIES:
        print(f"  {lab.split(' (')[0]:22s} RMS={rms(S-TI):.3f}  mean={np.mean(S-TI):+.3f}")


if __name__ == "__main__":
    main()
