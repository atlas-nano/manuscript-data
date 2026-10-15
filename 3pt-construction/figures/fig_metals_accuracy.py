"""Figure 3: entropy of seven Sutton-Chen liquid metals by method against thermodynamic integration. Reads ../data/metals.csv."""
import csv
import math
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, os.pardir, "data", "metals.csv")
OUT = os.path.join(HERE, "out")
os.makedirs(OUT, exist_ok=True)

rows = list(csv.DictReader(open(CSV)))
METALS = [r["metal"] for r in rows]


def col(name):
    return np.array([float(r[name]) for r in rows])


TI = col("S_TI")
SERIES = [("rig-HS (uncorrected)", col("rigHS"), "#8172b3", "v"),
          ("Lin-2003 ($+\\ln Z$)", col("Lin2003"), "#c44e52", "o"),
          ("Desjarlais (MF gas)", col("Desjarlais"), "#937860", "^"),
          ("R2PT ($\\delta{=}1.5$, Sun)", col("R2PT"), "#dd8452", "s"),
          ("3PT (parameter-free)", col("3PT"), "#55a868", "D")]

plt.rcParams.update({"font.size": 10, "axes.linewidth": 0.8,
                     "figure.dpi": 150, "savefig.bbox": "tight"})


def rms(x):
    return math.sqrt(np.mean(x ** 2))


def main():
    fig = plt.figure(figsize=(8.6, 4.1))
    gs = fig.add_gridspec(2, 2, width_ratios=[2.3, 1.0], height_ratios=[2.6, 1.0],
                          hspace=0.07, wspace=0.30)
    axT = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[1, 0], sharex=axT)
    axR = fig.add_subplot(gs[:, 1])
    x = np.arange(len(METALS))

    YT = (-0.20, 0.37)
    YB = (-0.45, -0.25)
    for ax in (axT, axB):
        ax.axhline(0, color="k", lw=0.9, zorder=2)
        for lab, S, c, mk in SERIES:
            ax.plot(x, S - TI, mk, color=c, ms=6, label=lab, zorder=3)
        ax.grid(axis="y", ls=":", lw=0.5, alpha=0.6)
    axT.axhspan(-0.05, 0.05, color="0.85", alpha=0.6, zorder=1, label=r"$\pm0.05\,k_B$")
    axT.set_ylim(*YT)
    axB.set_ylim(*YB)

    axT.spines["bottom"].set_visible(False)
    axB.spines["top"].set_visible(False)
    axT.tick_params(labelbottom=False, bottom=False)
    d = 0.012
    kw = dict(transform=axT.transAxes, color="k", clip_on=False, lw=0.9)
    axT.plot((-d, +d), (-d, +d), **kw)
    axT.plot((1 - d, 1 + d), (-d, +d), **kw)
    kw.update(transform=axB.transAxes)
    axB.plot((-d, +d), (1 - 3 * d, 1 + 3 * d), **kw)
    axB.plot((1 - d, 1 + d), (1 - 3 * d, 1 + 3 * d), **kw)

    axB.set_xticks(x)
    axB.set_xticklabels(METALS)
    axB.set_xlabel("Sutton–Chen metal (increasing $T$, $S_{\\rm TI}$)")
    axT.set_title("Deviation from thermodynamic integration (Sun 2017)", fontsize=9.5)
    fig.supylabel("")
    axT.set_ylabel(r"signed $\Delta S = S_{\rm method}-S_{\rm TI}$  ($k_B$/atom)")
    axT.yaxis.set_label_coords(-0.10, 0.30)
    axT.legend(fontsize=7.5, loc="upper right", framealpha=0.9, ncol=2,
               columnspacing=1.0, handletextpad=0.3)

    for i, (lab, S, c, mk) in enumerate(SERIES):
        r = rms(S - TI)
        axR.bar(i, r, color=c, width=0.7, zorder=3)
        axR.text(i, r + 0.006, f"{r:.2f}", ha="center", va="bottom",
                 fontsize=8.5, fontweight="bold")
    axR.set_xticks(range(len(SERIES)))
    axR.set_xticklabels(["rig", "+lnZ", "Desj", "R2PT", "3PT"], fontsize=8, rotation=0)
    axR.set_ylabel(r"RMS deviation vs TI  ($k_B$/atom)")
    axR.set_title("RMS", fontsize=9.5)
    axR.set_ylim(0, 0.40)
    axR.grid(axis="y", ls=":", lw=0.5, alpha=0.6)

    fig.savefig(os.path.join(OUT, "fig_metals_accuracy.pdf"))
    plt.close(fig)


if __name__ == "__main__":
    main()
