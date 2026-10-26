"""Draws Figure S19 from the phonon densities of states in data/simulation/phonons."""
import gzip, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from paths import data as datafile, OUT

THZ_CM = 33.35641


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def read_mesh(path):
    modes, w = [], None
    with _open(path) as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("weight:"):
                w = int(s.split(":")[1])
            elif s.startswith("frequency:"):
                modes.append((float(s.split(":")[1]), w))
    return modes


def read_natom(path):
    with _open(path) as fh:
        for line in fh:
            if line.strip().startswith("natom:"):
                return int(line.split(":")[1])


COL = {"LiBF4": "#E69F00", "LiPF6": "#CC79A7", "LiFSI": "#56B4E9", "LiTFSI": "#009E73", "LiTPFA": "#0072B2"}
LAB = {"LiBF4": "LiBF$_4$", "LiPF6": "LiPF$_6$", "LiFSI": "LiFSI", "LiTFSI": "LiTFSI", "LiTPFA": "LiTPFA"}
SRC = [("LiBF4", "libf4", 6), ("LiPF6", "lipf6", 8), ("LiFSI", "lifsi", 10),
       ("LiTFSI", "litfsi", 16), ("LiTPFA", "litpfa", 58)]
INK, INK2, MUTE, GRID = "#0b0b0b", "#52514e", "#8a8985", "#e6e5e1"
BIN = 5.0


def main():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.5, "axes.linewidth": 0.9, "axes.edgecolor": INK,
                         "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
                         "figure.facecolor": "white", "savefig.facecolor": "white"})
    data = []
    for salt, d, apf in SRC:
        natom = read_natom(datafile("simulation", "phonons", d, "thermal_properties.yaml"))
        nb = 3 * natom
        modes = read_mesh(datafile("simulation", "phonons", d, "mesh.yaml"))
        qw = sum(w for _, w in modes) / nb
        nfu = natom / apf
        nu = np.array([f for f, _ in modes]) * THZ_CM
        wt = np.array([w for _, w in modes], float) / qw / nfu
        pos = nu > 0
        edges = np.arange(0.0, nu[pos].max() + BIN, BIN)
        h, _ = np.histogram(nu[pos], bins=edges, weights=wt[pos])
        g = h / BIN
        imag = 100 * wt[~pos].sum() / wt.sum()
        total = h.sum()
        data.append((salt, 0.5 * (edges[1:] + edges[:-1]), g, total, apf, nu[pos].max()))
    xmax = 100 * np.ceil(max(x[5] for x in data) / 100) + 50
    fig, axes = plt.subplots(len(data), 1, figsize=(6.5, 6.9), dpi=300, sharex=True)
    fig.subplots_adjust(left=0.105, right=0.985, top=0.985, bottom=0.075, hspace=0.12)
    for ax, (salt, x, g, total, apf, _) in zip(axes, data):
        ax.fill_between(x, 0, g, color=COL[salt], alpha=0.25, lw=0, zorder=2)
        ax.plot(x, g, color=COL[salt], lw=1.0, zorder=3)
        top = g.max() * 1.45
        ax.set_ylim(0, top); ax.set_xlim(0, xmax)
        ax.grid(True, axis="x", color=GRID, lw=0.5, zorder=0)
        ax.tick_params(axis="y", labelsize=7.5)
        ax.plot([0.012, 0.042], [0.90, 0.90], transform=ax.transAxes, color=COL[salt], lw=3, solid_capstyle="butt")
        ax.text(0.052, 0.90, LAB[salt], transform=ax.transAxes, ha="left", va="center", fontsize=9.5, fontweight="bold", color=INK)
        ax.text(0.99, 0.90, f"3n = {3 * apf} modes per formula unit", transform=ax.transAxes, ha="right", va="center",
                fontsize=7.3, color=INK2)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    axes[-1].set_xlabel("ν (cm$^{-1}$)", fontsize=9)
    fig.text(0.018, 0.53, "g(ν) (modes cm$^{-1}$ per formula unit)", rotation=90, va="center", ha="center", fontsize=9, color=INK)
    fig.savefig(OUT / "figS19.png", dpi=300)
    print("figS19.png written")


if __name__ == "__main__":
    main()
