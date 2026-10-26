"""Draws Figure S11 from data/simulation/opes_pmf."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from paths import data, OUT
from pmf_panels import load

COL = {"LiTPFA": "#0072B2", "LiPF6": "#CC79A7", "LiBF4": "#E69F00"}
LAB = {"LiTPFA": "LiTPFA", "LiPF6": "LiPF$_6$", "LiBF4": "LiBF$_4$"}
INK, INK2, MUTE, GRID = "#0b0b0b", "#52514e", "#8a8985", "#e6e5e1"
SALTS = ("LiTPFA", "LiPF6", "LiBF4")
PANELS = [
    ("FM", "FM (fluoromethane)", "anchored Li–solvent row",
     {"LiTPFA": ("litpfa_fm.npz", "pm_opes_analytic"), "LiPF6": ("lipf6_fm.npz", None),
      "LiBF4": ("libf4_fm.npz", None)}),
    ("ACN", "ACN (acetonitrile)", "anchored Li–solvent row",
     {"LiTPFA": ("litpfa_acn.npz", None), "LiPF6": ("lipf6_acn.npz", None),
      "LiBF4": ("libf4_acn.npz", None)}),
    ("Me2O", "Me$_2$O (dimethyl ether)", "geometric Li–solvent rows",
     {"LiTPFA": ("litpfa_me2o.npz", None), "LiPF6": ("lipf6_me2o.npz", None),
      "LiBF4": ("libf4_me2o.npz", None)}),
    ("DME", "DME (dimethoxyethane)", "geometric Li–solvent rows",
     {"LiTPFA": ("litpfa_dme.npz", None), "LiPF6": ("lipf6_dme.npz", None),
      "LiBF4": ("libf4_dme.npz", None)}),
    ("DEE", "DEE (diethyl ether)", "geometric Li–solvent rows",
     {"LiTPFA": ("litpfa_dee.npz", None), "LiPF6": ("lipf6_dee.npz", None)}),
]
NOTE = {"DEE": "no LiBF$_4$ run (insoluble in DEE)"}
WALL = {"LiTPFA": 2.0, "LiPF6": 2.0, "LiBF4": 1.4}
DIRECT = {"LiTPFA": (4.4, 17.0), "LiPF6": (3.35, -27.0), "LiBF4": (2.25, -52.0)}


def main():
    panels = PANELS
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.5, "axes.linewidth": 0.9,
                         "axes.edgecolor": INK, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
                         "figure.facecolor": "white", "savefig.facecolor": "white"})
    fig = plt.figure(figsize=(7.2, 4.9), dpi=300)
    L, R, B, T, gx, gy = 0.075, 0.985, 0.085, 0.875, 0.075, 0.13
    wp, hp = (R - L - 2 * gx) / 3, (T - B - gy) / 2
    for idx in range(6):
        row, col = divmod(idx, 3)
        box = [L + col * (wp + gx), B + (1 - row) * (hp + gy), wp, hp]
        if idx >= len(panels):
            ax = fig.add_axes(box); ax.axis("off")
            continue
        key, title, footing, src = panels[idx]
        ax = fig.add_axes(box)
        ax.axvspan(10, 12, color=GRID, alpha=0.8, lw=0, zorder=0)
        ax.axhline(0, color=MUTE, lw=0.8, ls=":", zorder=1)
        ax.grid(True, axis="y", color=GRID, lw=0.5, zorder=1)
        lo, hi = [], []
        for salt in SALTS:
            if salt not in src:
                continue
            f, pre = src[salt]
            cen, m, s, cov = load(data("simulation", "opes_pmf", f), pre)
            k = np.isfinite(m) & np.isfinite(s) & (cen <= 12.05) & (cen >= WALL[salt])
            ax.fill_between(cen[k], (m - s)[k], (m + s)[k], color=COL[salt], alpha=0.22, lw=0, zorder=3)
            ax.plot(cen[k], m[k], color=COL[salt], lw=1.6, zorder=5)
            lo.append(np.nanmin((m - s)[k])); hi.append(np.nanmax((m + s)[k]))
        ylo = 10 * np.floor((min(lo) - 3) / 10)
        yhi = 10 * np.ceil((max(hi) + 3) / 10)
        ax.set_ylim(ylo, yhi); ax.set_xlim(1.2, 12.2)
        ax.set_xticks([2, 4, 6, 8, 10, 12])
        if row == 1 or idx + 3 >= len(panels):
            ax.set_xlabel("Li⁺–anion distance r (Å)", fontsize=8.0)
        if col == 0:
            ax.set_ylabel("W(r) (k$_B$T)", fontsize=8.0)
        letter = "abcdef"[idx]
        ax.set_title(f"({letter})  {title}", loc="left", fontsize=8.8, color=INK, pad=13)
        if key in NOTE:
            ax.text(0.97, 0.06, NOTE[key], transform=ax.transAxes, ha="right", va="bottom", fontsize=7.0, color=INK2)
        if idx == 0:
            for salt, (rx, wy) in DIRECT.items():
                ax.text(rx, wy, LAB[salt], color=INK, fontsize=7.6, fontweight="bold", ha="left", va="center", zorder=7,
                        bbox=dict(boxstyle="round,pad=0.15", fc="white", ec=COL[salt], lw=0.9))
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    x = L
    for salt in SALTS:
        fig.patches.append(plt.Rectangle((x, 0.962), 0.022, 0.016, transform=fig.transFigure, color=COL[salt], lw=0))
        fig.text(x + 0.028, 0.970, LAB[salt], fontsize=8.6, color=INK, va="center")
        x += 0.13
    fig.savefig(OUT / "figS11.png", dpi=300)
    print("figS11.png written")


if __name__ == "__main__":
    main()
