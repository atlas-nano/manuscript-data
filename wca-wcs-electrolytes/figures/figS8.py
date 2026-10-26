"""Draws Figure S8 from the Li-F radial distribution functions and ion trajectory in data/simulation/conc_ladder."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import rdf_cn as M
from paths import OUT

INK, INK2, GRID = "#0b0b0b", "#52514e", "#e6e5e1"
FS_LABEL, FS_TICK, FS_NOTE = 8.2, 7.8, 6.9
CONC_RAMP = ["#b0aad0", "#8880b4", "#5d5296", "#33257d"]
AN_COL, FM_COL = "#D55E00", "#2a78d6"


def despine(ax):
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


def panel_a(ax):
    for (tag, conc), col in zip(M.RUNGS, CONC_RAMP):
        d = M.fm_cn(tag)
        m = d["r"] <= 4.6
        ax.plot(d["r"][m], d["g"][m], color=col, lw=1.5, label=f"{conc:.3f} M", zorder=3)
    ax.set_ylabel("$g(r)$", fontsize=FS_LABEL)
    ax.set_xlim(1.2, 4.6)
    ax.set_ylim(0, 17)
    ax.set_yticks([5, 10, 15])
    ax.tick_params(labelsize=FS_TICK)
    ax.legend(loc="upper right", fontsize=FS_NOTE, frameon=False, handlelength=1.1, labelspacing=0.28,
              borderpad=0.1, handletextpad=0.5, title="LiTPFA/FM, 293 K", title_fontsize=FS_NOTE)
    ax.grid(True, axis="y", color=GRID, lw=0.5, zorder=0)
    despine(ax)


def panel_b(ax):
    n8 = M.fm_cn("litpfa_fm_n8")
    prod = [f for f in M.load_traj("litpfa_fm_n8") if f[0] > 50000]
    hist, nframes, n_li, n_an, V, _ = M.hist_over_frames(prod)
    _, an_n = M.gr_nr_from_hist(hist, nframes, n_li, n_an, V)
    _, nframes_full, _, _, _, closest_full = M.hist_over_frames(M.load_traj("litpfa_fm_n8"))
    zero_radius = next((M.EDGES[k] for k in range(M.NBINS) if hist[k] > 0), M.RMAX)
    m = n8["r"] <= 6.0
    am = M.CENTERS <= 6.0
    ax.axvspan(1.2, zero_radius, color=AN_COL, alpha=0.10, lw=0, zorder=0)
    ax.axvline(zero_radius, color=AN_COL, lw=0.8, ls=(0, (4, 2)), alpha=0.75, zorder=1)
    ax.plot(n8["r"][m], n8["n"][m], color=FM_COL, lw=1.6, zorder=3, label="Li–F(FM), solvent")
    ax.plot([n8["r_half"]], [n8["cn"]], marker="o", color=FM_COL, mec=INK, mew=0.5, ms=4.0, ls="none", zorder=4)
    ax.plot(M.CENTERS[am], an_n[am], color=AN_COL, lw=1.5, ls="--", zorder=3, label="Li–F(TPFA⁻), anion")
    ax.axhspan(n8["cn"] - M.CN_SIGMA, n8["cn"] + M.CN_SIGMA, color=FM_COL, alpha=0.14, zorder=1)
    ax.text(5.95, 4.95, f"N = {n8['cn']:.1f} ± {M.CN_SIGMA} at $r_{{1/2}}$ = {n8['r_half']:.2f} Å",
            ha="right", va="bottom", fontsize=FS_NOTE, color=FM_COL, fontweight="bold", zorder=5,
            bbox=dict(facecolor="white", edgecolor="none", pad=0.9))
    ax.text(zero_radius - 0.1, 9.3, f"no TPFA⁻ inside {zero_radius:.2f} Å", ha="right", va="bottom",
            fontsize=FS_NOTE, color=AN_COL, zorder=5)
    ax.text(0.982, 0.10, f"closest Li–F(TPFA⁻) approach\nover {nframes_full} frames: {closest_full:.2f} Å",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=FS_NOTE, color=AN_COL, linespacing=1.25,
            zorder=5, bbox=dict(facecolor="white", edgecolor="none", pad=0.8))
    ax.set_xlim(1.2, 6.0)
    ax.set_ylim(0, 17.5)
    ax.set_yticks([0, 5, 10, 15])
    ax.set_xlabel("Li–F separation $r$ (Å)", fontsize=FS_LABEL)
    ax.set_ylabel("$N(r)$, count per Li$^+$", fontsize=FS_LABEL)
    ax.tick_params(labelsize=FS_TICK)
    ax.grid(True, axis="y", color=GRID, lw=0.5, zorder=0)
    ax.legend(loc="upper left", fontsize=FS_NOTE, frameon=False, handlelength=1.5, labelspacing=0.28,
              borderpad=0.1, handletextpad=0.5, title="LiTPFA/FM, 0.232 M", title_fontsize=FS_NOTE)
    despine(ax)


def main():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.4, "axes.linewidth": 0.9,
                         "axes.edgecolor": INK, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
                         "figure.facecolor": "white", "savefig.facecolor": "white"})
    fig = plt.figure(figsize=(6.40, 2.566), dpi=300)
    panel_a(fig.add_axes([0.086, 0.421, 0.384, 0.538]))
    panel_b(fig.add_axes([0.591, 0.160, 0.384, 0.799]))
    fig.text(0.012, 0.93, "(a)", fontsize=9.2, color=INK, va="top")
    fig.text(0.515, 0.93, "(b)", fontsize=9.2, color=INK, va="top")
    fig.savefig(OUT / "figS8.png", dpi=300)
    print("figS8.png written")


if __name__ == "__main__":
    main()
