"""Draws panels A to E of Figure 3 from data/simulation."""
from __future__ import annotations

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

import fig3_layout as F5
from paths import OUT

SRC = F5.SRC
INK, INK2, MUTE, GRID = F5.INK, F5.INK2, F5.MUTE, F5.GRID
FS_LABEL, FS_TICK, FS_NOTE = F5.FS_LABEL, F5.FS_TICK, F5.FS_NOTE
FS_TITLE = 8.1
TITLES = {"a": "Li⁺–solvent binding energies",
          "b": "Li⁺ solvation structure in FM",
          "c": "MD simulation cell snapshots",
          "d": "Ion aggregation distribution in FM",
          "e": "Ion pairing thermodynamics in FM",
          "f": "Mixing free energy"}
def panel_d(ax, report):
    ns = SRC.nspecies_canonical()
    hc, rc = f"{SRC.NSPECIES_CUT:.1f}", f"{SRC.NSPECIES_ROBUST_CUT:.1f}"
    W = 0.74
    for salt in ("LiPF6", "LiTPFA"):
        d = ns[salt]
        col = SRC.SALT_COL[salt]
        for n, fr in d["pn"][hc].items():
            ax.bar(n, fr, width=W, color=col, alpha=0.92, zorder=3, lw=0)
        for n, fr in d["pn"][rc].items():
            if fr < 0.004:
                continue
            ax.add_patch(Rectangle((n - W / 2, 0), W, fr, facecolor="none", edgecolor=col,
                                   lw=1.0, ls=(0, (2.6, 1.6)), zorder=5))
        for n in d["pn"][hc]:
            vals = [d["pn_replicas"][hc][r].get(n, 0.0) for r in d["pn_replicas"][hc]]
            ax.plot([n] * len(vals), vals, marker="o", ms=2.8, mfc="white", mec=INK2, mew=0.7,
                    ls="none", zorder=6)
    med = {s: ns[s]["median_sweep"][hc] for s in ("LiPF6", "LiTPFA")}
    mode = max(ns["LiPF6"]["pn"][hc], key=ns["LiPF6"]["pn"][hc].get)
    ax.annotate(f"LiPF$_6$\nmedian {med['LiPF6']:.0f}", xy=(mode + 0.45, 0.72),
                xytext=(5.4, 0.80), ha="left", va="center", fontsize=FS_NOTE,
                color=SRC.SALT_COL["LiPF6"], fontweight="bold", linespacing=1.15, zorder=7,
                bbox=dict(facecolor="white", edgecolor="none", pad=0.9),
                arrowprops=dict(arrowstyle="-", color=SRC.SALT_COL["LiPF6"], lw=0.9,
                                shrinkA=1, shrinkB=1))
    ax.text(16.4, 1.035, f"LiTPFA\nmedian {med['LiTPFA']:.0f}", ha="right", va="bottom",
            fontsize=FS_NOTE, color=SRC.SALT_COL["LiTPFA"], fontweight="bold", linespacing=1.15,
            zorder=7)
    for x in (SRC.N_PAIRS, SRC.N_IONS):
        ax.axvline(x, color=INK2, lw=0.9, ls="--", zorder=2)
    ax.set_xlim(-0.6, 17.6)
    ax.set_ylim(0, 1.30)
    ax.set_xticks([0, 4, 8, 12, 16])
    ax.set_yticks([0, 0.5, 1.0])
    ax.set_xlabel("distinct ionic species $n$ (of 16)", fontsize=FS_LABEL)
    ax.set_ylabel("fraction of frames", fontsize=FS_LABEL)
    ax.tick_params(labelsize=FS_TICK)
    ax.grid(True, axis="y", color=GRID, lw=0.5, zorder=0)
    F5.despine(ax)


GUTTER = (0.72, 0.52, 0.62)
AXW, AXH = 1.68, 1.70
LPAD, GAP, RPAD = 0.06, 0.10, 0.04
AX_L = (LPAD + GUTTER[0],
        LPAD + GUTTER[0] + AXW + GAP + GUTTER[1],
        LPAD + GUTTER[0] + AXW + GAP + GUTTER[1] + AXW + GAP + GUTTER[2])
CELL_L = tuple(x - g for x, g in zip(AX_L, GUTTER))
CELL_R = (CELL_L[1], CELL_L[2], AX_L[2] + AXW + RPAD)
ROW_B = (3.32, 0.82)
HEAD_Y = (5.38, 2.82)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(OUT / "fig3_panels_ae"))
    ap.add_argument("--width", type=float, default=7.20)
    ap.add_argument("--height", type=float, default=5.52)
    a = ap.parse_args()
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.4, "axes.linewidth": 0.9,
                         "axes.edgecolor": INK, "axes.labelcolor": INK, "xtick.color": INK,
                         "ytick.color": INK, "figure.facecolor": "white",
                         "savefig.facecolor": "white"})
    fig = plt.figure(figsize=(a.width, a.height), dpi=300)
    titles = {}

    def X(inches):
        return inches / a.width

    def V(inches):
        return inches / a.height

    def head(letter, col, row):
        titles[letter] = fig.text(X(CELL_L[col]), V(HEAD_Y[row]),
                                  f"({letter})  {TITLES[letter]}", fontsize=FS_TITLE, color=INK,
                                  va="top")

    ax_a = fig.add_axes([X(AX_L[0]), V(ROW_B[0]), X(AXW), V(AXH)])
    ax_bs = fig.add_axes([X(CELL_L[1] + 0.06), V(ROW_B[0] + 0.04), X(1.66), V(AXH - 0.04)])
    ax_bk = fig.add_axes([X(CELL_L[1] + 1.78), V(ROW_B[0] + 0.20), X(0.34), V(AXH - 0.40)])
    F5.panel_a(ax_a, ax_bs, ax_bk, [])
    head("a", 0, 0)
    head("b", 1, 0)

    cw = CELL_R[2] - CELL_L[2] - RPAD
    cell = (cw - 0.08) / 2.0
    ax_ck = fig.add_axes([X(CELL_L[2]), V(ROW_B[0] + AXH - 0.30), X(cw), V(0.30)])
    ax_c1 = fig.add_axes([X(CELL_L[2]), V(ROW_B[0] + 0.18), X(cell), V(cell)])
    ax_c2 = fig.add_axes([X(CELL_L[2] + cell + 0.08), V(ROW_B[0] + 0.18), X(cell), V(cell)])
    F5.panel_b(ax_ck, ax_c1, ax_c2, [])
    head("c", 2, 0)

    ax_d = fig.add_axes([X(AX_L[0]), V(ROW_B[1]), X(AXW), V(AXH)])
    panel_d(ax_d, [])
    head("d", 0, 1)

    ax_e = fig.add_axes([X(AX_L[1]), V(ROW_B[1]), X(AXW), V(AXH)])
    F5.panel_d(ax_e, [])
    head("e", 1, 1)

    ax_f = fig.add_axes([X(AX_L[2]), V(ROW_B[1]), X(AXW), V(AXH)])
    F5.panel_e(ax_f, [])
    head("f", 2, 1)

    fig.canvas.draw()
    fig.savefig(a.out + ".png", dpi=300)
    print(f"{os.path.basename(a.out)}.png written")


if __name__ == "__main__":
    main()
