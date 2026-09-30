#!/usr/bin/env python3
"""Figure 1: one functional, two limits (Secs. 2 and 3).

The paper's claim is that one charge-equilibration solve is well posed for molecules AND for conductors. Top: the
functional and the constrained solve. Left: the molecular limit, one neutrality constraint per molecule; the
conventional failure is charge flowing between molecules. Right: the conductor limit, the finite-width
self-energy on the diagonal; the conventional failure is an indefinite operator as the hardness goes to zero.
Bottom: both in one cell and one projected Krylov solve, with the remaining run-time settings as a strip.

Schematic; nothing here is measured, so no number appears (Figure 3 carries the measurement).
The previous, settings-centred version is kept as fig_architecture_axes.py.bak_20260927.
"""
import os

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, FancyArrowPatch, FancyBboxPatch

import figstyle
from figstyle import BLUE, INK, MUTED, VERM

HERE = os.path.dirname(os.path.abspath(__file__))
figstyle.apply()

GOLD = "#b8860b"
BAD = "#9b1c1c"
fig, ax = plt.subplots(figsize=(figstyle.TEXTWIDTH, 3.55))
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.axis("off")


def box(x0, y0, x1, y1, edge=MUTED, face="none", lw=0.7, ls="-", z=1):
    ax.add_patch(FancyBboxPatch((x0, y0), x1 - x0, y1 - y0, boxstyle="round,pad=0,rounding_size=1.6",
                                linewidth=lw, linestyle=ls, edgecolor=edge, facecolor=face, zorder=z))


def arrow(p0, p1, colour=MUTED, rad=0.0, lw=0.8):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=8, lw=lw, color=colour,
                                 shrinkA=0, shrinkB=0, zorder=4, connectionstyle=f"arc3,rad={rad}"))


XIN, YIN = 6.5 / 100, 3.55 / 100      # inches per data unit in x and y


def disk(x, y, r, **kw):
    """a true circle of radius r INCHES at data (x, y)"""
    ax.add_patch(Ellipse((x, y), 2 * r / XIN, 2 * r / YIN, **kw))


def water(x, y, flip=1):
    """O with two H at 104.5 degrees, drawn to scale in inches"""
    import math
    disk(x, y, 0.075, fc="#e8eef5", ec=BLUE, lw=0.6, zorder=3)
    for sgn in (-1, 1):
        a = math.radians(-90 + sgn * 52.25) if flip == 1 else math.radians(90 + sgn * 52.25)
        disk(x + 0.095 * math.cos(a) / XIN, y + 0.095 * math.sin(a) / YIN, 0.045, fc="white", ec=BLUE, lw=0.5,
             zorder=3)


def cross(x, y, colour=BAD):
    ax.text(x, y, "×", ha="center", va="center", fontsize=9, color=colour, fontweight="bold", zorder=5)


# ---------------------------------------------------------------- the functional and the solve
box(4, 84.5, 96, 100, face="#f2f6fa", edge=INK, lw=0.8)
ax.text(50, 94.6, r"$E[\mathbf{q}]=\sum_i\left(\chi_i q_i+\frac{1}{2}\eta_i q_i^2\right)"
        r"+\frac{1}{2}\sum_{i\neq j} q_i\,J_{ij}(r_{ij})\,q_j$", ha="center", va="center", fontsize=7.6)
ax.text(50, 87.6, r"one projected Krylov solve,  $\mathbf{P}\mathbf{H}\mathbf{P}\,\boldsymbol{\Delta}"
        r"=\mathbf{P}\mathbf{b}$,  reciprocal space inside the matrix–vector product",
        ha="center", va="center", fontsize=6.0, color=MUTED)

# ---------------------------------------------------------------- molecular limit (left)
L0, L1 = 1, 48.5
box(L0, 27, L1, 82, edge=BLUE, lw=0.8)
ax.text(L0 + 2, 79.3, "molecular limit", ha="left", va="center", fontsize=7.6, color=BLUE, fontweight="bold")
ax.text(L1 - 2, 79.3, "insulators, liquids", ha="right", va="center", fontsize=6.3, color=MUTED)
# three molecules, each in its own fragment box
for k, x in enumerate((9.5, 23.5, 37.5)):
    box(x - 5.2, 58.5, x + 5.2, 73.5, edge=BLUE, lw=0.5, ls=(0, (2, 1.5)), z=2)
    water(x, 67.5, flip=1 if k != 1 else -1)
    ax.text(x, 61.0, r"$\sum q=0$", ha="center", va="center", fontsize=6.0, color=BLUE)
ax.text(24.75, 53.5, r"one neutrality constraint per molecule, $\mathbf{P}$",
        ha="center", va="center", fontsize=6.6, color=INK)
ax.text(24.75, 48.8, "charge moves within a molecule, not between them",
        ha="center", va="center", fontsize=6.3, color=INK)
cross(4.5, 38.0)
ax.text(7.2, 38.0, "conventional QEq, one global constraint:", ha="left", va="center", fontsize=6.1,
        color=BAD)
ax.text(7.2, 33.6, "charge flows between molecules, as in a metal", ha="left", va="center", fontsize=6.1, color=BAD)

# ---------------------------------------------------------------- conductor limit (right)
R0, R1 = 51.5, 99
box(R0, 27, R1, 82, edge=GOLD, lw=0.8)
ax.text(R0 + 2, 79.3, "conductor limit", ha="left", va="center", fontsize=7.6, color=GOLD, fontweight="bold")
ax.text(R1 - 2, 79.3, r"metals, $\eta\to 0$", ha="right", va="center", fontsize=6.3, color=MUTED)
# a slab: three rows of gold atoms, each wrapped in its Gaussian (drawn as a faint halo)
for row, y in enumerate((72.0, 67.2, 62.4)):
    for j in range(9):
        x = 57.0 + j * 4.6 + (2.3 if row % 2 else 0)
        if x > 95: continue
        disk(x, y, 0.11, fc="#fbf1d6", ec="none", zorder=2)
        disk(x, y, 0.055, fc="#f3d57a", ec=GOLD, lw=0.5, zorder=3)
ax.text(75.25, 55.0, r"finite-width self-energy on the diagonal:",
        ha="center", va="center", fontsize=6.6, color=INK)
ax.text(75.25, 50.3, r"$\eta_i\;\to\;\eta_i+K_e\,(2\alpha_w/\pi)^{1/2}$",
        ha="center", va="center", fontsize=7.0, color=INK)
ax.text(75.25, 45.6, r"$\mathbf{H}$ positive definite for any $\eta\geq 0$",
        ha="center", va="center", fontsize=6.3, color=INK)
cross(55.0, 38.0)
ax.text(57.7, 38.0, "conventional QEq, point charges:", ha="left", va="center", fontsize=6.1,
        color=BAD)
ax.text(57.7, 33.6, r"the operator turns indefinite as $\eta\to 0$", ha="left", va="center", fontsize=6.1, color=BAD)

# ---------------------------------------------------------------- both, in one cell and one solve
box(18, 12.5, 82, 22.5, face="#f2f6fa", edge=INK, lw=0.8)
ax.text(50, 19.6, "both limits in one cell, one solve", ha="center", va="center", fontsize=7.4, color=INK,
        fontweight="bold")
ax.text(50, 15.3, "molecules keep their constraints; the conductor keeps its self-energy",
        ha="center", va="center", fontsize=6.3, color=INK)
arrow((24.75, 27.0), (30.0, 22.8), colour=BLUE, rad=0.15)
arrow((75.25, 27.0), (70.0, 22.8), colour=GOLD, rad=-0.15)

# ---------------------------------------------------------------- the remaining settings
ax.text(50, 6.8, "run-time settings:   shielding kernel  cbrt · gaussian · slater      "
        "long range  mesh, inside the matrix–vector product", ha="center", va="center", fontsize=5.9, color=MUTED)
ax.text(50, 2.4, "charge dynamics  solve each step · predictor–corrector · extended Lagrangian",
        ha="center", va="center", fontsize=5.9, color=MUTED)

for ext in ("pdf", "png"):
    fig.savefig(os.path.join(HERE, f"fig_architecture.{ext}"), dpi=300)
print("wrote fig_architecture.pdf / .png")
