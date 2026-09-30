#!/usr/bin/env python3
"""Figure: where samQEq sits in a timestep, and what the operator is made of.

Companion to fig_architecture, which shows which model choices are run-time
settings. This one shows the code: the fix is called once per step before the
forces, the solve is a projected Krylov iteration whose only argument is a
matrix-vector product, and that product is assembled in one part and applied in
two. The host/device boundary is drawn where it actually falls.

Schematic. Nothing here is measured, so no number appears.
"""
import os

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

import figstyle
from figstyle import BLUE, GREEN, INK, MUTED, VERM

HERE = os.path.dirname(os.path.abspath(__file__))
figstyle.apply()

fig, ax = plt.subplots(figsize=(figstyle.TEXTWIDTH, 3.05))
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.axis("off")


def box(x, y, w, h, label, sub=None, fc="none", ec=MUTED, lw=0.9, ls="-",
        tc=INK, fs=8, z=3):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.6,rounding_size=1.6",
                                fc=fc, ec=ec, lw=lw, ls=ls, zorder=z))
    ax.text(x + w / 2, y + h / 2 + (1.9 if sub else 0), label, ha="center", va="center",
            fontsize=fs, color=tc, zorder=z + 1)
    if sub:
        ax.text(x + w / 2, y + h / 2 - 3.4, sub, ha="center", va="center",
                fontsize=6.6, color=MUTED, zorder=z + 1, style="italic")


def arrow(p, q, rad=0.0, color=MUTED, lw=1.0, ls="-"):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=9,
                                 lw=lw, ls=ls, color=color, zorder=5,
                                 connectionstyle=f"arc3,rad={rad}"))


# ----------------------------------------------------------------- timestep band
ax.add_patch(Rectangle((0.5, 68), 99, 30, fc="#f4f4f2", ec="none", zorder=0))
ax.text(1.8, 95.2, "O N E   L A M M P S   T I M E S T E P", fontsize=6.2, color=MUTED,
        ha="left", va="center", zorder=1)

box(2.0, 74, 19, 14, "neighbor list", sub="rebuilt as usual")
box(26.5, 74, 22, 14, "fix qeq/sam", sub="pre-force, every $N$ steps", ec=BLUE, lw=1.3)
box(54.0, 74, 20, 14, "pair + kspace", sub="forces from $q$")
box(79.5, 74, 18.5, 14, "integrate")

arrow((21.8, 81), (26.0, 81))
arrow((49.0, 81), (53.5, 81))
arrow((74.5, 81), (79.0, 81))
ax.text(51.2, 85.4, "$q$", fontsize=8, color=BLUE, ha="center", va="center", zorder=6)

# ----------------------------------------------------------------- the solve
ax.add_patch(FancyBboxPatch((20.0, 40), 35, 22,
                            boxstyle="round,pad=0.6,rounding_size=1.6",
                            fc="none", ec=BLUE, lw=1.3, zorder=3))
ax.text(37.5, 56.5, "projected Krylov solve", ha="center", va="center", fontsize=8,
        color=INK, zorder=4)
ax.text(37.5, 50.0, r"$\mathbf{x} \leftarrow \mathbf{P}\mathbf{x}$", ha="center", va="center",
        fontsize=8.5, color=INK, zorder=4)
ax.text(37.5, 44.4, r"$\mathbf{y} \leftarrow \mathbf{A}\mathbf{x}$", ha="center", va="center",
        fontsize=8.5, color=INK, zorder=4)
ax.text(57.0, 47.2, "per iterate", fontsize=6.6, color=MUTED, ha="left", va="center",
        style="italic")
arrow((32.0, 74), (32.0, 62.8), rad=0.0, color=BLUE)
arrow((43.0, 62.8), (43.0, 74), rad=0.0, color=BLUE)

# ----------------------------------------------------------------- the operator
ax.text(33.0, 33.0, r"$\mathbf{A} = \mathbf{P}\mathbf{H}\mathbf{P}$  is never formed. It is applied as:",
        fontsize=7.6, color=INK, ha="left", va="center")

def opbox(x, title, verb, note, ec):
    w, y, h = 29.0, 7.5, 18.5
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.6,rounding_size=1.6",
                                fc="none", ec=ec, lw=1.1, zorder=3))
    cx = x + w / 2
    ax.text(cx, y + 13.6, title, ha="center", va="center", fontsize=8, color=INK, zorder=4)
    ax.text(cx, y + 8.6, verb, ha="center", va="center", fontsize=6.6, color=ec,
            style="italic", zorder=4)
    ax.text(cx, y + 3.9, note, ha="center", va="center", fontsize=6.3, color=MUTED, zorder=4)


opbox(2.0, "short range", "assembled, CSR", "linear in pairs, not in $N^2$", VERM)
opbox(35.5, "particle mesh", "applied, not assembled", "no outer SCF loop", GREEN)
opbox(69.0, "projector", "applied, not assembled", "one sweep over fragments", GREEN)

ax.text(33.2, 16.7, "$\\oplus$", fontsize=10, color=MUTED, ha="center", va="center", zorder=4)
ax.text(66.7, 16.7, "$\\oplus$", fontsize=10, color=MUTED, ha="center", va="center", zorder=4)
arrow((25.0, 40), (25.0, 27.6), color=MUTED)

# ----------------------------------------------------------------- device boundary
ax.add_patch(Rectangle((0.8, 4.6), 98.0, 23.2, fc="none", ec=BLUE, lw=0.8,
                       ls=(0, (4, 2.6)), zorder=2))
ax.text(99.0, 1.6, "device path (Kokkos): the product, the kernels and the mesh; "
                   "with devsolve on, the whole iteration",
        fontsize=6.4, color=BLUE, ha="right", va="center")

for ext in ("pdf", "png"):
    fig.savefig(os.path.join(HERE, f"fig_callstructure.{ext}"), dpi=300, bbox_inches="tight")
print("wrote fig_callstructure.pdf / .png")
