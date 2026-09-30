"""Shared matplotlib style for the figures: Elsevier cas-sc text width, Okabe-Ito colours, markers and dash patterns."""
import matplotlib as mpl

BLUE = "#0072B2"
GREEN = "#009E73"
VERM = "#D55E00"
PURPLE = "#CC79A7"
INK = "#1a1a1a"
MUTED = "#6b6b6b"
GRID = "#d8d8d8"

TEXTWIDTH = 6.5

def apply():
    mpl.rcParams.update({
        "font.family": "serif",
        "font.serif": ["DejaVu Serif"],
        "mathtext.fontset": "dejavuserif",
        "font.size": 8,
        "axes.labelsize": 8,
        "axes.titlesize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "axes.edgecolor": MUTED,
        "axes.labelcolor": INK,
        "text.color": INK,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "xtick.labelcolor": INK,
        "ytick.labelcolor": INK,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.minor.width": 0.4,
        "ytick.minor.width": 0.4,
        "grid.color": GRID,
        "grid.linewidth": 0.5,
        "legend.frameon": False,
        "legend.handlelength": 2.0,
        "lines.linewidth": 1.2,
        "lines.markersize": 4.0,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
        "pdf.fonttype": 42,
    })

def panel_label(ax, letter, dx=-0.20, dy=1.04):
    ax.text(dx, dy, f"({letter})", transform=ax.transAxes,
            fontsize=9, fontweight="bold", va="bottom", ha="left")

def recessive(ax, which="both"):
    ax.grid(True, which=which, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
