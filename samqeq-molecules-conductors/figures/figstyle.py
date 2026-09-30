"""Shared matplotlib style for the figures: Elsevier cas-sc text width, Helvetica-type sans-serif text, a full frame with inward ticks, Okabe-Ito colours, markers and dash patterns."""
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
        "font.family": "sans-serif",
        "font.sans-serif": ["Liberation Sans", "Nimbus Sans", "DejaVu Sans"],
        "mathtext.fontset": "custom",
        "mathtext.rm": "Liberation Sans",
        "mathtext.it": "Liberation Sans:italic",
        "mathtext.bf": "Liberation Sans:bold",
        "mathtext.sf": "Liberation Sans",
        "mathtext.cal": "Liberation Sans",
        "font.size": 8,
        "axes.labelsize": 8,
        "axes.titlesize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "axes.edgecolor": "black",
        "axes.labelcolor": "black",
        "text.color": INK,
        "xtick.color": "black",
        "ytick.color": "black",
        "axes.linewidth": 0.6,
        "axes.grid": False,
        "axes.spines.top": True,
        "axes.spines.right": True,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.major.size": 3.5,
        "ytick.major.size": 3.5,
        "xtick.minor.size": 1.8,
        "ytick.minor.size": 1.8,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.minor.width": 0.4,
        "ytick.minor.width": 0.4,
        "legend.frameon": False,
        "legend.handlelength": 2.0,
        "lines.linewidth": 1.0,
        "lines.markersize": 4.0,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
        "pdf.fonttype": 42,
    })

def panel_label(ax, letter, dx=-0.20, dy=1.04):
    ax.text(dx, dy, f"({letter})", transform=ax.transAxes,
            fontsize=9, fontweight="bold", va="bottom", ha="left")

def recessive(ax, which="both"):
    ax.grid(False)
