"""Defines the plot style used by Figure S3."""
from __future__ import annotations

import matplotlib as mpl
import matplotlib.pyplot as plt

SALT_ORDER = ["LiBF4", "LiPF6", "LiFSI", "LiTFSI", "LiTPFA"]

SALT_LABEL = {
    "LiBF4":  r"LiBF$_4$",
    "LiPF6":  r"LiPF$_6$",
    "LiFSI":  "LiFSI",
    "LiTFSI": "LiTFSI",
    "LiTPFA": "LiTPFA",
}

SOLVENT_ORDER = ["ACN", "DME", "DEE", "Me2O", "FM"]

SOLVENT_LABEL = {
    "ACN":  "ACN",
    "DME":  "DME",
    "DEE":  "DEE",
    "Me2O": r"Me$_2$O",
    "FM":   "FM",
}

SALT_COLOR = {
    "LiBF4":  "#eda100",
    "LiPF6":  "#e87ba4",
    "LiFSI":  "#2a78d6",
    "LiTFSI": "#008300",
    "LiTPFA": "#4a3aa7",
}

SALT_MARKER = {
    "LiBF4":  "s",
    "LiPF6":  "o",
    "LiFSI":  "^",
    "LiTFSI": "v",
    "LiTPFA": "D",
}

SOLVENT_MARKER = {"ACN": "o", "DME": "s", "DEE": "^", "Me2O": "v", "FM": "D"}


CONC_RAMP = ["#b0aad0", "#8880b4", "#5d5296", "#33257d"]


CONC_MARKER = ["s", "o", "D", "^"]

MARKER_SIZE_TRIM = {"s": 0.90, "o": 1.00, "D": 1.05, "^": 1.10, "v": 1.10}


def conc_markers(n: int):
    if not 1 <= n <= len(CONC_MARKER):
        raise ValueError(f"conc_markers: n must be 1..{len(CONC_MARKER)}, got {n}")
    return list(CONC_MARKER[:n])


def conc_colors(n: int):
    if not 1 <= n <= len(CONC_RAMP):
        raise ValueError(f"conc_colors: n must be 1..{len(CONC_RAMP)}, got {n}")
    if n == len(CONC_RAMP):
        return list(CONC_RAMP)
    step = (len(CONC_RAMP) - 1) / (n - 1) if n > 1 else 0
    return [CONC_RAMP[round(i * step)] for i in range(n)]


INK        = "#1a1a19"
INK_2      = "#55554e"
GRID       = "#dcdcd6"
MARKER_EDGE = "#2b2b28"

MARKER_SIZE      = 6.0
MARKER_EDGE_W    = 0.6
LINE_W           = 1.6
ERRORBAR_W       = 1.0
ERRORBAR_CAPSIZE = 2.6

FS_TICK   = 9.5
FS_LABEL  = 10.5
FS_LEGEND = 9
FS_PANEL  = 12.5

W_DOUBLE = 7.087
W_SINGLE = 3.465


def use_style():
    mpl.rcParams.update({
        "figure.dpi":          150,
        "savefig.dpi":         600,
        "savefig.bbox":        "tight",
        "savefig.pad_inches":  0.02,
        "pdf.fonttype":        42,
        "ps.fonttype":         42,
        "font.family":         "sans-serif",
        "font.sans-serif":     ["DejaVu Sans", "Arial", "Helvetica"],
        "font.size":           FS_TICK,
        "axes.labelsize":      FS_LABEL,
        "axes.labelcolor":     INK,
        "axes.edgecolor":      INK_2,
        "axes.linewidth":      0.8,
        "axes.spines.top":     False,
        "axes.spines.right":   False,
        "axes.grid":           False,
        "grid.color":          GRID,
        "grid.linewidth":      0.6,
        "xtick.color":         INK_2,
        "ytick.color":         INK_2,
        "xtick.labelsize":     FS_TICK,
        "ytick.labelsize":     FS_TICK,
        "xtick.direction":     "out",
        "ytick.direction":     "out",
        "xtick.major.size":    3.0,
        "ytick.major.size":    3.0,
        "xtick.major.width":   0.8,
        "ytick.major.width":   0.8,
        "legend.frameon":      False,
        "legend.fontsize":     FS_LEGEND,
        "legend.handletextpad": 0.5,
        "legend.labelspacing":  0.35,
        "legend.borderpad":     0.2,
        "lines.linewidth":     LINE_W,
        "lines.markersize":    MARKER_SIZE,
    })


def style_axes(ax, ygrid: bool = True):
    if ygrid:
        ax.set_axisbelow(True)
        ax.yaxis.grid(True, color=GRID, linewidth=0.6)
    ax.tick_params(labelcolor=INK_2)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK_2)


def panel_label(ax, letter: str, dx: float = -0.16, dy: float = 1.06):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=FS_PANEL,
            fontweight="bold", color=INK, va="top", ha="left")


def salt_legend(ax, salts, **kw):
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], linestyle="none", marker=SALT_MARKER[s],
                      markerfacecolor=SALT_COLOR[s], markeredgecolor=MARKER_EDGE,
                      markeredgewidth=MARKER_EDGE_W, markersize=MARKER_SIZE,
                      label=SALT_LABEL[s])
               for s in salts]
    return ax.legend(handles=handles, **kw)


def save(fig, stem: str, outdir=None, tight: bool = True):
    import pathlib
    outdir = pathlib.Path(outdir or pathlib.Path(__file__).parent / "out")
    outdir.mkdir(parents=True, exist_ok=True)
    paths = []
    kw = {} if tight else dict(bbox_inches=fig.bbox_inches, pad_inches=0.0)
    for ext in ("pdf", "png"):
        p = outdir / f"{stem}.{ext}"
        fig.savefig(p, **kw)
        paths.append(p)
    return paths
