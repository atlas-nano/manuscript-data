"""Defines the shared style, palettes, paths and helpers of the figure scripts."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

from paths import ROOT, DATA, EXP, SIM, TABLES, OUT, data, panel_image, register_fonts

register_fonts()
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans", "Arimo"],
    "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
    "font.size": 7, "axes.labelsize": 7, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.minor.width": 0.4,
    "ytick.minor.width": 0.4, "xtick.direction": "in", "ytick.direction": "in", "xtick.major.size": 2.5,
    "ytick.major.size": 2.5, "xtick.minor.size": 1.5, "ytick.minor.size": 1.5, "lines.linewidth": 1.0,
    "legend.frameon": False, "savefig.dpi": 300, "pdf.fonttype": 42, "svg.fonttype": "none",
    "axes.unicode_minus": True,
})
CM = 1 / 2.54
W_FULL, W_MID, W_COL = 18.4 * CM, 12.1 * CM, 5.7 * CM
DPI = 300

SOLV = {"ACN": "#0072B2", "DME": "#AA3377", "Me2O": "#009E73", "DEE": "#E69F00"}
FM_FACE, FM_EDGE = "#F5C518", "#000000"
NEUTRAL = "#4D4D4D"
ELEC = {"WCA-WCS": "#d92120", "IL-LGE": "#3f51a3", "LHCE": "#e69f00"}
ELEC_MARK = {"WCA-WCS": "o", "IL-LGE": "s", "LHCE": "^"}
SALT_MARK = {"LiBF4": "o", "LiPF6": "s", "LiFSI": "^", "LiTFSI": "D", "LiTPFA": "*"}
SALT_TEX = {"LiBF4": r"LiBF$_4$", "LiPF6": r"LiPF$_6$", "LiFSI": "LiFSI", "LiTFSI": "LiTFSI", "LiTPFA": "LiTPFA"}
SOLV_TEX = {"ACN": "ACN", "DME": "DME", "Me2O": r"Me$_2$O", "DEE": "DEE", "FM": "FM"}
SALTS = ["LiBF4", "LiPF6", "LiFSI", "LiTFSI", "LiTPFA"]
SOLVENTS = ["ACN", "DME", "Me2O", "DEE"]

TABLE_S2 = {
    ("ACN", "LiBF4"): (1.58, 1.91, 0.7939, 0.2319, 0.5619, 1.348, 0.2275),
    ("ACN", "LiPF6"): (1.59, 2.16, 2.369, 0.8997, 1.47, 0.6903, 0.6318),
    ("ACN", "LiFSI"): (1.53, 2.18, 2.659, 1.004, 1.654, 0.5257, 0.7166),
    ("ACN", "LiTFSI"): (1.46, 1.79, 1.675, 0.6727, 1.003, 0.7873, 0.5155),
    ("ACN", "LiTPFA"): (1.42, 1.24, 1.889, 1.034, 0.8543, 0.3857, 0.7100),
    ("DME", "LiBF4"): (1.28, 1.27, 0.03676, 0.02338, 0.01338, 1.257, 0.01442),
    ("DME", "LiPF6"): (0.96, 1.00, 0.2723, 0.1161, 0.1561, 0.8439, 0.1389),
    ("DME", "LiFSI"): (0.95, 1.07, 0.5285, 0.2043, 0.3243, 0.7457, 0.2617),
    ("DME", "LiTFSI"): (0.92, 0.98, 0.3479, 0.144, 0.204, 0.776, 0.1831),
    ("DME", "LiTPFA"): (0.89, 0.79, 0.5221, 0.3111, 0.2111, 0.5789, 0.3108),
    ("DEE", "LiPF6"): (1.74, 1.75, 0.001513, -0.004243, 0.005757, 1.744, 4.336e-4),
    ("DEE", "LiFSI"): (1.89, 1.95, 0.002568, -0.02872, 0.03128, 1.919, 6.688e-4),
    ("DEE", "LiTFSI"): (1.99, 1.98, 0.001146, 0.005573, -0.004427, 1.984, 2.885e-4),
    ("DEE", "LiTPFA"): (0.99, 0.92, 0.358, 0.214, 0.144, 0.776, 0.1875),
    ("Me2O", "LiBF4"): (3.42, 3.48, 0.0714, 0.005698, 0.0657, 3.414, 0.01035),
    ("Me2O", "LiPF6"): (3.43, 3.43, 0.3527, 0.1764, 0.1764, 3.254, 0.05142),
    ("Me2O", "LiFSI"): (3.17, 3.21, 0.5648, 0.2624, 0.3024, 2.908, 0.08852),
    ("Me2O", "LiTFSI"): (3.60, 3.55, 0.4305, 0.2403, 0.1903, 3.36, 0.06021),
    ("Me2O", "LiTPFA"): (2.35, 1.99, 1.739, 1.05, 0.6895, 1.3, 0.4007),
    ("FM", "LiTPFA"): (5.37, 4.32, 2.556, 1.803, 0.7532, 3.567, 0.2638),
}
TABLE_S4 = {
    "LiBF4": {"ACN": 0.5992, "DME": 0.6332, "Me2O": 0.6145, "DEE": 0.5735, "FM": 0.7165},
    "LiPF6": {"ACN": 0.6749, "DME": 0.7133, "Me2O": 0.6922, "DEE": 0.6460, "FM": 0.8071},
    "LiFSI": {"ACN": 0.7523, "DME": 0.7951, "Me2O": 0.7716, "DEE": 0.7201, "FM": 0.8996},
    "LiTFSI": {"ACN": 0.8655, "DME": 0.9147, "Me2O": 0.8877, "DEE": 0.8284, "FM": 1.0350},
    "LiTPFA": {"ACN": 1.2811, "DME": 1.3540, "Me2O": 1.3139, "DEE": 1.2262, "FM": 1.5320},
}


def k_full():
    import pandas as pd
    d = pd.read_csv(data("experimental", "Fig2B.csv"))
    k = {(r.solvent, r.salt): r.k_r_minus_over_r_plus for r in d.itertuples()}
    for (s, a), v in k.items():
        assert abs(v - TABLE_S4[a][s]) < 6e-5, (s, a, v)
    return k


def t_ratio(s, a):
    d7, d19 = TABLE_S2[(s, a)][:2]
    return d7 / (d7 + d19)


def T_transf(s, a):
    r = TABLE_S2[(s, a)]
    return r[3] / r[2]


def letter(fig, x, y, s, size=10):
    fig.text(x, y, s, fontsize=size, fontweight="bold", ha="left", va="top")


def embed(ax, img, note=None):
    if img is not None:
        ax.imshow(np.asarray(img), interpolation="lanczos")
    ax.set_axis_off()


def flatten(path):
    if path is None:
        return None
    im = Image.open(path).convert("RGBA")
    bg = Image.new("RGBA", im.size, "white")
    bg.alpha_composite(im)
    return bg.convert("RGB")


def save(fig, name):
    png = os.path.join(OUT, name + ".png")
    fig.savefig(png, dpi=DPI)
    plt.close(fig)
    w, h = Image.open(png).size
    print(f"{name}: {w}x{h} px")
    return png
