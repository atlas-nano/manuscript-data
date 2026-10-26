"""Draws Figure 4 from data/experimental/Fig4A.csv to Fig4F.csv."""
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from figstyle import *

A = pd.read_csv(data("experimental", "Fig4A.csv"))
Bd = pd.read_csv(data("experimental", "Fig4B.csv"))
Dd = pd.read_csv(data("experimental", "Fig4D.csv"))
E = pd.read_csv(data("experimental", "Fig4E.csv"))
F = pd.read_csv(data("experimental", "Fig4F.csv"))
fig = plt.figure(figsize=(W_FULL, 14.0 * CM))
ORDER = ["LHCE", "IL-LGE", "WCA-WCS"]

letter(fig, 0.005, 0.995, "A")
ax = fig.add_axes([0.075, 0.715, 0.47, 0.265])
for s in ORDER:
    d = A[(A.series == s) & (A.view == "main")]
    for _, g in d.groupby("svg_path_index", sort=False):
        ax.plot(g.time_h, g.voltage_V, color=ELEC[s], lw=0.45, solid_joinstyle="miter")
ax.set_xlim(0, 75); ax.set_ylim(-0.8, 0.6)
ax.set_xlabel("Time (h)", labelpad=1); ax.set_ylabel("Voltage (V)")
ax.legend([Line2D([], [], color=ELEC[s], lw=1.2) for s in ["WCA-WCS", "IL-LGE", "LHCE"]], ["WCA-WCS", "IL-LGE", "LHCE"],
          loc="upper left", fontsize=6, handlelength=1.6, borderaxespad=0.3)
ins = ax.inset_axes([0.07, 0.07, 0.40, 0.33])
for s in ORDER:
    d = A[(A.series == s) & (A.view == "inset_final_hour")]
    for _, g in d.groupby("svg_path_index", sort=False):
        ins.plot(g.time_h, g.voltage_V, color=ELEC[s], lw=0.7)
ins.set_xlim(74, 75); ins.set_ylim(-0.3, 0.3); ins.set_xticks([74.0, 74.2, 74.4, 74.6, 74.8, 75.0])
ins.set_yticks([-0.2, 0, 0.2]); ins.tick_params(labelsize=5, pad=1, length=1.5)
ins.set_facecolor("white")

letter(fig, 0.60, 0.995, "B")
ax = fig.add_axes([0.665, 0.715, 0.32, 0.265])
cds = [1, 2, 4, 8]
RAMP = {1: "#cbc9e2", 2: "#9e9ac8", 4: "#6a51a3", 8: "#3f007d"}
groups = ["LHCE", "IL-LGE", "WCA-WCS"]
w = 0.2
for gi, e in enumerate(groups):
    for ci, c in enumerate(cds):
        r = Bd[(Bd.electrolyte == e) & (Bd.current_density_mA_cm2 == c)].iloc[0]
        x = gi + (ci - 1.5) * w
        ax.bar(x, r.bar_height_mV_from_vector, width=w, color=RAMP[c], edgecolor="white", linewidth=0.5,
               label=f"{c} mA cm$^{{-2}}$" if gi == 0 else None)
        ax.text(x, r.bar_height_mV_from_vector + 6, f"{r.printed_label_mV:g}", ha="center", va="bottom", fontsize=5)
ax.set_xticks(range(3)); ax.set_xticklabels(groups); ax.set_xlim(-0.55, 2.55); ax.set_ylim(0, 450)
ax.set_ylabel("Average overpotential (mV)"); ax.set_xlabel("Electrolyte", labelpad=1)
ax.legend(loc="upper right", fontsize=5.5, handlelength=1.0, borderaxespad=0.2)

letter(fig, 0.005, 0.635, "C")
src = flatten(panel_image("fig4C_sem_source.png"))
BOX = {("WCA-WCS", 0): (346, 1174, 960, 1583), ("WCA-WCS", 1): (348, 1610, 960, 2018),
       ("IL-LGE", 0): (994, 1175, 1613, 1585), ("IL-LGE", 1): (988, 1613, 1607, 2015),
       ("LHCE", 0): (1643, 1174, 2255, 1582), ("LHCE", 1): (1643, 1605, 2255, 2013)}
INSET = 7
cw, ch, cx0, cy0 = 0.150, 0.135, 0.07, 0.48
for j, e in enumerate(["WCA-WCS", "IL-LGE", "LHCE"]):
    for i in range(2):
        x0, y0, x1, y1 = BOX[(e, i)]
        box = (x0 + INSET, y0 + INSET, x1 - INSET, y1 - INSET)
        im = src.crop(box) if src is not None else Image.new("RGB", (box[2] - box[0], box[3] - box[1]), "white")
        a = fig.add_axes([cx0 + j * (cw + 0.008), cy0 - i * (ch + 0.012), cw, ch])
        a.imshow(np.asarray(im), interpolation="lanczos", aspect="equal"); a.set_anchor("C")
        a.set_xticks([]); a.set_yticks([])
        for sp in a.spines.values():
            sp.set_edgecolor(ELEC[e]); sp.set_linewidth(1.6)
fig.text(0.045, cy0 + ch / 2, "Surface", rotation=90, ha="center", va="center", fontsize=7, fontweight="bold")
fig.text(0.045, cy0 - ch / 2 - 0.012, "Cross-section", rotation=90, ha="center", va="center", fontsize=7, fontweight="bold")

letter(fig, 0.60, 0.635, "D")
ax = fig.add_axes([0.665, 0.355, 0.27, 0.26]); ax2 = ax.twinx()
for e in ["WCA-WCS", "IL-LGE", "LHCE"]:
    m = ELEC_MARK[e]
    c = Dd[(Dd.electrolyte == e) & (Dd.quantity == "discharge_capacity")].dropna(subset=["value"])
    q = Dd[(Dd.electrolyte == e) & (Dd.quantity == "coulombic_efficiency")].dropna(subset=["value"])
    ax.plot(c.cycle, c.value, marker=m, ls="none", ms=2.4, mfc=ELEC[e], mec=ELEC[e], mew=0.3)
    ax2.plot(q.cycle, q.value, marker=m, ls="none", ms=2.4, mfc="white", mec=ELEC[e], mew=0.5)
ax.set_xlim(-2, 103); ax.set_ylim(-5, 250); ax2.set_ylim(88, 102)
ax.set_xlabel("Cycle", labelpad=1); ax.set_ylabel(r"Specific discharge capacity (mAh g$^{-1}$)", fontsize=6.5)
ax2.set_ylabel("Coulombic efficiency (%)", fontsize=6.5)
ax.legend([Line2D([], [], marker=ELEC_MARK[e], ls="none", mfc=ELEC[e], mec=ELEC[e], ms=3.5) for e in ELEC], list(ELEC),
          loc="center right", bbox_to_anchor=(1.0, 0.71), fontsize=6, handletextpad=0.1, borderaxespad=0.2)
ax.text(0.98, 0.24, "Li||NMC811, 2.8–4.4 V\nWCA-WCS 2.73, IL-LGE 2.43,\nLHCE 3.09 mAh cm$^{-2}$\nC/10,D/10 x3; C/3,D/3 x3;\n1C,D/3 x3; 2C,D/3 x3; 5C,D/3 (25 °C)",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=5, fontweight="bold")

def vc(ax, df, order, cols, lab_fmt):
    for s in order:
        for _, g in df[df.series == s].groupby(["branch", "svg_path_index"], sort=False):
            ax.plot(g.specific_capacity_mAh_g, g.voltage_V, color=cols[s], lw=0.9)
    ax.set_ylim(2.8, 4.4); ax.set_yticks(np.arange(2.8, 4.41, 0.2))
    ax.set_xlabel(r"Specific capacity (mAh g$^{-1}$)", labelpad=1); ax.set_ylabel("Voltage (V)")
    ax.legend([Line2D([], [], color=cols[s], lw=1.2) for s in order], [lab_fmt(s) for s in order], loc="upper right",
              fontsize=6, handlelength=1.6, borderaxespad=0.3)

letter(fig, 0.005, 0.30, "E")
ax = fig.add_axes([0.075, 0.06, 0.42, 0.225])
colE = {"20 C": "#d92120", "0 C": "#9ecae1", "-20 C": "#4292c6", "-40 C": "#2160a8", "-60 C": "#08306b"}
vc(ax, E, ["20 C", "0 C", "-20 C", "-40 C", "-60 C"], colE, lambda s: s.replace(" C", " °C").replace("-", "−"))
ax.set_xlim(-5, 325); ax.set_xticks(range(0, 301, 50))
ax.text(0.98, 0.05, "Li||NMC811, 2.8–4.4 V\n≈3 mAh cm$^{-2}$ (nominal)\nC/3,D/3", transform=ax.transAxes, fontsize=5.5,
        fontweight="bold", va="bottom", ha="right")
letter(fig, 0.53, 0.30, "F")
ax = fig.add_axes([0.60, 0.06, 0.385, 0.225])
colF = {"C/3-C/3": "#9ecae1", "C/2-C/3": "#4292c6", "1C-C/3": "#2160a8", "2C-C/3": "#08306b"}
vc(ax, F, ["C/3-C/3", "C/2-C/3", "1C-C/3", "2C-C/3"], colF, lambda s: s.replace("-", "–"))
ax.set_xlim(-3, 180)
ax.text(0.06, 0.05, "Li||NMC811, 2.8–4.4 V\n≈3 mAh cm$^{-2}$ (nominal)", transform=ax.transAxes, fontsize=5.5,
        fontweight="bold", va="bottom")
ax.text(0.97, 0.05, "−40 °C", transform=ax.transAxes, fontsize=7, fontweight="bold", ha="right", va="bottom")
save(fig, "Fig4")
