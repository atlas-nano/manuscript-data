"""Draws panel F of Figure 3 from data/simulation/pairing_cycle.json and the xPT tables."""
import json, os, statistics as st
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

import fig_sources as SRC
from paths import OUT as OUTDIR

OUT = str(OUTDIR / ("fig3f_plate" if os.environ.get("NOHEADER", "1") == "1" else "fig3f"))
SYSTEMS = [(s, v) for v in ("FM", "ACN") for s in ("LiTPFA", "LiPF6", "LiBF4")]
PAIRED = {("LiPF6", "FM"), ("LiBF4", "FM")}
GGLA_CUT = -10.0
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e6e5e1"
C_DISS, C_SOLV, C_TDS = "#d97b66", "#4d86c6", "#a9a8a3"
FS_LABEL, FS_TICK, FS_NOTE = 8.2, 7.8, 6.9
MINUS = "−"


def data():
    A, RAW = SRC._per_box(SRC.TABLE_A), SRC._per_box(SRC.TABLE_RAW)
    diss, tds = SRC.gated_diss_tds()
    pc = json.load(open(SRC.data("simulation", "pairing_cycle.json")))
    rows = []
    for salt, solv in SYSTEMS:
        key = f"{salt.lower()}_{solv.lower()}"
        a, r = A[key], RAW[key]
        U = [x["c_ggla"] for x in a]
        n_paired = sum(u < GGLA_CUT for u in U)
        d = dict(solvent=solv, salt=salt, n_cells=len(a), n_paired=n_paired,
                 U_LA_median=st.median(U),
                 dG_diss=diss[salt], minus_TdS_ideal=-tds[solv],
                 dG_solv_total0=st.median([x["dG_solv"] for x in a]),
                 dG_mix_total0=st.median([x["dG_mix"] for x in a]),
                 dG_solv_B=st.median([x["dG_solv"] + x["c_ggla"] for x in a]),
                 dG_mix_B=st.median([x["dG_mix"] + x["c_ggla"] for x in a]),
                 dG_mix_B_raw=st.median([x["dG_mix"] + x["c_ggla"] for x in r]),
                 dG_mix_B_bootstrap_sd=pc[f"{salt}/{solv}"]["totals"]["B"]["s_med"])
        d["closure_B"] = d["dG_solv_B"] + d["dG_diss"] + d["minus_TdS_ideal"] - d["dG_mix_B"]
        d["plotted"] = "dG_diss only (contact pair in 8/8 cells)" if (salt, solv) in PAIRED else "full stack + net"
        rows.append(d)
    return rows


def panel(ax, rows):
    BOX = dict(facecolor="white", edgecolor="none", pad=1.0)
    W = 0.62
    for i, d in enumerate(rows):
        paired = (d["salt"], d["solvent"]) in PAIRED
        ax.bar(i, d["dG_diss"], width=W, color=C_DISS, zorder=3, label="$\\Delta G_{\\mathrm{diss}}$" if i == 0 else None)
        if paired:
            continue
        ax.bar(i, d["dG_solv_B"], width=W, color=C_SOLV, zorder=3, label="$\\Delta G_{\\mathrm{solv}}$" if i == 0 else None)
        ax.bar(i, d["minus_TdS_ideal"], width=W, bottom=d["dG_solv_B"], color=C_TDS, zorder=3,
               edgecolor="white", linewidth=0.6, label="$-T\\Delta S_{\\mathrm{ideal}}$" if i == 0 else None)
        g = d["dG_mix_B"]
        ax.plot([i], [g], marker="D", ms=4.4, color=INK, ls="none", zorder=6,
                label="net $\\Delta G_{\\mathrm{mix}}$" if i == 0 else None)
        if (d["salt"], d["solvent"]) == ("LiTPFA", "FM"):
            ax.annotate("≈ 0", (i, g), textcoords="offset points", xytext=(0, -10), ha="center", va="top",
                        fontsize=FS_NOTE, color=INK, zorder=7, bbox=BOX)
        else:
            ax.annotate(f"{g:+.0f}".replace("-", MINUS), (i, g), textcoords="offset points",
                        xytext=(0, 10 if g >= 0 else -10), ha="center", va="bottom" if g >= 0 else "top",
                        fontsize=FS_NOTE, color=INK, zorder=7, bbox=BOX)
    ax.text(1.5, -14, "contact\npair in\n8/8 cells", ha="center", va="top", fontsize=FS_NOTE, color=INK2,
            linespacing=1.1, zorder=7)
    ax.axhline(0, color=INK, lw=1.0, zorder=4)
    ax.axvline(2.5, color=INK2, lw=0.7, ls=":", zorder=2)
    ax.set_xticks(np.arange(len(rows)))
    ax.set_xticklabels([{"LiPF6": "LiPF\u2086", "LiBF4": "LiBF\u2084"}.get(d["salt"], SRC.LAB[d["salt"]]) for d in rows], fontsize=FS_NOTE, rotation=35, ha="right",
                       rotation_mode="anchor")
    for xf, lab in ((1.0, "in FM"), (4.0, "in ACN")):
        ax.text((xf + 0.7) / 6.0, -0.295, lab, transform=ax.transAxes, ha="center", va="top",
                fontsize=FS_NOTE + 0.5, color=INK2, fontweight="bold")
    ax.set_xlim(-0.7, len(rows) - 0.3)
    ax.set_ylim(-268, 206)
    ax.set_yticks([-150, -100, -50, 0, 50, 100, 150])
    ax.set_yticklabels([f"{t}".replace("-", MINUS) for t in (-150, -100, -50, 0, 50, 100, 150)])
    ax.set_ylabel("Free energy (kcal mol$^{-1}$)", fontsize=FS_LABEL)
    ax.tick_params(labelsize=FS_TICK)
    ax.grid(True, axis="y", color=GRID, lw=0.5, zorder=0)
    ax.legend(fontsize=FS_NOTE, frameon=False, loc="lower left", ncol=2, handlelength=1.0, columnspacing=0.8,
              handletextpad=0.4, borderpad=0.0, labelspacing=0.22)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


def main():
    rows = data()
    plt.rcParams.update({"font.family": "DejaVu Sans", "mathtext.fontset": "dejavusans", "font.size": 8.4,
                         "axes.linewidth": 0.9, "axes.edgecolor": INK, "axes.labelcolor": INK,
                         "xtick.color": INK, "ytick.color": INK, "figure.facecolor": "white",
                         "savefig.facecolor": "white", "axes.unicode_minus": True})
    Wd, Ht = 2.40, 2.50
    fig = plt.figure(figsize=(Wd, Ht), dpi=300)
    ax = fig.add_axes([0.60 / Wd, 0.70 / Ht, 1.68 / Wd, 1.70 / Ht])
    panel(ax, rows)
    if os.environ.get("NOHEADER", "1") != "1":
        fig.text(0.04 / Wd, (Ht - 0.04) / Ht, "F", fontsize=9.5, fontweight="bold", color=INK, va="top")
        fig.text(0.24 / Wd, (Ht - 0.05) / Ht, "Mixing free energy", fontsize=8.1, color=INK, va="top")
    fig.savefig(OUT + ".png", dpi=300)
    print(f"{os.path.basename(OUT)}.png written")


if __name__ == "__main__":
    main()
