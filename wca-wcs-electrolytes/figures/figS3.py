"""Draws Figure S3 from data/simulation/pairing_cycle.json and the xPT tables."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

import plotstyle as S
import fig_sources as SRC
import fig3f as F3F
from paths import OUT as OUTDIR

OUT = str(OUTDIR / "figS3")
MK = {"dissolves": "o", "at_zero": "o", "forbidden": "X"}


def classify(g, sd):
    if abs(g) <= 2 * sd:
        return "at_zero"
    return "dissolves" if g < 0 else "forbidden"


def main():
    rows = F3F.data()
    be = {k: v[0] for k, v in SRC.BE_TABLE_S8.items()}
    S.use_style()
    fig, ax = plt.subplots(figsize=(5.20, 3.95))
    fig.subplots_adjust(left=0.150, right=0.975, top=0.955, bottom=0.150)
    ax.set_xlim(137.0, 182.0)
    ax.set_ylim(-53.5, -13.5)
    ax.add_patch(Rectangle((138.5, -32.0), 17.5, 13.0, facecolor=SRC.SALT_COL["LiTPFA"], alpha=0.09, lw=0, zorder=0))
    ax.text(139.4, -31.4, "WCA–WCS design target", fontsize=S.FS_LEGEND, color=SRC.SALT_COL["LiTPFA"],
            fontweight="bold", ha="left", va="bottom", zorder=6)
    off = {"FM": (0, 11), "ACN": (0, -11)}
    va = {"FM": "bottom", "ACN": "top"}
    for d in rows:
        salt, solv = d["salt"], d["solvent"]
        g, sd = d["dG_mix_B"], d["dG_mix_B_bootstrap_sd"]
        state = classify(g, sd)
        x, y = d["dG_diss"], be[solv]
        hollow = state == "at_zero"
        ax.scatter([x], [y], s=96, marker=MK[state], facecolor="white" if hollow else SRC.SALT_COL[salt],
                   edgecolor=SRC.SALT_COL[salt] if hollow else S.MARKER_EDGE,
                   linewidth=1.6 if hollow else S.MARKER_EDGE_W, zorder=5)
        majority_free = (d["n_cells"] - d["n_paired"]) * 2 >= d["n_cells"]
        lbl = f"{SRC.LAB[salt]}/{solv}"
        num = None
        if majority_free:
            num = "≈ 0" if hollow else f"{g:+.0f}".replace("-", "−")
            lbl += "\n" + ("$\\approx$ 0" if hollow else num)
        ax.annotate(lbl, (x, y), textcoords="offset points", xytext=off[solv], ha="center", va=va[solv],
                    fontsize=S.FS_LEGEND, color=S.INK_2, zorder=6)
    ax.set_xlabel("$\\Delta G_{diss}$ (kcal mol$^{-1}$)   (← easier)", fontsize=S.FS_LABEL)
    ax.set_ylabel("Li$^+$–solvent binding energy (kcal mol$^{-1}$)\n(weaker →)", fontsize=S.FS_LABEL)
    S.style_axes(ax)
    ax.grid(True, axis="x", color=S.GRID, lw=0.5, zorder=0)
    ax.text(139.0, -35.4, "marker: sign of $\\Delta G_{mix}$", fontsize=S.FS_LEGEND, color=S.INK_2, ha="left",
            va="bottom", zorder=6)
    for i, (state, txt) in enumerate((("dissolves", "≤ 0"), ("forbidden", "> 0"))):
        xd = 139.6 + 13.6 * i
        ax.scatter([xd], [-37.6], s=44, marker=MK[state], facecolor="none", edgecolor=S.MARKER_EDGE,
                   linewidth=S.MARKER_EDGE_W, zorder=7)
        ax.text(xd + 1.1, -37.6, txt, fontsize=S.FS_LEGEND, color=S.INK_2, ha="left", va="center", zorder=7)
    kw = dict(bbox_inches=fig.bbox_inches, pad_inches=0.0)
    fig.savefig(OUT + ".png", dpi=300, **kw)
    print("figS3.png written")


if __name__ == "__main__":
    main()
