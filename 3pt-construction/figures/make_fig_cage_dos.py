#!/usr/bin/env python3
r"""Figure 4 — gas|cage|solid decomposition of the velocity density of states,
four panels:
  (a) TIP4P/2005 water, TRANSLATIONAL cage
  (b) TIP4P/2005 water, ROTATIONAL cage (libration band)
  (c) liquid Al  (900 K)  — translational cage
  (d) liquid Au (1300 K)  — translational cage (soft, heavy: low-frequency,
                            cage-dominated DoS — contrasts Al)

Cage columns (cage_G1_<channel>) come from py-xPT .pwr files written with
show_2pt_split + the cage corrections (trans: cage_entropy; rot: cage_entropy_rot).
x-axis in cm^-1 (consistent with the text); the THz equivalent of each cage peak
(1 THz = 33.356 cm^-1) is given in parentheses.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
PWR  = HERE.parent / "data" / "pwr"
WATER = PWR / "tip4p2005.298K.1728.3PT.pwr"
METAL = PWR
CM_PER_THZ = 33.35641  # 1 THz = 33.356 cm^-1

GAS_C, CAGE_C, SOL_C = "#4C9BE6", "#E6A23C", "#9B59B6"


def load_pwr(path):
    names = None
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                t = line.lstrip("#").split()
                if t and t[0].startswith("freq"):
                    names = t
            elif line.strip():
                break
    return names, np.loadtxt(path)


def panel(ax, path, channel, title, xmax, ann_xy=(0.40, 0.74)):
    names, d = load_pwr(path)
    ci = names.index
    f = d[:, 0]                        # cm^-1
    DoS = d[:, ci(f"DoS_G1_{channel}")]
    gas = np.clip(d[:, ci(f"gas_G1_{channel}")], 0, None)
    sol = np.clip(d[:, ci(f"sol_G1_{channel}")], 0, None)
    cage = np.clip(d[:, ci(f"cage_G1_{channel}")], 0, None)
    dnu = f[1] - f[0]
    Ig, Ic, Is = (np.trapezoid(x, dx=dnu) for x in (gas, cage, sol))
    m = f <= xmax
    ax.stackplot(f[m], gas[m], cage[m], sol[m],
                 colors=[GAS_C, CAGE_C, SOL_C], alpha=0.85,
                 labels=[f"gas  ∫={Ig:.2f}", f"cage  ∫={Ic:.2f}", f"solid  ∫={Is:.2f}"])
    ax.plot(f[m], DoS[m], "k-", lw=1.1, label="total DoS")
    ymax = float(np.max(DoS[m])) * 1.08
    # Report the cage-band CENTROID (first moment), not the argmax peak: the cage
    # band is broad and flat-topped, so its argmax is noise-sensitive (it flickers
    # ~14 cm^-1 between equivalent recomputations), whereas the centroid is robust
    # and reproducible. Annotate at the centroid frequency.
    fcen = float(np.trapezoid(f[m] * cage[m], dx=dnu) / np.trapezoid(cage[m], dx=dnu))
    ycen = float(np.interp(fcen, f[m], cage[m]))
    ax.annotate(f"cage centroid\n{fcen:.0f} cm$^{{-1}}$ ({fcen/CM_PER_THZ:.1f} THz)",
                xy=(fcen, ycen), xycoords="data",
                xytext=ann_xy, textcoords="axes fraction",
                fontsize=7.5, color="#d62728",
                ha="left", va="top",
                arrowprops=dict(arrowstyle="->", color="#d62728", lw=0.8))
    ax.set_title(title, fontsize=9.5)
    ax.set_xlim(0, xmax)
    ax.set_ylim(0, ymax)
    ax.legend(fontsize=7, frameon=False, loc="upper right")
    ax.grid(ls=":", lw=0.4, alpha=0.5)
    print(f"  {title[:34]:34s} cage_centroid={fcen:.0f} cm-1 ({fcen/CM_PER_THZ:.2f} THz)  "
          f"int gas/cage/sol = {Ig:.2f}/{Ic:.2f}/{Is:.2f}")
    return fcen, (Ig, Ic, Is)


def main():
    plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.8,
                         "figure.dpi": 150, "savefig.bbox": "tight"})
    fig, axes = plt.subplots(2, 2, figsize=(9.2, 6.6))
    (aa, ab), (ac, ad) = axes

    panel(aa, WATER, "trans", "(a) TIP4P/2005 water — translational cage", 400.0)
    panel(ab, WATER, "angul", "(b) TIP4P/2005 water — rotational cage", 900.0,
          ann_xy=(0.04, 0.92))   # peak is on the right; park label upper-left
    panel(ac, METAL/"Al_SC_900K_3pt_split.pwr", "trans", "(c) liquid Al (900 K)", 400.0)
    panel(ad, METAL/"Au_SC_1300K_3pt_split.pwr", "trans", "(d) liquid Au (1300 K)", 300.0)

    for ax in (ac, ad):
        ax.set_xlabel(r"frequency  (cm$^{-1}$)")
    for ax in (aa, ac):
        ax.set_ylabel("DoS  (DoF / particle / cm$^{-1}$)")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(HERE / f"fig_cage_dos_partition.{ext}")
    plt.close(fig)
    print("wrote fig_cage_dos_partition.png/.pdf")


if __name__ == "__main__":
    main()
