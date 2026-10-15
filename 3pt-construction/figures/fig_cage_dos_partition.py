"""Figure 4: gas, cage and solid decomposition of the density of states for TIP4P/2005 water and liquid Al and Au. Reads ../data/cage_dos_partition.npz."""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
NPZ = os.path.join(HERE, os.pardir, "data", "cage_dos_partition.npz")
OUT = os.path.join(HERE, "out")
os.makedirs(OUT, exist_ok=True)

CM_PER_THZ = 33.35641
GAS_C, CAGE_C, SOL_C = "#4C9BE6", "#E6A23C", "#9B59B6"
D = np.load(NPZ)


def panel(ax, key, title, xmax, ann_xy=(0.40, 0.74)):
    f, DoS, gas, cage, sol = (D[f"{key}_{n}"] for n in ("nu", "total", "gas", "cage", "solid"))
    Ig, Ic, Is = D[f"{key}_integrals"]
    dnu = f[1] - f[0]
    ax.stackplot(f, gas, cage, sol, colors=[GAS_C, CAGE_C, SOL_C], alpha=0.85,
                 labels=[f"gas  ∫={Ig:.2f}", f"cage  ∫={Ic:.2f}", f"solid  ∫={Is:.2f}"])
    ax.plot(f, DoS, "k-", lw=1.1, label="total DoS")
    ymax = float(np.max(DoS)) * 1.08
    fcen = float(np.trapezoid(f * cage, dx=dnu) / np.trapezoid(cage, dx=dnu))
    ycen = float(np.interp(fcen, f, cage))
    ax.annotate(f"cage centroid\n{fcen:.0f} cm$^{{-1}}$ ({fcen/CM_PER_THZ:.1f} THz)",
                xy=(fcen, ycen), xycoords="data", xytext=ann_xy, textcoords="axes fraction",
                fontsize=7.5, color="#d62728", ha="left", va="top",
                arrowprops=dict(arrowstyle="->", color="#d62728", lw=0.8))
    ax.set_title(title, fontsize=9.5)
    ax.set_xlim(0, xmax)
    ax.set_ylim(0, ymax)
    ax.legend(fontsize=7, frameon=False, loc="upper right")
    ax.grid(ls=":", lw=0.4, alpha=0.5)


def main():
    plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.8, "figure.dpi": 150, "savefig.bbox": "tight"})
    fig, axes = plt.subplots(2, 2, figsize=(9.2, 6.6))
    (aa, ab), (ac, ad) = axes
    panel(aa, "water_trans", "(a) TIP4P/2005 water — translational cage", 400.0)
    panel(ab, "water_rot", "(b) TIP4P/2005 water — rotational cage", 900.0, ann_xy=(0.04, 0.92))
    panel(ac, "Al", "(c) liquid Al (900 K)", 400.0)
    panel(ad, "Au", "(d) liquid Au (1300 K)", 300.0)
    for ax in (ac, ad):
        ax.set_xlabel(r"frequency  (cm$^{-1}$)")
    for ax in (aa, ac):
        ax.set_ylabel("DoS  (DoF / particle / cm$^{-1}$)")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_cage_dos_partition.pdf"))
    plt.close(fig)


if __name__ == "__main__":
    main()
