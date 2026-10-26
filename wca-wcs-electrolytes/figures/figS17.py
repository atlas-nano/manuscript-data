"""Draws Figure S17 from data/experimental/figS17/figS17_rate_cycling.csv."""
import pandas as pd
from matplotlib.lines import Line2D
from figstyle import *

CE_MIN_FRAC = 0.05
d = pd.read_csv(data("experimental", "figS17", "figS17_rate_cycling.csv"))
ORDER = ["WCA-WCS", "IL-LGE", "LHCE"]
assert sorted(d.electrolyte.unique()) == sorted(ORDER) and len(d) == 81

load, shown = {}, {}
for e in ORDER:
    g = d[d.electrolyte == e]
    f = g[g.charge_rate_C == 0.1]
    assert list(f.cycle) == [1, 2, 3]
    load[e] = f.discharge_current_density_mA_cm2.abs().mean() / 0.1
    keep = g.specific_discharge_capacity_mAh_g >= CE_MIN_FRAC * f.specific_discharge_capacity_mAh_g.mean()
    assert (keep.astype(int) == g.ce_shown).all(), e
    shown[e] = g[keep]

fig = plt.figure(figsize=(17.03 * CM, 7.57 * CM))
ax = fig.add_axes([0.075, 0.13, 0.85, 0.74]); ax2 = ax.twinx()
for e in ORDER:
    g, q, m = d[d.electrolyte == e], shown[e], ELEC_MARK[e]
    ax.plot(g.cycle, g.specific_discharge_capacity_mAh_g, marker=m, ls="none", ms=3.2, mfc=ELEC[e], mec=ELEC[e], mew=0.3)
    ax2.plot(q.cycle, q.coulombic_efficiency_percent, marker=m, ls="none", ms=3.2, mfc="white", mec=ELEC[e], mew=0.6)
ax.set_xlim(0.3, 27.7); ax.set_xticks([1, 5, 10, 15, 20, 25]); ax.set_ylim(-5, 300); ax2.set_ylim(40, 108)
ax.set_xlabel("Cycle", labelpad=1); ax.set_ylabel(r"Specific discharge capacity (mAh g$^{-1}$)")
ax2.set_ylabel("Coulombic efficiency (%)")

r = d[d.electrolyte == "WCA-WCS"].sort_values("cycle")
for e in ORDER[1:]:
    assert (d[d.electrolyte == e].sort_values("cycle").charge_rate_C.values == r.charge_rate_C.values).all()
seg, start = [], 0
rates, cyc = list(r.charge_rate_C), list(r.cycle)
for i in range(1, len(rates) + 1):
    if i == len(rates) or rates[i] != rates[start]:
        seg.append((cyc[start], cyc[i - 1], rates[start])); start = i
for a, b, c in seg:
    if a > 1:
        ax.axvline(a - 0.5, color=NEUTRAL, lw=0.5, ls=(0, (3, 2)), zorder=0)
    ax.text((a + b) / 2, 1.015, f"{c:g}C", transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=6.5,
            fontweight="bold")
fig.text(0.075, 0.965, "Charge rate", fontsize=6.5, ha="left", va="top")

h = [Line2D([], [], marker=ELEC_MARK[e], ls="none", mfc=ELEC[e], mec=ELEC[e], ms=3.8) for e in ORDER]
h += [Line2D([], [], marker="o", ls="none", mfc=NEUTRAL, mec=NEUTRAL, ms=3.8),
      Line2D([], [], marker="o", ls="none", mfc="white", mec=NEUTRAL, mew=0.6, ms=3.8)]
ax.legend(h, ORDER + ["Discharge capacity", "C.E."], loc="lower left", bbox_to_anchor=(0.0, 0.03), fontsize=6,
          handletextpad=0.1, borderaxespad=0.3, ncol=1, frameon=True, facecolor="white", edgecolor="none", framealpha=1)
ax.text(0.47, 0.10, "Li||NMC811, 2.8–4.4 V, 25 °C\n" + ", ".join(f"{e} {load[e]:.2f}" for e in ORDER)
        + " mAh cm$^{-2}$\ncharge as labelled (1C = 200 mA g$^{-1}$); discharge C/3\n(C/10 in cycles 1–3)",
        transform=ax.transAxes, ha="left", va="bottom", fontsize=5.5, fontweight="bold", zorder=5,
        bbox=dict(boxstyle="square,pad=0.3", fc="white", ec="none"))
fig.savefig(os.path.join(OUT, "figS17.pdf"))
save(fig, "figS17")
