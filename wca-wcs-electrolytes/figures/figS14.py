"""Draws Figure S14 from data/experimental/figS14AB.csv and figS14C.csv."""
import pandas as pd
from matplotlib import cm
from figstyle import *

ab = pd.read_csv(data("experimental", "figS14AB.csv"))
C = pd.read_csv(data("experimental", "figS14C.csv"), dtype={"printed_label_text": str})
temps = [-78, -40, -20, 0, 20, 40]
MARK = dict(zip(temps, ["o", "s", "^", "D", "v", "P"]))
RAMP = {"WCA-WCS": cm.Reds, "IL-LGE": cm.Blues}
col = lambda e, t: RAMP[e](np.linspace(1.0, 0.42, 6)[temps.index(t)])

fig = plt.figure(figsize=(W_FULL, 14.0 * CM))
for lab, e, rect, ly in (("A", "WCA-WCS", [0.085, 0.59, 0.36, 0.325], 0.99),
                         ("B", "IL-LGE", [0.085, 0.06, 0.36, 0.325], 0.50)):
    ax = fig.add_axes(rect)
    for t in temps:
        s = ab[(ab.electrolyte == e) & (ab.temperature_C == t)].sort_values("frequency_Hz", ascending=False)
        ax.plot(s.Z_real_ohm, s.minus_Z_imag_ohm, ls="-", lw=0.5, color=col(e, t), marker=MARK[t], ms=2.6,
                mfc="white", mew=0.7, label=f"{t} °C", clip_on=True)
    ax.set_xlim(0, 3000); ax.set_ylim(0, 3000); ax.set_aspect("equal")
    ax.set_xlabel(r"$Z^{\prime}$ ($\Omega$)"); ax.set_ylabel(r"$-Z^{\prime\prime}$ ($\Omega$)")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=3, fontsize=6, handlelength=1.6, columnspacing=1.0,
              title=e, title_fontproperties={"weight": "bold", "size": 7}, borderaxespad=0.2)
    letter(fig, 0.01, ly, lab)

letter(fig, 0.50, 0.99, "C")
W = 0.38
specs = [("R_electrolyte", r"$R_{\mathrm{electrolyte}}$ ($\Omega$)", (1, 1e3)),
         ("R_CT", r"$R_{\mathrm{CT}}$ ($\Omega$)", (1, 1e4)),
         ("W_SEI", r"$W_{\mathrm{SEI}}$ ($\Omega$)", (10, 1e7))]
def lab_fmt(v):
    s = f"{v:.3g}"
    if v < 1000: return s
    e = int(np.floor(np.log10(v))); m = float(f"{v / 10 ** e:.3g}")
    return f"{m:g}" + r"$\times$10$^{%d}$" % e
axs = []
for i, (q, yl, lim) in enumerate(specs):
    ax = fig.add_axes([0.585, 0.70 - i * 0.305, 0.40, 0.27]); axs.append(ax)
    for j, e in enumerate(["WCA-WCS", "IL-LGE"]):
        for k, t in enumerate(temps):
            r = C[(C.fitted_quantity == q) & (C.electrolyte == e) & (C.temperature_C == t)].iloc[0]
            x = k + (j - 0.5) * W
            if r.plotted:
                ax.bar(x, r.value_ohm, width=W, color=ELEC[e], edgecolor="white", linewidth=0.6, bottom=lim[0] * 1e-3,
                       label=e if k == 0 else None, zorder=2)
                ax.text(x, r.value_ohm * 1.25, lab_fmt(r.value_ohm), ha="center", va="bottom", fontsize=5, zorder=3,
                        rotation=90 if r.value_ohm >= 1000 else 0)
            else:
                ax.text(x, lim[0] * 1.5, "n.r.", ha="center", va="bottom", fontsize=5, color=NEUTRAL, rotation=90)
    ax.set_yscale("log"); ax.set_ylim(lim[0], lim[1] * 3)
    ax.set_xticks(range(6)); ax.set_xlim(-0.6, 5.6)
    ax.set_xticklabels([str(t) for t in temps] if i == 2 else [])
    ax.set_ylabel(yl)
axs[0].legend(loc="upper right", fontsize=6, handlelength=1.2)
axs[2].set_xlabel("Temperature (°C)")
p = save(fig, "figS14")
out = ab.copy(); out["shown_in_axes"] = (out.Z_real_ohm <= 3000) & (out.minus_Z_imag_ohm <= 3000) & (out.minus_Z_imag_ohm >= 0)
