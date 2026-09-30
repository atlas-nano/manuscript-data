"""Figure 4: multi-node GPU strong scaling and the time per step by task. Reads ../data/scaling.csv."""
import csv
import os

import matplotlib.pyplot as plt

import figstyle
from figstyle import BLUE, GREEN, INK, MUTED, PURPLE, VERM

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, os.pardir, "data", "scaling.csv")
figstyle.apply()

rows = list(csv.DictReader(open(CSV)))
sizes = sorted({int(r["natoms"]) for r in rows})

def series(n):
    s = sorted((int(r["gpus"]), float(r["steps_per_s"])) for r in rows
               if int(r["natoms"]) == n)
    return [g for g, _ in s], [v for _, v in s]

fig, (axL, axR) = plt.subplots(1, 2, figsize=(figstyle.TEXTWIDTH, 2.55),
                               gridspec_kw={"width_ratios": [1.02, 1.0], "wspace": 0.30})

style = {sizes[0]: (BLUE, "o", "-"), sizes[1]: (VERM, "s", "-")}
for n in sizes:
    g, v = series(n)
    c, mk, ls = style[n]
    axL.plot(g, v, ls, color=c, marker=mk, ms=4.5, lw=1.4, zorder=4,
             label=f"{n/1000:.0f}k atoms")
    axL.plot(g, [v[0] * gg / g[0] for gg in g], ls=(0, (3, 2.2)), color=c, lw=0.9,
             alpha=0.55, zorder=2)

    for gg, vv in zip(g, v):
        axL.annotate(f"{vv:.2f}", (gg, vv), textcoords="offset points",
                     xytext=(0, 7), fontsize=6.4, color=c, ha="center")

axL.set_xscale("log", base=2)
axL.set_yscale("log")
axL.set_yticks([1, 2, 3, 5, 10])
axL.get_yaxis().set_major_formatter(plt.matplotlib.ticker.ScalarFormatter())
axL.get_yaxis().set_minor_formatter(plt.matplotlib.ticker.NullFormatter())
axL.set_ylim(1.0, 16)
axL.set_xticks([4, 8, 16])
axL.set_xticklabels(["4\n(1 node)", "8\n(2)", "16\n(4)"])
axL.set_xlabel("NVIDIA A100 devices")
axL.set_ylabel("timesteps s$^{-1}$")
axL.set_title("(a)  strong scaling", loc="left", fontsize=8)
axL.legend(frameon=False, loc="upper left", handlelength=1.8)
axL.text(0.30, 0.035, "dashed: ideal scaling from the 4-device point", transform=axL.transAxes,
         ha="left", va="bottom", fontsize=6.4, color=MUTED, style="italic")
axL.grid(True, which="major", color=figstyle.GRID, lw=0.5, zorder=0)
axL.set_axisbelow(True)

CATS = [("modify", "charge solve (Modify)", PURPLE),
        ("kspace", "force-time mesh (Kspace)", GREEN),
        ("_rest", "pair, bond, comm, neighbor", MUTED)]
labels, bottoms = [], []
order = [r for r in rows if int(r["natoms"]) == sizes[1]]
order.sort(key=lambda r: int(r["gpus"]))
x = range(len(order))
base = [0.0] * len(order)
for key, lab, col in CATS:
    if key == "_rest":
        vals = [sum(float(r[k]) for k in ("pair", "bond", "comm", "neigh", "output", "other"))
                for r in order]
    else:
        vals = [float(r[key]) for r in order]
    axR.bar(x, vals, bottom=base, width=0.62, color=col, label=lab,
            edgecolor="white", linewidth=0.9, zorder=3)
    base = [b + v for b, v in zip(base, vals)]
axR.set_xticks(list(x))
axR.set_xticklabels([f"{r['gpus']}" for r in order])
axR.set_xlabel("A100 devices")
axR.set_ylabel("share of the timestep (%)")
axR.set_ylim(0, 100)
axR.set_title(f"(b)  where the step goes, {sizes[1]/1000:.0f}k atoms", loc="left", fontsize=8)
axR.legend(frameon=False, fontsize=6.3, loc="center left", bbox_to_anchor=(1.01, 0.5),
           handlelength=1.1, labelspacing=0.9)
for i, r in enumerate(order):
    axR.text(i, float(r["modify"]) / 2, f"{float(r['modify']):.0f}%", ha="center",
             va="center", fontsize=6.6, color="white", zorder=5)
axR.grid(True, axis="y", color=figstyle.GRID, lw=0.5, zorder=0)
axR.set_axisbelow(True)

fig.savefig(os.path.join(HERE, "fig_gpuscale.pdf"), dpi=300, bbox_inches="tight")
fig.savefig(os.path.join(HERE, "fig_gpuscale.png"), dpi=300, bbox_inches="tight")
print("wrote fig_gpuscale.pdf / .png")
