"""Draws Figure S5 from data/experimental/figS5.csv and Fig2B.csv."""
import pandas as pd
from matplotlib.lines import Line2D
from figstyle import *

K = k_full()
SALTC = {"LiBF4": "#333333", "LiPF6": "#D55E00", "LiFSI": "#0072B2", "LiTFSI": "#CC79A7", "LiTPFA": "#E6A800"}
SALTM = {"LiBF4": "s", "LiPF6": "o", "LiFSI": "^", "LiTFSI": "v", "LiTPFA": "D"}
rows = []
for s in SOLVENTS + ["FM"]:
    for a in SALTS:
        if s == "FM" and a != "LiTPFA":
            continue
        k = K[(s, a)]; ti = k / (k + 1)
        if (s, a) not in TABLE_S2:
            rows.append(dict(solvent=s, salt=a, k_table_S4=k, t_ideal=ti, t_table_S2=np.nan, deviation_percent_signed=np.nan,
                             plotted=False)); continue
        t = t_ratio(s, a)
        rows.append(dict(solvent=s, salt=a, k_table_S4=k, t_ideal=ti, t_table_S2=t,
                         deviation_percent_signed=100 * (t / ti - 1), plotted=(s != "FM")))
d = pd.DataFrame(rows)
lt = d[(d.salt == "LiTPFA")].deviation_percent_signed.abs()

fig, ax = plt.subplots(figsize=(W_COL, 5.6 * CM))
fig.subplots_adjust(left=0.2, right=0.97, bottom=0.17, top=0.80)
xpos = {s: i for i, s in enumerate(SOLVENTS)}
off = {a: (j - 2) * 0.09 for j, a in enumerate(SALTS)}
ax.axhline(0, color="#9a9a9a", lw=0.6, zorder=1)
for r in d[d.plotted].itertuples():
    ax.plot(xpos[r.solvent] + off[r.salt], r.deviation_percent_signed, marker=SALTM[r.salt], ms=4.2, ls="none",
            mfc=SALTC[r.salt], mec="k" if r.salt == "LiTPFA" else SALTC[r.salt], mew=0.5, zorder=3)
ax.set_xticks(range(4)); ax.set_xticklabels([SOLV_TEX[s] for s in SOLVENTS])
ax.set_xlim(-0.5, 3.5); ax.set_ylim(-15, 40)
ax.set_ylabel(r"Deviation from $t^{\mathrm{ideal}}_{\mathrm{Li^{+}}}$ (%)")
ax.set_xlabel("Solvent")
h = [Line2D([], [], marker=SALTM[a], ls="none", mfc=SALTC[a], mec="k" if a == "LiTPFA" else SALTC[a], mew=0.5, ms=4.2)
     for a in SALTS]
fig.legend(h, [SALT_TEX[a] for a in SALTS], loc="upper center", ncol=3, bbox_to_anchor=(0.56, 1.0), handletextpad=0.1,
           columnspacing=0.8, fontsize=6)
save(fig, "figS5")
