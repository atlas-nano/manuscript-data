"""Draws Figure 1 from data/experimental/Fig1B.csv, Fig1C.csv and Fig1D.csv."""
import pandas as pd
from matplotlib.lines import Line2D
from figstyle import *

B = pd.read_csv(data("experimental", "Fig1B.csv"))
C = pd.read_csv(data("experimental", "Fig1C.csv"))
D = pd.read_csv(data("experimental", "Fig1D.csv"))
FLU = "dielectric_fluidity_eps_r_over_eta_mPa-1_s-1"
DS = "charge_diffusivity_D_sigma_1e-9_m2_s-1"


def mk(ax, x, y, s, a, filled=True, ms=None, z=3):
    if s == "FM":
        ax.plot(x, y, marker="*", ms=ms or 10, mfc=FM_FACE, mec=FM_EDGE, mew=0.7, ls="none", zorder=6); return
    m = SALT_MARK[a]; size = {"*": 7.0, "D": 3.6}.get(m, 4.2) if ms is None else ms; c = SOLV[s]
    ax.plot(x, y, marker=m, ms=size, mfc=c if filled else "none", mec=c, mew=0.4 if filled else 0.8, ls="none", zorder=z)


A_W = 0.965 * W_FULL
A_H = A_W * 741 / 1739
H_TOP = 0.15 * CM + A_H + 0.30 * CM
H_FIG = 11.4 * CM + H_TOP
S = 11.4 * CM / H_FIG


def Y(y):
    return y * S


fig = plt.figure(figsize=(W_FULL, H_FIG))
axA = fig.add_axes([0.03, 1 - (0.15 * CM + A_H) / H_FIG, 0.965, A_H / H_FIG])
imgA = flatten(panel_image("fig1A_schematic.png"))
embed(axA, imgA if imgA is not None else Image.new("RGB", (1739, 741), "white"))
letter(fig, 0.005, 1 - 0.1 * CM / H_FIG, "A")

letter(fig, 0.005, Y(0.985), "B")
x0, wB = 0.075, 0.182
subs = SOLVENTS + ["LiTPFA (all solvents)"]
axs = []
for i, sp in enumerate(subs):
    ax = fig.add_axes([x0 + i * wB + 0.004, Y(0.60), wB - 0.008, Y(0.30)], sharey=axs[0] if axs else None)
    axs.append(ax)
    d = B[B.subpanel == sp]
    ax.axhline(8.8, color="#9a9a9a", lw=0.6, zorder=1)
    for _, r in d.iterrows():
        mk(ax, r.k_r_minus_over_r_plus, r.conductivity_mS_per_cm, r.solvent, r.salt,
           z=4 if (sp == "LiTPFA (all solvents)" and r.solvent == "ACN") else 3)
        if pd.notna(r.reported_uncertainty_mS_per_cm) and r.reported_uncertainty_mS_per_cm > 1e-5:
            ax.errorbar(r.k_r_minus_over_r_plus, r.conductivity_mS_per_cm, yerr=r.reported_uncertainty_mS_per_cm,
                        fmt="none", ecolor="k", elinewidth=0.5, capsize=1, zorder=2)
    if sp in SOLVENTS:
        ax.set_xlim(0.4, 1.7); ax.set_xticks([0.5, 1.0, 1.5]); title, tc = SOLV_TEX[sp], SOLV[sp]
        ax.text(0.05, 0.97, rf"$\varepsilon_r$ = {d.eps_r_label.iloc[0]}", transform=ax.transAxes, fontsize=5.5,
                va="top", fontweight="bold")
    else:
        ax.set_xlim(1.15, 1.65); ax.set_xticks([1.2, 1.4, 1.6]); title, tc = "LiTPFA", "#000000"
        ax.text(1.532, 25.8, "FM", fontsize=6, fontweight="bold", ha="center")
    ax.set_ylim(-1, 30); ax.set_yticks([0, 10, 20, 30])
    ax.set_title(title, color=tc, fontsize=7.5, fontweight="bold", pad=3)
    if i:
        plt.setp(ax.get_yticklabels(), visible=False)
    if sp == "ACN":
        ax.text(1.66, 7.6, "LP57\n(8.8 mS cm$^{-1}$)", fontsize=5.3, color="#6e6e6e", va="top", ha="right")
    if sp == "DEE":
        ins = ax.inset_axes([0.30, 0.42, 0.62, 0.40])
        for _, r in d[d.salt != "LiTPFA"].iterrows():
            mk(ins, r.k_r_minus_over_r_plus, r.conductivity_mS_per_cm, r.solvent, r.salt, ms=3.2)
        ins.set_xlim(0.52, 0.88); ins.set_ylim(-0.004, 0.032); ins.set_xticks([0.6, 0.7, 0.8])
        ins.set_yticks([0, 0.01, 0.02, 0.03]); ins.tick_params(labelsize=5, length=1.5, pad=1)
        ins.set_title("Zoomed in", fontsize=5.3, pad=1.5)
axs[0].set_ylabel(r"Conductivity (mS cm$^{-1}$)")
fig.text(x0 + 2.5 * wB, Y(0.535), r"$k$ = $r_-$/$r_+$", ha="center", fontsize=7)
handles = [Line2D([], [], marker=SALT_MARK[a], ls="none", color="k", ms=4.5 if a != "LiTPFA" else 7) for a in SALTS]
fig.legend(handles, [SALT_TEX[a] for a in SALTS], loc="upper center", bbox_to_anchor=(x0 + 2.5 * wB, Y(0.995)), ncol=5,
           handletextpad=0.2, columnspacing=1.0, fontsize=6.5)

letter(fig, 0.005, Y(0.47), "C")
colsC = [("eps_r", "Dielectric constant,\n$\\varepsilon_r$", (0, 42), [0, 10, 20, 30, 40]),
         ("viscosity_mPa_s", "Viscosity,\n$\\eta$ (mPa s)", (0, 1.75), [0, 0.5, 1.0, 1.5]),
         (FLU, "Dielectric fluidity,\n$\\varepsilon_r\\eta^{-1}$ (mPa$^{-1}$ s$^{-1}$)", (-5, 125), [0, 40, 80, 120])]
axc = []
for j, (col, lab, xl, xt) in enumerate(colsC):
    ax = fig.add_axes([0.115 + j * 0.135, Y(0.155), 0.118, Y(0.30)], sharey=axc[0] if axc else None)
    axc.append(ax)
    for _, r in C.iterrows():
        mk(ax, r[col], r[DS], r.solvent, r.salt, filled=(r.salt == "LiTPFA"))
    ax.set_yscale("log"); ax.set_ylim(1e-4, 10); ax.set_xlim(*xl); ax.set_xticks(xt)
    ax.set_xlabel(lab, fontsize=6.5, linespacing=1.1)
    if j:
        plt.setp(ax.get_yticklabels(), visible=False)
axc[0].set_ylabel(r"Charge diffusivity (10$^{-9}$ m$^2$ s$^{-1}$)")
fig.text(0.012, Y(0.30), "Decreasing\nion transport", rotation=90, fontsize=6, style="italic", ha="center", va="center",
         color=NEUTRAL)
fig.add_artist(matplotlib.patches.FancyArrowPatch((0.032, Y(0.44)), (0.032, Y(0.16)), transform=fig.transFigure,
                                                  arrowstyle="-|>", mutation_scale=7, color=NEUTRAL, lw=0.8))
sh = [Line2D([], [], marker="o", ls="none", mfc="none", mec=SOLV[s], mew=0.9, ms=4.5) for s in SOLVENTS] + \
     [Line2D([], [], marker="*", ls="none", mfc=FM_FACE, mec=FM_EDGE, mew=0.6, ms=8)]
fig.legend(sh, [SOLV_TEX[s] for s in SOLVENTS] + ["FM (LiTPFA)"], loc="lower center", bbox_to_anchor=(0.31, Y(-0.005)),
           ncol=5, handletextpad=0.2, columnspacing=1.0, fontsize=6.5)

letter(fig, 0.555, Y(0.47), "D")
ax = fig.add_axes([0.64, Y(0.11), 0.34, Y(0.345)])
ramp = {1.4: "#3f007d", 1.0: "#6a51a3", 0.5: "#8f8cc4", 0.25: "#b3b1d9"}
mks = {1.4: "^", 1.0: "s", 0.5: "o", 0.25: "*"}
for c in [1.4, 1.0, 0.5, 0.25]:
    d = D[np.isclose(D.concentration_M, c)].sort_values("temperature_C")
    ax.plot(d.temperature_C, d.conductivity_mS_per_cm, color=ramp[c], marker=mks[c], ms=3.6 if c != 0.25 else 5,
            lw=1.0, label={1.4: "1.4 M", 1.0: "1.0 M", 0.5: "0.5 M", 0.25: "0.25 M"}[c])
ax.set_xlim(-70, 65); ax.set_ylim(0, 45); ax.set_xticks(range(-60, 61, 20))
ax.set_xlabel("Temperature (°C)"); ax.set_ylabel(r"Conductivity (mS cm$^{-1}$)")
ax.legend(title="x M LiTPFA : FM", loc="lower center", bbox_to_anchor=(0.6, 0.0), ncol=2, fontsize=6, title_fontsize=6.5, handlelength=2.2)
save(fig, "Fig1")
