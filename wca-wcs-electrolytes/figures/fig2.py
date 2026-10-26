"""Draws Figure 2 from data/experimental/Fig1C.csv, Fig2B.csv, Fig2C.csv and Fig2D.csv."""
import pandas as pd
from matplotlib.lines import Line2D
from figstyle import *
from PIL import ImageDraw

K = k_full()
xs = pd.read_csv(data("experimental", "Fig1C.csv")).drop_duplicates("solvent").set_index("solvent")

rowsB = []
for sp in SOLVENTS + ["LiTPFA (all solvents)"]:
    pairs = [(sp, a) for a in SALTS] if sp in SOLVENTS else [(s, "LiTPFA") for s in ["FM", "ACN", "Me2O", "DME", "DEE"]]
    for s, a in pairs:
        if (s, a) not in TABLE_S2:
            rowsB.append(dict(subpanel=sp, solvent=s, salt=a, k_r_minus_over_r_plus=K[(s, a)], D_7Li_1e9_m2_s=np.nan,
                              D_19F_1e9_m2_s=np.nan, Li_transport_number_t=np.nan, t_2dp=np.nan,
                              ideal_t_k_over_k_plus_1=K[(s, a)] / (K[(s, a)] + 1), plotted=False))
            continue
        d7, d19 = TABLE_S2[(s, a)][:2]
        t = d7 / (d7 + d19)
        rowsB.append(dict(subpanel=sp, solvent=s, salt=a, k_r_minus_over_r_plus=K[(s, a)], D_7Li_1e9_m2_s=d7,
                          D_19F_1e9_m2_s=d19, Li_transport_number_t=t, t_2dp=round(t, 2),
                          ideal_t_k_over_k_plus_1=K[(s, a)] / (K[(s, a)] + 1), plotted=True))
B = pd.DataFrame(rowsB)

FLU = "dielectric_fluidity_eps_r_over_eta_mPa-1_s-1"
rowsC = []
for (s, a), r in TABLE_S2.items():
    rowsC.append({"solvent": s, "salt": a, "inverse_Haven_ratio_H_inv": r[6], "H_inv_Dsigma_over_sum": r[2] / (r[0] + r[1]),
                  "eps_r": xs.loc[s, "eps_r"], "viscosity_mPa_s": xs.loc[s, "viscosity_mPa_s"], FLU: xs.loc[s, FLU]})
C = pd.DataFrame(rowsC)

rowsD = []
for (s, a), r in TABLE_S2.items():
    T = r[3] / r[2]
    inc = 0 < T < 1
    note = "Included" if inc else "Excluded: transference outside (0,1)"
    rowsD.append(dict(solvent=s, salt=a, concentration_M=0.25, Li_transport_number_t_x=t_ratio(s, a),
                      Li_transference_number_T_y=T, plotted=note))
D = pd.DataFrame(rowsD)

bx = {(r.solvent, r.salt): r.Li_transport_number_t for r in B[B.plotted].itertuples()}


def mk(ax, x, y, s, a, filled=True, ms=None, z=3):
    if s == "FM":
        ax.plot(x, y, marker="*", ms=ms or 10, mfc=FM_FACE, mec=FM_EDGE, mew=0.7, ls="none", zorder=6)
        return
    m = SALT_MARK[a]
    size = {"*": 7.0, "D": 3.6}.get(m, 4.2) if ms is None else ms
    c = SOLV[s]
    if filled:
        ax.plot(x, y, marker=m, ms=size, mfc=c, mec=c, mew=0.4, ls="none", zorder=z)
    else:
        ax.plot(x, y, marker=m, ms=size, mfc="none", mec=c, mew=0.8, ls="none", zorder=z)


A_H = 6.0 * CM
H_FIG = 12.4 * CM
S = 11.4 * CM / H_FIG


def Y(y):
    return y * S


def T(y):
    return 1 - (1 - y) * S


fig = plt.figure(figsize=(W_FULL, H_FIG))
axA = fig.add_axes([0.03, 1 - (0.15 * CM + A_H) / H_FIG, 0.33, A_H / H_FIG])
imgA = panel_image("fig2A_schematic.png")
if imgA is not None:
    embed(axA, flatten(imgA))
else:
    axA.axis("off")
letter(fig, 0.005, T(0.985), "A")

letter(fig, 0.385, T(0.985), "B")
x0, wB, gap = 0.445, 0.106, 0.0
subs = SOLVENTS + ["LiTPFA (all solvents)"]
axs = []
for i, sp in enumerate(subs):
    ax = fig.add_axes([x0 + i * wB + 0.004, T(0.60), wB - 0.008, Y(0.30)], sharey=axs[0] if axs else None)
    axs.append(ax)
    d = B[(B.subpanel == sp) & B.plotted]
    if sp in SOLVENTS:
        kk = np.linspace(0.4, 1.7, 100)
        ax.set_xlim(0.4, 1.7); ax.set_xticks([0.5, 1.0, 1.5])
        title, tc = SOLV_TEX[sp], SOLV[sp]
    else:
        kk = np.linspace(1.15, 1.65, 50)
        ax.set_xlim(1.15, 1.65); ax.set_xticks([1.2, 1.4, 1.6])
        title, tc = "LiTPFA", "#000000"
    ax.plot(kk, kk / (kk + 1), ls=(0, (3, 2)), color="k", lw=0.8, zorder=1)
    ax.axhline(0.5, color=NEUTRAL, lw=0.7, ls=":", zorder=1)
    ax.axhline(0.41, color="#9a9a9a", lw=0.6, zorder=1)
    for r in d.itertuples():
        mk(ax, r.k_r_minus_over_r_plus, r.Li_transport_number_t, r.solvent, r.salt)
    ax.set_ylim(0.3, 0.7)
    ax.set_title(title, color=tc, fontsize=7.5, fontweight="bold", pad=3)
    if i:
        plt.setp(ax.get_yticklabels(), visible=False)
    if sp == "ACN":
        ax.text(1.68, 0.355, "LP57\n$T_{\\mathrm{Li^+}}$ = 0.41", fontsize=5.5, color="#6e6e6e", ha="right", va="center")
    if sp == "DME":
        ax.text(0.80, 0.62, r"$t^{\mathrm{ideal}}_{\mathrm{Li^+}}$ = $k$/($k$ + 1)", fontsize=5.5, rotation=58,
                ha="center", va="center")
    if sp.startswith("LiTPFA"):
        ax.text(1.532, 0.575, "FM", fontsize=6, fontweight="bold", ha="center")
        ax.text(1.40, 0.43, "$t_+$ = $t_-$ = 0.5;\n$D_-$ = $D_+$", fontsize=5.5, ha="center", color=NEUTRAL)
axs[0].set_ylabel(r"$\mathrm{Li^{+}\ transport\ number,}\ t_{\mathrm{Li^{+}}}$")
fig.text(x0 + 2.5 * wB, T(0.535), r"$k$ = $r_-$/$r_+$", ha="center", fontsize=7)
handles = [Line2D([], [], marker=SALT_MARK[a], ls="none", color="k", ms=4.5 if a != "LiTPFA" else 7) for a in SALTS]
fig.legend(handles, [SALT_TEX[a] for a in SALTS], loc="upper center", bbox_to_anchor=(x0 + 2.5 * wB, T(0.995)), ncol=5,
           handletextpad=0.2, columnspacing=1.0, fontsize=6.5)

letter(fig, 0.005, Y(0.47), "C")
colsC = [("eps_r", "Dielectric constant,\n$\\varepsilon_r$", (0, 42), [0, 10, 20, 30, 40]),
         ("viscosity_mPa_s", "Viscosity,\n$\\eta$ (mPa s)", (0, 1.75), [0, 0.5, 1.0, 1.5]),
         ("dielectric_fluidity_eps_r_over_eta_mPa-1_s-1", "Dielectric fluidity,\n$\\varepsilon_r\\eta^{-1}$ (mPa$^{-1}$ s$^{-1}$)",
          (-5, 125), [0, 40, 80, 120])]
axc = []
for j, (col, lab, xl, xt) in enumerate(colsC):
    ax = fig.add_axes([0.115 + j * 0.135, Y(0.155), 0.118, Y(0.30)], sharey=axc[0] if axc else None)
    axc.append(ax)
    for _, r in C.iterrows():
        mk(ax, r[col], r.inverse_Haven_ratio_H_inv, r.solvent, r.salt, filled=(r.salt == "LiTPFA"))
    ax.set_yscale("log"); ax.set_ylim(1e-4, 2); ax.set_xlim(*xl); ax.set_xticks(xt)
    ax.set_xlabel(lab, fontsize=6.5, linespacing=1.1)
    if j:
        plt.setp(ax.get_yticklabels(), visible=False)
axc[0].set_ylabel(r"Inverse Haven ratio, $H^{-1}$")
fig.text(0.012, Y(0.30), "Increasing\nion correlation", rotation=90, fontsize=6, style="italic", ha="center", va="center",
         color=NEUTRAL)
fig.add_artist(matplotlib.patches.FancyArrowPatch((0.032, Y(0.44)), (0.032, Y(0.16)), transform=fig.transFigure,
                                                  arrowstyle="-|>", mutation_scale=7, color=NEUTRAL, lw=0.8))
sh = [Line2D([], [], marker="o", ls="none", mfc="none", mec=SOLV[s], mew=0.9, ms=4.5) for s in SOLVENTS] + \
     [Line2D([], [], marker="*", ls="none", mfc=FM_FACE, mec=FM_EDGE, mew=0.6, ms=8)]
fig.legend(sh, [SOLV_TEX[s] for s in SOLVENTS] + ["FM (LiTPFA)"], loc="lower center", bbox_to_anchor=(0.31, Y(-0.005)),
           ncol=5, handletextpad=0.2, columnspacing=1.0, fontsize=6.5)

letter(fig, 0.535, Y(0.47), "D")
ax = fig.add_axes([0.62, Y(0.11), 0.36, Y(0.345)])
ax.plot([0, 1], [0, 1], ls=(0, (4, 2.5)), color="k", lw=0.8, zorder=1)
for r in D[D.plotted.str.startswith("Included")].itertuples():
    mk(ax, r.Li_transport_number_t_x, r.Li_transference_number_T_y, r.solvent, r.salt, filled=(r.salt == "LiTPFA"))
ax.set_xlim(0, 1); ax.set_ylim(0, 1)
ax.set_xticks(np.arange(0, 1.01, 0.2)); ax.set_yticks(np.arange(0, 1.01, 0.2))
ax.set_xlabel(r"$\mathrm{Li^{+}\ transport\ number,}\ t_{\mathrm{Li^{+}}}$")
ax.set_ylabel(r"$\mathrm{Li^{+}\ transference\ number,}\ T_{\mathrm{Li^{+}}}$")
ax.text(0.05, 0.93, r"$T_{\mathrm{Li^+}}$ > $t_{\mathrm{Li^+}}$", style="italic", va="top", fontsize=6.5)
ax.text(0.95, 0.07, r"$T_{\mathrm{Li^+}}$ < $t_{\mathrm{Li^+}}$", style="italic", ha="right", fontsize=6.5)
save(fig, "Fig2")
