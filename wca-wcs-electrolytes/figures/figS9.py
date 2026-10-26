"""Draws Figure S9 from data/simulation/opes2d_surfaces.npz."""
import pandas as pd
from figstyle import *

NPZ = data("simulation", "opes2d_surfaces.npz")
KT_KCAL = 0.0019872043 * 293.15
VMAX = 18.0
d = np.load(NPZ)
he01 = d["h_edges_01"]; hc01 = 0.5 * (he01[:-1] + he01[1:]); cn005 = d["cn_grid_005"]
Fb = d["F_bias_h_cn_fine"].copy(); Fb -= np.nanmin(Fb)
he02 = d["h_edges_02"]; hc02 = 0.5 * (he02[:-1] + he02[1:]); cnc01 = d["cn_centres_01"]
Fr = d["F_h_cn"].copy(); Fr -= np.nanmin(Fr)
Fr_m = np.where(d["nwalk_h_cn"] >= 1, Fr, np.nan)
REG = (("F_cn_R", "#000000", (0, (4, 2)), "R, 11–14 Å (near-bulk layer)"),
       ("F_cn_L1", "#0072B2", "-", "L1, 3–6 Å"),
       ("F_cn_C", "#E69F00", "-", "C, 1.4–3 Å (contact)"))
F1 = {}
for key, *_ in REG:
    f = d[key].astype(float); nw = d["nwalk_" + key[2:]]
    f = np.where(nw >= 2, f, np.nan); f -= np.nanmin(f); F1[key] = (f, nw)

SNAP = panel_image("figS9A_render.png")
snap = Image.open(SNAP).convert("RGB") if SNAP is not None else Image.new("RGB", (1298, 1869), "white")
HT = 6.6
fig = plt.figure(figsize=(W_FULL, HT * CM))
AXY, AXH = 0.19, 0.70
HA = 5.4; WA = HA * snap.width / snap.height
axs = fig.add_axes([0.012, (HT - HA) / 2 / HT - 0.01, WA / 18.4, HA / HT]); axs.imshow(np.asarray(snap), interpolation="lanczos"); axs.set_axis_off()
cmap = plt.get_cmap("viridis").copy(); cmap.set_over("#fde725"); cmap.set_bad("white")
lev = np.linspace(0, VMAX, 37)
axa = fig.add_axes([0.265, AXY, 0.170, AXH]); axb = fig.add_axes([0.445, AXY, 0.170, AXH], sharey=axa)
cax = fig.add_axes([0.625, AXY, 0.010, AXH]); axc = fig.add_axes([0.800, AXY, 0.185, AXH])
for ax, X, Y, F, L in ((axa, cn005, hc01, Fb, "A"), (axb, cnc01, hc02, Fr_m, "B")):
    cs = ax.contourf(X, Y, np.clip(F, 0, VMAX + 1e-6), levels=lev, cmap=cmap, extend="max")
    ax.set_xlim(0, 6); ax.set_ylim(1.0, 15.0); ax.set_xticks(range(0, 7))
    ax.set_xlabel("")
    for y in (3.0, 6.0, 11.0, 14.0):
        ax.axhline(y, color="#9a9a9a", lw=0.6, ls=(0, (3, 3)))
    for y, t in ((2.0, "C"), (4.5, "L1"), (12.5, "R")):
        ax.text(5.85, y, t, color="black", ha="right", va="center", fontsize=6.5, fontweight="bold", bbox=dict(fc="white", ec="none", pad=0.6))
    ax.set_facecolor("white")
letter(fig, 0.005, 0.985, "A")
letter(fig, 0.212, 0.985, "B"); letter(fig, 0.437, 0.985, "C"); letter(fig, 0.745, 0.985, "D")
fig.text(0.44, 0.035, r"Li$^{+}$–F(FM) coordination, $n$", ha="center", va="bottom", fontsize=7)
axa.set_ylabel(r"Li$^{+}$ height above the metal face, $h$ (Å)")
axa.set_title("Bias-derived F(h, n)", fontsize=7, loc="center")
axb.set_title("Reweighted F(h, n)", fontsize=7, loc="center")
plt.setp(axb.get_yticklabels(), visible=False)
axa.text(0.15, 14.75, "zero net\nelectrode charge", color="black", fontsize=5.5, va="top", bbox=dict(fc="white", ec="none", pad=0.6))
cb = fig.colorbar(cs, cax=cax, ticks=np.arange(0, VMAX + 1, 3)); cb.set_label(r"$F$ ($k_\mathrm{B}T$)", fontsize=7, labelpad=1)
cax.yaxis.set_ticks_position("right"); cax.yaxis.set_label_position("right")
cax2 = cax.twinx(); cax2.set_ylim(0, VMAX * KT_KCAL); cax2.set_ylabel(r"$F$ (kcal mol$^{-1}$)", fontsize=7)
cax2.spines["right"].set_position(("outward", 24)); cax2.set_yticks([0, 2, 4, 6, 8, 10])
cax.yaxis.tick_right(); cax.yaxis.set_label_position("right")
for key, col, ls, lab in REG:
    axc.plot(cn005, F1[key][0], color=col, ls=ls, lw=1.1, label=lab)
axc.set_xlim(0, 6); axc.set_ylim(0, 40); axc.set_yticks([0, 5, 10, 15, 20, 25]); axc.set_xticks(range(0, 7))
axc.set_xlabel(""); fig.text(0.8925, 0.035, r"Li$^{+}$–F(FM) coordination, $n$", ha="center", va="bottom", fontsize=7); axc.set_ylabel(r"$F(n\,|\,\mathrm{region})$ ($k_\mathrm{B}T$)")
axc.legend(frameon=False, fontsize=5.5, loc="upper center", handlelength=2.2)
png = os.path.join(OUT, "figS9.png")
fig.savefig(png, dpi=600); plt.close(fig)
print("figS9.png written")
