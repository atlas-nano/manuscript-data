"""Central result: molecules and a conductor in one solve (conductor_demo/, revised for the referee, 2026-09-27).

(a) smallest eigenvalue of the constrained operator (dense) vs the gold hardness, four arms; (b) water->metal charge
and the largest single-molecule charge vs an offset of the water electronegativities (eta = 5.172 eV); (c) charge
above the slab mid-plane (Gaussian densities) opposite a unit point charge at height z; (d) the self-energy width:
face charge relative to the Gauss limit in a field, and the sub-surface layer relative to the face.
Data: ../conductor_demo/{dense_exp1,offset,exp2b,wscan}.csv. Previous version: fig_conductor_v1.py.bak_20260927.
"""
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import figstyle as fs

HERE = Path(__file__).resolve().parent
DEMO = HERE.parent / "conductor_demo"
fs.apply()

ARMS = [  # key, label, colour, marker, dash
    ("qeq",       "global, no self-energy",       fs.VERM,   "s", (0, (4, 2))),
    ("qeq_gself", "global, self-energy",          fs.PURPLE, "D", (0, (1, 1.5))),
    ("mol",       "per molecule, no self-energy", fs.GREEN,  "^", (0, (6, 2, 1, 2))),
    ("samqeq",    "per molecule, self-energy",    fs.BLUE,   "o", "solid"),
]


def load(name):
    with open(DEMO / name) as f:
        return list(csv.DictReader(f))


def eta_ticks(ax):
    t = [5, 2, 1, 0.3, 0.1, 0.05]
    ax.set_xscale("log"); ax.invert_xaxis()
    ax.set_xticks(t); ax.set_xticklabels([f"{v:g}" for v in t]); ax.minorticks_off()


dense, off, e2, ws = load("dense_exp1.csv"), load("offset.csv"), load("exp2b.csv"), load("wscan.csv")
fig, axs = plt.subplots(2, 2, figsize=(fs.TEXTWIDTH, 4.3), constrained_layout=True)

# (a) dense lambda_min vs eta
ax = axs[0, 0]
ax.axhspan(-4.5, 0, color=fs.GRID, alpha=0.45, lw=0)
ax.text(4.6, -3.0, "indefinite", fontsize=6.5, color=fs.MUTED, ha="left")
for key, lab, c, m, ls in ARMS:
    r = [x for x in dense if x["arm"] == key]
    ax.plot([float(x["eta"]) for x in r], [float(x["lmin"]) for x in r], color=c, marker=m, ls=ls, ms=3.4, label=lab)
eta_ticks(ax)
ax.set_xlabel(r"gold hardness $\eta_{\rm Au}$ (eV)"); ax.set_ylabel(r"$\lambda_{\min}$ (eV)")
ax.set_ylim(-4.2, 2.4); ax.axhline(0, color=fs.MUTED, lw=0.6)
fs.panel_label(ax, "a", dx=-0.19)

# (b) offset scan at eta = 5.172
ax = axs[0, 1]
ax.axhline(0, color=fs.MUTED, lw=0.6)
for key, lab, c, m, ls in ARMS:
    r = [x for x in off if x["arm"] == key and x["eta"] == "5.172"]
    if not r: continue
    x = [float(v["dchi"]) for v in r]
    ax.plot(x, [float(v["qwater"]) for v in r], color=c, marker=m, ls=ls, ms=3.2)
    ax.plot(x, [float(v["maxnet"]) for v in r], color=c, marker=m, ls="none", ms=3.2, mfc="white")
ax.set_xlabel(r"offset of the water electronegativities $\Delta\chi$ (eV)")
ax.set_ylabel(r"charge ($e$)")
ax.text(-0.25, 3.95, r"filled: $Q_{\rm water}$ (water $\to$ metal)" "\n" r"open: largest $|Q|$ of one molecule",
        fontsize=6.0, color=fs.MUTED, va="top")
fs.panel_label(ax, "b", dx=-0.17)

# (c) density-based screening vs height
ax = axs[1, 0]
ax.axhline(-1, color=fs.MUTED, lw=0.6, ls=(0, (2, 2)))
ax.text(2.2, -0.97, "ideal conductor", fontsize=6.3, color=fs.MUTED, ha="left", va="bottom")
st = {("gself", "0.3"): (fs.BLUE, "o", "solid"), ("gself", "5.172"): (fs.BLUE, "o", (0, (3, 1.5))),
      ("nogself", "0.3"): (fs.GREEN, "^", "solid"), ("nogself", "5.172"): (fs.GREEN, "^", (0, (3, 1.5)))}
for (arm, eta), (c, m, ls) in st.items():
    r = [x for x in e2 if x["arm"] == arm and x["eta"] == eta]
    lab = ("self-energy" if arm == "gself" else "no self-energy") + rf", $\eta$ = {float(eta):g} eV"
    ax.plot([float(x["z"]) for x in r], [float(x["Q_density_pair"]) for x in r], color=c, marker=m, ls=ls, ms=3.0,
            mfc=c if eta == "0.3" else "white", label=lab)
ax.set_xlabel("height of the unit charge $z$ (Å)"); ax.set_ylabel("charge above mid-plane ($e$)")
ax.set_ylim(-2.1, -0.85)
ax.legend(loc="center right", bbox_to_anchor=(1.0, 0.47), fontsize=5.4, handlelength=2.4, labelspacing=0.25)
fs.panel_label(ax, "c", dx=-0.24)

# (d) width scan
ax = axs[1, 1]
ax.axhline(0, color=fs.MUTED, lw=0.6)
ax.axhline(1, color=fs.MUTED, lw=0.6, ls=(0, (2, 2)))
for eta, mf in (("5.172", "white"), ("0.3", None)):
    r = [x for x in ws if x["eta"] == eta]
    w = [float(x["w"]) for x in r]
    ls = "solid" if eta == "0.3" else (0, (3, 1.5))
    ax.plot(w, [float(x["face_gauss"]) for x in r], color=fs.BLUE, marker="o", ms=3.0, mfc=mf or fs.BLUE, ls=ls,
            label=rf"face / Gauss, $\eta$ = {float(eta):g} eV")
    ax.plot(w, [float(x["L4"]) / float(x["L5"]) for x in r], color=fs.VERM, marker="s", ms=3.0, mfc=mf or fs.VERM,
            ls=ls, label=rf"sub-surface / face, $\eta$ = {float(eta):g} eV")
ax.axvline(2.378, color=fs.MUTED, lw=0.6, ls=(0, (1, 1.5)))
ax.text(2.33, 0.62, "own width", fontsize=6.0, color=fs.MUTED, ha="right")
ax.set_xlabel(r"self-energy width $w$ (Å)"); ax.set_ylabel("ratio")
ax.set_ylim(-1.1, 1.35)
ax.text(1.55, 1.25, "face / Gauss", fontsize=6.3, color=fs.BLUE, ha="center")
ax.text(0.33, -0.5, "sub-surface / face", fontsize=6.3, color=fs.VERM, ha="left")
ax.text(0.3, -1.0, r"solid: $\eta$ = 0.3 eV   dashed: $\eta$ = 5.172 eV", fontsize=6.0, color=fs.MUTED)
fs.panel_label(ax, "d", dx=-0.17)

h, l = axs[0, 0].get_legend_handles_labels()
fig.legend(h, l, loc="outside upper center", ncol=4, fontsize=6.5, handlelength=2.6)
fig.savefig(HERE / "fig_conductor.pdf")
fig.savefig(HERE / "fig_conductor.png", dpi=220)
print("wrote fig_conductor.pdf/.png")
