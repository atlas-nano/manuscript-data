"""Figure 5: cage efficiency against the kernel non-Markovianity for Lennard-Jones, metal and water systems. Reads ../data/cage_efficiency.json."""
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
JSN = os.path.join(HERE, os.pardir, "data", "cage_efficiency.json")
OUT = os.path.join(HERE, "out")
os.makedirs(OUT, exist_ok=True)

P = json.load(open(JSN))
DEF = P["deficit"]
RES = P["resolution"]
LABEL = {"tip4p2005": "TIP4P/2005", "spcew": "SPC/E"}

fig, ax = plt.subplots(figsize=(8, 6.4))
xs = np.linspace(0.0, 5.2, 300)
ax.plot(xs, 1 - np.exp(-P["fit_k"] * xs), "k-", lw=1.6,
        label=f"η = 1−e^(−{P['fit_k']:.2f}γ/Ω₀)  (resolved LJ + metals, R²={P['fit_R2']:.2f})")

lad = {}
for p in P["ladder"]:
    lad.setdefault(p["system"], []).append(p)
for s, v in lad.items():
    v = sorted(v, key=lambda q: q["dt_fs"])
    nres = sum(q["dt_Omega0"] <= RES for q in v)
    v = v[:min(len(v), nres + 1)]
    c = "crimson" if v[0]["cls"] == "metals" else "tab:blue"
    x = [q["x"] for q in v]
    e = [q["eta"] for q in v]
    ax.plot(x, e, "-", color=c, lw=0.8, alpha=0.55, zorder=2)
    for q in v:
        ok = q["dt_Omega0"] <= RES
        ax.scatter(q["x"], q["eta"], s=10 + 4 * np.log2(q["dt_fs"] + 1) ** 2, marker="o",
                   facecolor=c if ok else "white", edgecolor=c, lw=0.7, alpha=0.8, zorder=3)
    ax.annotate(f"{v[-1]['dt_fs']} fs", (x[-1], e[-1]), fontsize=6, color=c, xytext=(3, -8),
                textcoords="offset points")
ax.plot([], [], "-o", color="crimson", ms=4, mfc="white", lw=0.8,
        label="Δt ladder, metals (resolved spacings + first unresolved, hollow)")
ax.plot([], [], "-o", color="tab:blue", ms=4, mfc="white", lw=0.8,
        label="Δt ladder, LJ (resolved spacings + first unresolved, hollow)")

okp = [p for p in P["lj"] if p["dt_Omega0"] <= RES]
bad = [p for p in P["lj"] if p["dt_Omega0"] > RES]
ax.scatter([p["x"] for p in okp], [p["eta"] for p in okp], c="tab:blue", marker="s", s=30, alpha=.75,
           edgecolor="k", lw=.2, zorder=4, label="LJ 8 fs (Δt·Ω₀ ≤ 0.066)")
ax.scatter([p["x"] for p in bad], [p["eta"] for p in bad], facecolor="white", marker="s", s=30,
           edgecolor="tab:blue", lw=.8, zorder=4, label="LJ 8 fs (Δt·Ω₀ > 0.066)")
MET = P["metals"]
ax.scatter([p["x"] for p in MET], [p["eta"] for p in MET], c="crimson", marker="D", s=55, edgecolor="k",
           lw=.3, zorder=5, label="metals (TI), 1 fs")
for p in MET:
    ax.annotate(p["name"], (p["x"], p["eta"]), fontsize=6.5, xytext=(3, 2), textcoords="offset points")

for w in P["water"]:
    lab = LABEL[w["model"]]
    ax.scatter(w["x_trans"], w["eta_trans"], c="teal", marker="*", s=260, edgecolor="k", lw=.6, zorder=6)
    ax.annotate(lab + " trans", (w["x_trans"], w["eta_trans"]), xytext=(8, -12), textcoords="offset points",
                fontsize=7.5)
    ax.errorbar(w["x_total"], w["eta_total"],
                yerr=[[w["eta_total"] - w["eta_total_low"]], [w["eta_total_high"] - w["eta_total"]]],
                fmt="P", color="purple", ms=10, capsize=3, zorder=6, mec="k", mew=.5)
    ax.annotate(lab + " total", (w["x_total"], w["eta_total"]), xytext=(8, -3), textcoords="offset points",
                fontsize=7.5)
ax.scatter([], [], c="teal", marker="*", s=160, edgecolor="k", label=f"water trans (ref-free, deficit={DEF:.3f})")
ax.errorbar([], [], fmt="P", color="purple", label="water total (FEP ref ± σ)")
ax.axhline(1, color="0.6", ls=":", lw=1)
ax.set_xlabel(r"non-Markovianity  $\gamma/\Omega_0$")
ax.set_ylabel(r"cage efficiency  $\Delta S_{cage}/$deficit")
ax.set_xlim(-0.1, 5.2)
ax.set_ylim(0.0, 1.5)
ax.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.11), ncol=2, framealpha=0.9)
ax.grid(alpha=.3)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig_universal_nonmarkovian.pdf"))
