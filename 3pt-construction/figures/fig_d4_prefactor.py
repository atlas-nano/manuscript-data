#!/usr/bin/env python3
"""SI figure: the d=4 cage-prefactor test (ms2a §III.F).

(a) The non-memory (soft-mode anharmonic) entropy per mode that each candidate
    prefactor p∈{1/4,1/3,1/2} REQUIRES to land the cage efficiency η=ΔS_cage(p)/
    (S_TI−S_rigHS) on the 3D-calibrated universal curve η=1.050(1−e^{−2.386·𝓜}),
    plotted vs the memory-area metric 𝓜.  p=1/4 leaves a flat residual at the
    per-mode anharmonic (liquid-over-crystal) scale ~0.1-0.2 k_B/mode (shaded);
    p=1/2 falls below it with M-dependent scatter.
(b) Distribution of that residual per prefactor (box = IQR, whiskers = range),
    against the anharmonic scale.

Data: 4D-LJ S_q runs (T*=3.0 & 4.0 isotherms, P*=3.2 isobar, decisive point),
analyzed by the d=4 py-xPT engine.
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data" / "lj_4d"

D = 4
ETA = lambda M: 1.050 * (1.0 - np.exp(-2.386 * M))   # 3D-calibrated universal curve
BAND = (0.10, 0.23)                                  # per-mode anharmonic (liquid-over-crystal) scale
PS = [(0.25, "p=1/4", "C3", "o"), (1/3, "p=1/3", "C0", "s"), (0.5, "p=1/2", "C7", "^")]

def collect():
    """Return (memarea, deficit, cage14) for cage-active fluid 4D states."""
    M, defi, cg = [], [], []
    src = [(str(DATA / "isotherm_T3_result.npz"), "iso"),
           (str(DATA / "isotherm_T4_result.npz"), "iso"),
           (str(DATA / "isobar_P32_result.npz"), "bar")]
    for f, kind in src:
        d = np.load(f); o = d['rows']; ci = {c: i for i, c in enumerate(str(d['cols']).split(','))}
        for r in o:
            mem = r[ci['memarea']]; Srig = r[ci['S_rigHS']]; STI = r[ci['S_TI']]
            c14 = (r[ci['S_3PT']] - Srig) if 'S_3PT' in ci else r[ci['dScage14']]
            de = STI - Srig
            if not np.isfinite(mem) or de <= 0.05 or c14 <= 0.05 or r[ci['rho']] < 0.6:
                continue
            M.append(mem); defi.append(de); cg.append(c14)
    # The decisive point (rho*=1.10, T*=1.0) requires the 4D production
    # trajectory, which is not shipped in this repository (see README). It is
    # omitted here; the figure regenerates from the isotherm/isobar result files.
    return np.array(M), np.array(defi), np.array(cg)

M, defi, cg = collect()
order = np.argsort(M)
M, defi, cg = M[order], defi[order], cg[order]
ec = ETA(M)

fig, ax = plt.subplots(1, 2, figsize=(9.2, 3.8), gridspec_kw=dict(width_ratios=[2.0, 1.0]))

# ── (a) required soft-mode/mode vs memarea ───────────────────────────────────
ax[0].axhspan(*BAND, color="0.85", zorder=0, label="anharmonic scale ($\\sim$0.1$-$0.2 $k_B$/mode)")
for p, lab, col, mk in PS:
    sm = (defi - (4.0*p)*cg/ec) / D
    ax[0].plot(M, sm, mk+"-", color=col, ms=4.5, lw=1.0, label=lab)
ax[0].axhline(0, color="k", lw=0.5)
ax[0].set_xlabel(r"memory area $\mathcal{M}=\int|F_K-F_M|/\int|F_K|$")
ax[0].set_ylabel(r"required soft-mode $\Delta s_{\rm anh}$ / mode  [$k_B$]")
ax[0].set_title(r"(a) non-memory residual each $p$ requires", fontsize=10)
ax[0].legend(fontsize=8, loc="upper left", framealpha=0.9)
ax[0].set_ylim(-0.02, 0.26)

# ── (b) distribution per prefactor ───────────────────────────────────────────
ax[1].axhspan(*BAND, color="0.85", zorder=0)
data, labels, colors = [], [], []
for p, lab, col, mk in PS:
    sm = (defi - (4.0*p)*cg/ec) / D
    data.append(sm); labels.append(lab); colors.append(col)
bp = ax[1].boxplot(data, labels=labels, widths=0.6, patch_artist=True, showfliers=False)
for patch, col in zip(bp['boxes'], colors):
    patch.set_facecolor(col); patch.set_alpha(0.45)
for med in bp['medians']:
    med.set_color("k")
ax[1].axhline(0, color="k", lw=0.5)
ax[1].set_ylabel(r"$\Delta s_{\rm anh}$ / mode  [$k_B$]")
ax[1].set_title("(b) selected: $p=1/4$", fontsize=10)
ax[1].set_ylim(-0.02, 0.26)

fig.tight_layout()
out = str(HERE / "fig_d4_prefactor")
for ext in ("pdf", "png"):
    fig.savefig(f"{out}.{ext}", dpi=200, bbox_inches="tight")
print(f"saved {out}.pdf/.png   (n={len(M)} 4D states; memarea {M.min():.2f}–{M.max():.2f})")
for p, lab, col, mk in PS:
    sm = (defi - (4.0*p)*cg/ec) / D
    inb = np.mean((sm >= BAND[0]) & (sm <= BAND[1]))*100
    print(f"  {lab}: mean {sm.mean():.3f}  CV {100*sm.std()/abs(sm.mean()):.0f}%  in-band {inb:.0f}%")
