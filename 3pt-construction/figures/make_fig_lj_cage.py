#!/usr/bin/env python
"""Regenerate fig1_error_summary and fig2_signed_error for the 3PT/cage paper
from the full 71-point LJ grid (lj_fullgrid_cage.csv, this work).

fig1  signed mean bias +/- std by method over the 43 homogeneous-liquid points,
      with RMS annotated -> rig-HS low, Lin-2003 high, R2PT tuned, parameter-free 3PT ~0.
fig2  point-by-point signed dS* = S*_method - S*_MBWR vs rho* over the liquid
      grid -> shows the systematic bias is removed by the cage at every state.
"""
import csv, math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
CSV  = DATA / "lj_fullgrid_cage.csv"

# five entropy functionals (rig-HS low, Lin-2003 high, Desjarlais MF gas, R2PT tuned,
# 3PT parameter-free) — the same five compared in the theory section.
METHODS = [("d_rigHS", "rig-HS\n(uncorrected)", "#8172b3"),
           ("d_Lin",   "Lin-2003\n($+\\ln Z$)", "#c44e52"),
           ("d_Desj",  "Desjarlais\n(MF gas)", "#937860"),
           ("d_R2PT",  "R2PT\n($\\delta{=}1.5$)", "#dd8452"),
           ("d_3PT",   "3PT\n(parameter-free)", "#55a868")]

plt.rcParams.update({"font.size": 10, "axes.linewidth": 0.8,
                     "figure.dpi": 150, "savefig.bbox": "tight"})


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(HERE / f"{name}.{ext}")
    plt.close(fig)
    print("  wrote", name + ".png/.pdf")


def load():
    rows = list(csv.DictReader(open(CSV)))
    for r in rows:
        for k in r:
            if k == "ref_source":      # string column; leave as-is
                continue
            r[k] = float(r[k]) if r[k] != "" else float("nan")   # R2PT blank at near-solid rho*>1.02
    # homogeneous-liquid regime: well-defined diffusive gas fraction
    liq = [r for r in rows if 0.70 <= r["rho_star"] <= 1.02]
    return rows, liq


def rms(xs):  return math.sqrt(sum(x * x for x in xs) / len(xs))


def fig_bias(liq):
    fig, ax = plt.subplots(figsize=(6.2, 3.7))
    xs = np.arange(len(METHODS))
    for i, (key, lab, col) in enumerate(METHODS):
        d = [r[key] for r in liq]
        m, s, rr = np.mean(d), np.std(d), rms(d)
        ax.bar(i, m, yerr=s, color=col, width=0.62, capsize=4,
               error_kw=dict(lw=1.0, ecolor="0.3"), zorder=3)
        va = "bottom" if m >= 0 else "top"
        off = 0.012 if m >= 0 else -0.012
        ax.text(i, m + np.sign(m) * s + off, f"RMS {rr:.2f}", ha="center",
                va=va, fontsize=9, fontweight="bold")
    ax.axhline(0, color="k", lw=0.8, zorder=2)
    ax.set_xticks(xs)
    ax.set_xticklabels([m[1] for m in METHODS], fontsize=9)
    ax.set_ylabel(r"signed mean $\langle \Delta S^* \rangle$ vs MBWR  ($k_B$/atom)")
    ax.set_title(f"Bias by method, {len(liq)} homogeneous-liquid states\n"
                 r"($0.70\leq\rho^*\leq1.02$); error bars $=\pm$ std",
                 fontsize=9.5)
    ax.set_ylim(-0.45, 0.30)
    ax.grid(axis="y", ls=":", lw=0.5, alpha=0.6)
    fig.tight_layout()
    save(fig, "fig1_error_summary")


def fig_signed(liq, legpos="top", name="fig2_signed_error"):
    """legpos: 'top' or 'bottom' -> horizontal legend (ncol=3) above/below the axes,
    leaving the data field full-width (no dead whitespace)."""
    fig, ax = plt.subplots(figsize=(6.0, 4.0))
    liq = sorted(liq, key=lambda r: r["rho_star"])
    for key, lab, col in METHODS:
        rho = [r["rho_star"] for r in liq]
        d   = [r[key] for r in liq]
        ax.plot(rho, d, "o", ms=5, color=col, alpha=0.85,
                label=lab.replace("\n", " "), zorder=3)
    ax.axhline(0, color="k", lw=0.9, zorder=2)
    ax.axhspan(-0.05, 0.05, color="0.85", alpha=0.5, zorder=1,
               label=r"$\pm0.05\,k_B$")
    ax.set_xlabel(r"reduced density  $\rho^*$")
    ax.set_ylabel(r"signed $\Delta S^* = S^*_{\rm method} - S^*_{\rm MBWR}$  ($k_B$/atom)")
    ax.set_ylim(-0.40, 0.22)
    # horizontal legend (3 cols x 2 rows) parked OUTSIDE the axes; savefig bbox="tight"
    # (rcParams) crops to include it.  Title is placed on the opposite side from the legend.
    if legpos == "top":
        ax.legend(fontsize=8, loc="lower center", bbox_to_anchor=(0.5, 1.01),
                  ncol=3, framealpha=0.9, columnspacing=1.3, handletextpad=0.4)
        ax.set_xlabel(r"reduced density  $\rho^*$")  # title omitted; legend occupies the top
    else:  # bottom
        ax.set_title("Per-state deviation across the liquid grid", fontsize=9.5)
        ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.16),
                  ncol=3, framealpha=0.9, columnspacing=1.3, handletextpad=0.4)
    ax.grid(ls=":", lw=0.5, alpha=0.6)
    fig.tight_layout()
    save(fig, name)


if __name__ == "__main__":
    rows, liq = load()
    print(f"loaded {len(rows)} grid points, {len(liq)} liquid")
    for key, lab, _ in METHODS:
        d = [r[key] for r in liq]
        print(f"  {lab.splitlines()[0]:10s} RMS={rms(d):.3f}  mean={np.mean(d):+.3f}")
    fig_bias(liq)
    fig_signed(liq, legpos="top",    name="fig2_signed_error")         # canonical (top legend)
    fig_signed(liq, legpos="bottom", name="fig2_signed_error_legbot")  # alternative layout
