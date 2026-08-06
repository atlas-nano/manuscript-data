#!/usr/bin/env python3
"""Regenerate main-text Figures 5 and 6 and SI Figure S7 from the deposited data.

This script reads the deposited CSVs directly and assembles the full figures, so the deposit can
redraw its own figures with no external inputs.

REFERENCE MODE (--ref) selects which IAPWS-95 reference is used:

  original   iapws_sex_kb.csv as published. Reproduces the figures in the submitted manuscript,
             including the state points where the reference implementation silently flashed to a
             two-phase mixture and returned the saturated-liquid entropy.
  corrected  iapws_sex_kb_corrected.csv. The 12 points beyond the liquid spinodal are dropped, the
             9 metastable points are re-referenced onto the mechanically stable stretched-liquid
             branch, and the rest are unchanged. This is what the revised Section 4.5 describes.

Run --ref original first as a control: it must reproduce the published artwork before the corrected
version means anything.

Usage:
  python3 reproduce_figures.py --ref original  --outdir control
  python3 reproduce_figures.py --ref corrected --outdir corrected
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from matplotlib.cm import ScalarMappable

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)                     # the 2pt-excess-entropy/ deposit root
DATA = os.path.join(ROOT, "data", "entropy")
REFDIR = os.path.join(ROOT, "data", "reference")

MODELS = [("spce", "SPC/E"), ("tip4p-2005", "TIP4P/2005"),
          ("tip4p-ice", "TIP4P/ICE"), ("tip4p-2005f", "TIP4P/2005f")]

# the deposited filenames are not systematic, so they are spelled out rather than built
SEX = {"spce": "spce_excess_entropy.csv", "tip4p-2005": "tip4p-2005_excess_entropy.csv",
       "tip4p-ice": "tip4p-ice_excess_entropy.csv", "tip4p-2005f": "tip4p-2005f_excess_entropy.csv"}
S2 = {"spce": "spce_total_s2_total_entropy.csv", "tip4p-2005": "tip4p-2005_s2_total_entropy.csv",
      "tip4p-ice": "tip4p-ice_s2_total_entropy.csv",
      "tip4p-2005f": "tip4p-2005f_s2_total_entropy.csv"}
SEX_ADJ = {"spce": "spce_excess_entropy_adj.csv",
           "tip4p-2005": "tip4p-2005_excess_entropy_adj.csv",
           "tip4p-ice": "tip4p-ice_excess_entropy_adj.csv",
           "tip4p-2005f": "tip4p-2005f_excess_entropy_adj.csv"}
S2_ADJ = {"spce": "spce_total_s2_entropy_adj.csv", "tip4p-2005": "tip4p-2005_s2_entropy_adj.csv",
          "tip4p-ice": "tip4p-ice_s2_entropy_adj.csv",
          "tip4p-2005f": "tip4p-2005f_s2_entropy_adj.csv"}

ISOTHERMS = [220, 300, 360]
MAD_Z = 3.0          # same robust cutoff the original notebooks used


def _mad_mask(v):
    """Outliers within one (Density, Temperature) group, by robust z-score."""
    if v.size < 5 or np.allclose(v.std(ddof=1), 0.0):
        return np.zeros(v.size, dtype=bool)
    mad = np.median(np.abs(v - np.median(v)))
    if mad == 0:
        return np.zeros(v.size, dtype=bool)
    return np.abs(0.6745 * (v - np.median(v)) / mad) > MAD_Z


def agg(path):
    """Mean and standard deviation per state point, after MAD outlier removal."""
    df = pd.read_csv(path)
    df.columns = df.columns.str.strip()
    out = df.groupby(["Density", "Temperature"])["thermoValue"].transform(
        lambda s: pd.Series(_mad_mask(s.to_numpy()), index=s.index))
    kept = df[~out]
    return (kept.groupby(["Density", "Temperature"])["thermoValue"]
                .agg(mean="mean", std="std", n="count").reset_index())


def load_reference(mode):
    """Return the IAPWS reference, with dropped points removed when mode == 'corrected'."""
    if mode == "original":
        df = pd.read_csv(os.path.join(REFDIR, "iapws_sex_kb.csv"))
        df.columns = df.columns.str.strip()
        return df.rename(columns={"thermoValue": "ref"})[["Density", "Temperature", "ref"]]

    df = pd.read_csv(os.path.join(REFDIR, "iapws_sex_kb_corrected.csv"))
    df.columns = df.columns.str.strip()
    keep = ~df["status"].str.startswith("drop")
    print(f"  reference: dropped {(~keep).sum()} points beyond the spinodal, "
          f"kept {keep.sum()} of {len(df)}")
    n_re = (df["status"] == "re-referenced").sum()
    shift = (df.loc[df["status"] == "re-referenced", "thermoValue"]
             - df.loc[df["status"] == "re-referenced", "originalValue"]).abs().max()
    print(f"  reference: {n_re} points re-referenced, largest shift {shift:.3f} kB")
    return df[keep].rename(columns={"thermoValue": "ref"})[["Density", "Temperature", "ref"]]


# ----------------------------------------------------------------- Figure 5
def figure_5(ref, outdir):
    """4 models x 3 isotherms: 2PT S_ex, S2, and the IAPWS reference against density."""
    fig, axes = plt.subplots(4, 3, figsize=(9.0, 10.4), sharex=True, sharey=True)
    for i, (key, label) in enumerate(MODELS):
        sex = agg(os.path.join(DATA, "excess", SEX[key]))
        s2 = agg(os.path.join(DATA, "2body_total", S2[key]))
        for j, T in enumerate(ISOTHERMS):
            ax = axes[i, j]
            a, b = sex[sex.Temperature == T], s2[s2.Temperature == T]
            r = ref[ref.Temperature == T].sort_values("Density")
            ax.errorbar(a.Density, a["mean"], yerr=a["std"], fmt="s-", ms=4, lw=1.3,
                        capsize=2, color="tab:blue", label="$S_{ex}$")
            ax.errorbar(b.Density, b["mean"], yerr=b["std"], fmt="o-", ms=4, lw=1.3,
                        capsize=2, color="tab:orange", label="$S_2$")
            ax.plot(r.Density, r["ref"], "^--", ms=4, lw=1.3, color="tab:green",
                    label="$S_{ex}$ (IAPWS)")
            ax.grid(alpha=0.3)
            ax.set_ylim(-14, -4)
            if i == 0:
                ax.set_title(f"T = {T}K", fontsize=12)
            if j == 0:
                ax.set_ylabel("Entropy (S/$k_B$)", fontsize=10)
            if j == 2:
                ax.yaxis.set_label_position("right")
                ax.text(1.06, 0.5, label, transform=ax.transAxes, rotation=270,
                        va="center", ha="left", fontsize=11, fontweight="bold")
    # one x label for the whole grid, in its own band below the panels
    fig.supxlabel("Density (g/$cm^3$)", fontsize=12, y=0.055)
    h, l = axes[0, 0].get_legend_handles_labels()
    order = [l.index(k) for k in ("$S_{ex}$", "$S_2$", "$S_{ex}$ (IAPWS)")]
    h, l = [h[i] for i in order], [l[i] for i in order]
    fig.legend(h, l, loc="lower center", ncol=3, frameon=True, fontsize=11,
               bbox_to_anchor=(0.5, 0.005))
    fig.tight_layout(rect=[0, 0.075, 0.97, 1])
    p = os.path.join(outdir, "iapws_entropy_contributions_updated.png")
    fig.savefig(p, dpi=200)
    plt.close(fig)
    print(f"  wrote {p}")


# ----------------------------------------------------------------- Figure 6
def figure_6(ref, outdir):
    """4 models x 2 estimators: residual against the reference, on the offset-corrected values."""
    rows = []
    for key, label in MODELS:
        sex = agg(os.path.join(DATA, "excess", "adjusted", SEX_ADJ[key])).rename(
            columns={"mean": "sex"})
        s2 = agg(os.path.join(DATA, "2body_total", "adjusted", S2_ADJ[key])).rename(
            columns={"mean": "s2"})
        m = (ref.merge(sex[["Density", "Temperature", "sex"]], on=["Density", "Temperature"])
                .merge(s2[["Density", "Temperature", "s2"]], on=["Density", "Temperature"]))
        m["resid_sex"] = m["sex"] - m["ref"]
        m["resid_s2"] = m["s2"] - m["ref"]
        m["model"] = label
        rows.append(m)
    allm = pd.concat(rows)

    vmax = max(abs(allm["resid_sex"]).max(), abs(allm["resid_s2"]).max())
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
    # Dropped state points must not read as zero residual on a diverging scale, which is what an
    # unpainted (white) cell looks like. Paint them a neutral grey instead.
    cmap = matplotlib.colormaps["RdBu_r"].with_extremes(bad="0.82")
    full_grid = pd.MultiIndex.from_product(
        [sorted(allm.Temperature.unique()), sorted(allm.Density.unique())],
        names=["Temperature", "Density"])

    fig, axes = plt.subplots(4, 2, figsize=(7.2, 12.6))
    for i, (key, label) in enumerate(MODELS):
        m = allm[allm.model == label]
        for j, col in enumerate(["resid_sex", "resid_s2"]):
            ax = axes[i, j]
            piv = (m.set_index(["Temperature", "Density"])[col]
                    .reindex(full_grid)          # restore the dropped cells as NaN
                    .unstack("Density"))
            ax.pcolormesh(np.arange(piv.shape[1] + 1), np.arange(piv.shape[0] + 1),
                          np.ma.masked_invalid(piv.to_numpy()), cmap=cmap, norm=norm)
            ax.set_xticks(np.arange(piv.shape[1]) + 0.5)
            ax.set_xticklabels([f"{d:g}" if d in (0.9, 1.0, 1.1, 1.2, 1.3) else ""
                                for d in piv.columns], fontsize=9)
            ax.set_yticks(np.arange(piv.shape[0]) + 0.5)
            ax.set_yticklabels([f"{t:g}" for t in piv.index], fontsize=9)
            if i == 0:
                ax.set_title("2PT S$_{ex}$ - IAPWS-95" if j == 0 else "S$_2$ - IAPWS-95",
                             fontsize=11, fontweight="bold")
            if j == 0:
                ax.set_ylabel("Temperature (K)", fontsize=10)
            if j == 1:
                ax.text(1.04, 0.5, label, transform=ax.transAxes, rotation=270,
                        va="center", ha="left", fontsize=10, fontweight="bold")
            if i == 3:
                ax.set_xlabel("Density (g/$cm^3$)", fontsize=10)
    fig.tight_layout(rect=[0, 0.10, 0.95, 1])
    if allm.shape[0] // 4 < 70:
        fig.text(0.5, 0.082, "grey: no homogeneous-liquid reference exists "
                 "(state point beyond the liquid spinodal)", ha="center", fontsize=8.5)
    cax = fig.add_axes([0.28, 0.048, 0.44, 0.011])
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal")
    cb.set_label("Residual Entropy (S/$k_B$)", fontsize=10, fontweight="bold")
    p = os.path.join(outdir, "Residual_Entropy_updated.png")
    fig.savefig(p, dpi=200)
    plt.close(fig)
    print(f"  wrote {p}")

    stats = {"MAE_sex": allm.resid_sex.abs().mean(), "MAE_s2": allm.resid_s2.abs().mean(),
             "MaxAbs_sex": allm.resid_sex.abs().max(), "MaxAbs_s2": allm.resid_s2.abs().max(),
             "n_states": len(allm) // 4}
    print("  residual summary: " + "  ".join(f"{k}={v:.4g}" for k, v in stats.items()))
    return allm


# ----------------------------------------------------------------- SI Figure S7
def figure_s7(ref, outdir):
    """The IAPWS-95 excess entropy itself, isotherm by isotherm."""
    fig, ax = plt.subplots(figsize=(8, 6))
    for T in sorted(ref.Temperature.unique()):
        s = ref[ref.Temperature == T].sort_values("Density")
        ax.plot(s.Density, s["ref"], marker="o", label=f"T = {T:g} K")
    ax.set_xlabel("Density (g/$cm^3$)")          # published raster had these units inverted
    ax.set_ylabel("Excess Entropy ($S_{ex}/k_B$)")
    ax.grid(True)
    ax.legend(title="Isotherms", loc="lower right", fontsize=9)
    fig.tight_layout()
    p = os.path.join(outdir, "iapws_excess_entropy.png")
    fig.savefig(p, dpi=200)
    plt.close(fig)
    print(f"  wrote {p}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", choices=["original", "corrected"], default="corrected")
    ap.add_argument("--outdir", default=None)
    a = ap.parse_args()
    outdir = os.path.join(HERE, a.outdir or a.ref)
    os.makedirs(outdir, exist_ok=True)
    print(f"reference mode: {a.ref}  ->  {outdir}")
    ref = load_reference(a.ref)
    figure_5(ref, outdir)
    figure_6(ref, outdir)
    figure_s7(ref, outdir)


if __name__ == "__main__":
    main()
