"""Provides the panel drawing functions of Figure 3."""
from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

import fig_sources as SRC
from paths import data, panel_image
from pmf_panels import load as load_pmf

INK, INK2, MUTE, GRID = "#0b0b0b", "#52514e", "#8a8985", "#e6e5e1"
FS_TITLE, FS_LABEL, FS_TICK, FS_NOTE = 9.2, 8.2, 7.8, 6.9
SOLV_COL = {"FM": "#2a78d6", "ACN": "#eb6834", "Me₂O": "#4a3aa7", "DME": "#009E73", "DEE": "#6a7f8c"}
LI_COL, PAIR_COL, FREE_COL, BOX_COL = "#4a4a48", "#D55E00", "#56B4E9", "#3f3f96"

PMF_SRC = {"LiTPFA": ("litpfa_fm.npz", "pm_opes_analytic"),
           "LiPF6": ("lipf6_fm.npz", None),
           "LiBF4": ("libf4_fm.npz", None)}
PANEL_E_STRUCT = panel_image("fig3B_solvation_render.png")
PANEL_E_LEGEND = panel_image("fig3B_legend.png")

DGMIX_SYSTEMS = [(s, v) for v in ("FM", "ACN") for s in ("LiTPFA", "LiPF6", "LiBF4")]


def despine(ax, keep=("left", "bottom")):
    for sp in ("top", "right", "left", "bottom"):
        ax.spines[sp].set_visible(sp in keep)


def title(ax, txt, pad=4):
    ax.set_title(txt, loc="left", fontsize=FS_TITLE, color=INK, pad=pad)


def autocrop_white(path, thresh=248):
    from PIL import Image
    a = np.array(Image.open(path).convert("RGB"))
    ink = (a < thresh).any(axis=2)
    ys, xs = np.where(ink)
    if not len(xs):
        return a
    pad = int(round(0.02 * max(xs.max() - xs.min(), ys.max() - ys.min())))
    return a[max(0, ys.min() - pad):ys.max() + 1 + pad, max(0, xs.min() - pad):xs.max() + 1 + pad]


def autocrop_alpha(path):
    from PIL import Image
    a = np.array(Image.open(path).convert("RGBA"))
    ys, xs = np.where(a[..., 3] > 8)
    if not len(xs):
        return a
    return a[ys.min():ys.max() + 1, xs.min():xs.max() + 1]


def panel_a(ax, ax_struct, ax_key, report):
    be = SRC.BE_TABLE_S8
    order = sorted(be, key=lambda s: be[s][0], reverse=True)
    xs = np.arange(len(order))
    ax.axhline(0, color=MUTE, lw=0.8, ls=":", zorder=1)
    ax.grid(True, axis="y", color=GRID, lw=0.5, zorder=0)
    for x, s in zip(xs, order):
        ax.bar(x, be[s][0], width=0.62, color=SOLV_COL[s], zorder=3)
        for v in be[s][1:]:
            if v is not None:
                ax.plot([x], [v], marker="o", ms=2.6, mfc="white", mec=INK2, mew=0.7, zorder=5)
    ax.set_xticks(xs)
    ax.set_xticklabels(order, fontsize=FS_TICK)
    ax.set_ylim(-58, 0)
    ax.set_ylabel("Li⁺–solvent binding\nenergy (kcal mol⁻¹)", fontsize=FS_LABEL)
    ax.tick_params(labelsize=FS_TICK)
    ax.text(0.025, 0.05, "bars n = 1;  ○ n = 2–4", transform=ax.transAxes, ha="left", va="bottom",
            fontsize=FS_NOTE, color=INK2)
    despine(ax)

    if PANEL_E_STRUCT is not None:
        ax_struct.imshow(autocrop_white(PANEL_E_STRUCT), interpolation="bilinear")
    if PANEL_E_LEGEND is not None:
        ax_key.imshow(autocrop_alpha(PANEL_E_LEGEND), interpolation="bilinear")
    for a in (ax_struct, ax_key):
        a.set_axis_off()
    ax_struct.text(0.5, -0.035, "Li⁺ in FM, dilute limit", transform=ax_struct.transAxes,
                   ha="center", va="top", fontsize=FS_NOTE, color=INK2)


def _snapshot_stats():
    import json
    return json.load(open(data("simulation", "snapshot_stats.json")))


def common_crop(paths, thresh=246, pad_frac=0.012):
    from PIL import Image
    arr = [np.array(Image.open(p).convert("RGB")) for p in paths]
    x0, y0, x1, y1 = 10 ** 9, 10 ** 9, -1, -1
    for a in arr:
        ys, xs = np.where((a < thresh).any(axis=2))
        x0, x1 = min(x0, int(xs.min())), max(x1, int(xs.max()))
        y0, y1 = min(y0, int(ys.min())), max(y1, int(ys.max()))
    pad = int(round(pad_frac * max(x1 - x0, y1 - y0)))
    H, W = arr[0].shape[:2]
    y0, y1 = max(0, y0 - pad), min(H - 1, y1 + pad)
    x0, x1 = max(0, x0 - pad), min(W - 1, x1 + pad)
    return [a[y0:y1 + 1, x0:x1 + 1] for a in arr]


def panel_b(ax_key, ax_left, ax_right, report):
    st = _snapshot_stats()
    order = [("litpfa", "LiTPFA", ax_left), ("lipf6", "LiPF$_6$", ax_right)]
    pngs = [panel_image(f"fig3C_{s}_fm_cell_render.png") for s, _, _ in order]
    imgs = common_crop(pngs) if all(p is not None for p in pngs) else [None] * len(pngs)

    for (salt, lab, ax), im in zip(order, imgs):
        d = st["systems"][salt]
        if im is not None:
            ax.imshow(im, interpolation="bilinear")
        ax.set_axis_off()
        npair = d["n_paired_anions"]
        ax.text(0.5, -0.015, lab, transform=ax.transAxes, ha="center", va="top",
                fontsize=FS_NOTE + 0.8, color=INK, fontweight="bold")
        ax.text(0.5, -0.115, f"{d['selection']['nsp']} species · {npair}/{d['n_pairs']} paired",
                transform=ax.transAxes, ha="center", va="top", fontsize=FS_NOTE, color=INK2)

    from matplotlib.lines import Line2D
    ax_key.set_axis_off()
    ax_key.set_xlim(0, 1)
    ax_key.set_ylim(0, 1)
    handles = [Line2D([], [], marker="o", ms=4.6, mfc=c, mec=INK, mew=0.35, ls="none", label=l)
               for c, l in ((LI_COL, "Li$^+$"), (PAIR_COL, "Li$^+$-contacted anion"),
                            (FREE_COL, "free anion"))]
    leg = ax_key.legend(handles=handles, ncol=3, loc="upper left", frameon=False,
                        fontsize=FS_NOTE, handlelength=0.8, handletextpad=0.35,
                        columnspacing=0.9, borderpad=0.0, borderaxespad=0.0)
    note = ax_key.text(0.0, 0.02, "8 ion pairs in FM, 293 K · solvent omitted",
                       ha="left", va="bottom", fontsize=FS_NOTE - 0.5, color=INK2)


def panel_d(ax, report):
    ylo, yhi = -56, 40
    ax.axvspan(10, 12, color=GRID, alpha=0.8, lw=0, zorder=0)
    ax.axhline(0, color=MUTE, lw=0.8, ls=":", zorder=1)
    ax.grid(True, axis="y", color=GRID, lw=0.5, zorder=1)
    ax.text(11.0, -22.0, "", ha="center", va="center", fontsize=FS_NOTE,
            color=INK2, linespacing=1.15)
    for salt in ("LiTPFA", "LiPF6", "LiBF4"):
        f, pre = PMF_SRC[salt]
        cen, m, s, _ = load_pmf(data("simulation", "opes_pmf", f), pre)
        k = np.isfinite(m) & np.isfinite(s) & (cen <= 12.05)
        ax.fill_between(cen[k], (m - s)[k], (m + s)[k], color=SRC.SALT_COL[salt], alpha=0.22,
                        lw=0, zorder=3)
        ax.plot(cen[k], m[k], color=SRC.SALT_COL[salt], lw=1.7, zorder=5)

    ann = {"LiTPFA": (4.4, 30.5), "LiPF6": (4.4, -35.0)}
    for salt, (x, y) in ann.items():
        g, sd = SRC.DG_ASSOC_FM[salt]
        txt = (f"{SRC.LAB[salt]}\n$\\Delta G_{{assoc}}$ {g:+.1f}"
               + (f" ± {sd:.1f}" if sd else "") + " k$_B$T")
        a = ax.text(x, y, txt, fontsize=FS_NOTE, color=SRC.SALT_COL[salt], fontweight="bold",
                    ha="left", va="center", zorder=6, linespacing=1.2)
    ax.annotate("", xy=(4.25, SRC.DG_ASSOC_FM["LiTPFA"][0]), xytext=(4.25, SRC.DG_ASSOC_FM["LiPF6"][0]),
                arrowprops=dict(arrowstyle="<->", color=INK2, lw=0.9, shrinkA=0, shrinkB=0), zorder=6)
    ax.text(4.45, -20.0, f"$\\Delta\\Delta G$  {SRC.DDG_ASSOC_FM:+.1f} k$_B$T", ha="left",
            va="center", fontsize=FS_NOTE, color=INK2, zorder=6,
            bbox=dict(facecolor="white", edgecolor="none", pad=1.2))
    ax.text(0.985, 0.035, "no $\\Delta G_{assoc}$ is reported for LiBF$_4$", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=FS_NOTE, color=SRC.SALT_COL["LiBF4"])
    ax.set_xlim(1.5, 12.1)
    ax.set_ylim(ylo, yhi)
    ax.set_xticks([2, 4, 6, 8, 10, 12])
    ax.set_xlabel("Li⁺–anion centre-of-mass\ndistance $r$ (Å)", fontsize=FS_LABEL,
                  linespacing=1.15)
    ax.set_ylabel("$W(r)$ (k$_B$T)", fontsize=FS_LABEL)
    ax.tick_params(labelsize=FS_TICK)
    despine(ax)


def panel_e(ax, report):
    BOX = dict(facecolor="white", edgecolor="none", pad=1.0)
    diss, tds = SRC.gated_diss_tds()
    rows = SRC.dgmix_systems(DGMIX_SYSTEMS)
    xs = np.arange(len(DGMIX_SYSTEMS))
    for i, (salt, solv) in enumerate(DGMIX_SYSTEMS):
        v = rows[(salt, solv)]
        ax.bar(i, diss[salt], width=0.62, color="#d97b66", zorder=3,
               label="$\\Delta G_{diss}$" if i == 0 else None)
        ax.bar(i, v["gsolv"], width=0.62, color="#4d86c6", zorder=3,
               label="$\\Delta G_{solv}$" if i == 0 else None)
        ax.bar(i, -tds[solv], width=0.62, bottom=v["gsolv"], color=MUTE, zorder=3,
               label="$-T\\Delta S_{ideal}$" if i == 0 else None)
        ax.plot([i], [v["gmix"]], marker="D", ms=4.4, color=INK, ls="none", zorder=6,
                label="net $\\Delta G_{mix}$" if i == 0 else None)
        if solv == "FM" and salt == "LiTPFA":
            lo, hi = sorted((v["gmix"], v["gmix_raw"]))
            ax.add_patch(Rectangle((i - 0.31, lo), 0.62, hi - lo, facecolor=INK, alpha=0.20,
                                   lw=0, zorder=5))
            ax.annotate("≈ 0", (i, v["gmix"]), textcoords="offset points", xytext=(0, -10),
                        ha="center", va="top", fontsize=FS_NOTE, color=INK, zorder=7, bbox=BOX)
        else:
            ax.annotate(f"{v['gmix']:+.0f}", (i, v["gmix"]), textcoords="offset points",
                        xytext=(0, 10 if v["gmix"] >= 0 else -10), ha="center",
                        va="bottom" if v["gmix"] >= 0 else "top", fontsize=FS_NOTE, color=INK,
                        zorder=7, bbox=BOX)
    ax.axhline(0, color=INK, lw=1.0, zorder=4)
    ax.axvline(2.5, color=INK2, lw=0.7, ls=":", zorder=2)
    ax.set_xticks(xs)
    ax.set_xticklabels([SRC.LAB[s] for s, _ in DGMIX_SYSTEMS], fontsize=FS_NOTE, rotation=35,
                       ha="right", rotation_mode="anchor")
    for xf, lab in ((1.0, "in FM"), (4.0, "in ACN")):
        ax.text((xf + 0.7) / 6.0, -0.295, lab, transform=ax.transAxes, ha="center", va="top",
                fontsize=FS_NOTE + 0.5, color=INK2, fontweight="bold")
    ax.set_xlim(-0.7, len(DGMIX_SYSTEMS) - 0.3)
    ax.set_ylim(-268, 206)
    ax.set_yticks([-150, -100, -50, 0, 50, 100, 150])
    ax.set_ylabel("kcal mol⁻¹ (per formula unit)", fontsize=FS_LABEL)
    ax.tick_params(labelsize=FS_TICK)
    ax.grid(True, axis="y", color=GRID, lw=0.5, zorder=0)
    ax.legend(fontsize=FS_NOTE, frameon=False, loc="lower left", ncol=2, handlelength=1.0,
              columnspacing=0.8, handletextpad=0.4, borderpad=0.0, labelspacing=0.22)
    despine(ax)
