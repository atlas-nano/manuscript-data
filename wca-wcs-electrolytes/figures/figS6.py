"""Composes Figure S6 from the three electrophoretic NMR spectra in data/experimental/enmr."""
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np

from paths import EXP, OUT as OUTDIR

SRC = EXP / "enmr"
OUT = OUTDIR / "figS6"

PANELS = [("enmr_1H.png", "$^{1}$H", "solvent", "$v$ = −4.9 ± 0.8 µm s$^{-1}$"),
          ("enmr_7Li.png", "$^{7}$Li", "Li$^+$",  "$v$ = +25.2 ± 1.7 µm s$^{-1}$"),
          ("enmr_19F.png", "$^{19}$F", "TPFA$^-$", "$v$ = −24.6 ± 1.4 µm s$^{-1}$")]


def trim(a, pad=6):
    if a.shape[2] == 4:
        al = a[..., 3:4]
        a = a[..., :3] * al + (1.0 - al)
    ink = (a.sum(axis=2) < 2.97)
    rows, cols = np.where(ink.any(axis=1))[0], np.where(ink.any(axis=0))[0]
    r0, r1 = max(rows[0] - pad, 0), min(rows[-1] + pad + 1, a.shape[0])
    c0, c1 = max(cols[0] - pad, 0), min(cols[-1] + pad + 1, a.shape[1])
    return a[r0:r1, c0:c1]


def main():
    fig = plt.figure(figsize=(6.9, 2.62))
    for k, (img, nuc, species, vel) in enumerate(PANELS):
        ax = fig.add_axes([0.018 + k * 0.3283, 0.035, 0.305, 0.775])
        ax.imshow(trim(mpimg.imread(SRC / img)), interpolation="antialiased")
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_title(f"({'abc'[k]})  {nuc}   {species}", fontsize=8.4, loc="left", pad=4.5)
        ax.text(0.5, -0.055, vel, transform=ax.transAxes, ha="center", va="top", fontsize=7.6)
    fig.text(0.5, 0.955,
             "Electrophoretic NMR of 0.25 M LiTPFA in acetonitrile · phase shift against "
             "applied voltage, 1 V mm$^{-1}$",
             ha="center", va="top", fontsize=8.6)
    fig.savefig(f"{OUT}.png", dpi=300)
    print("figS6.png written")


if __name__ == "__main__":
    main()
