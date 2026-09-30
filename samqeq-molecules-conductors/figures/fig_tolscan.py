"""Figure S1: the solved liquid dipole and the solve cost against the charge-solve tolerance. The data are listed below."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import figstyle
from figstyle import BLUE, VERM, INK, MUTED

figstyle.apply()

DATA = [
    (1.0e-1,  2, 2.715739812, 6.66e-02),
    (1.0e-2,  5, 2.798166254, 5.06e-03),
    (1.0e-3,  8, 2.798719915, 5.09e-04),
    (1.0e-4, 11, 2.798709480, 5.40e-05),
    (1.0e-5, 13, 2.798710246, 8.90e-06),
    (1.0e-6, 16, 2.798710179, 8.91e-07),
    (1.0e-7, 19, 2.798710187, 8.76e-08),
    (1.0e-8, 22, 2.798710187, 8.92e-09),
    (1.0e-9, 25, 2.798710187, 7.16e-10),
    (1.0e-10, 28, 2.798710187, 5.33e-11),
    (1.0e-11, 31, 2.798710187, 3.55e-12),
]

tol = np.array([d[0] for d in DATA])
mv = np.array([d[1] for d in DATA], dtype=float)
dip = np.array([d[2] for d in DATA])
conv = dip[-1]
DECK = 1.0e-5

fig, (axa, axb) = plt.subplots(
    1, 2, figsize=(figstyle.TEXTWIDTH, 2.35), gridspec_kw={"wspace": 0.30})

axa.semilogx(tol[1:], dip[1:], "o-", color=BLUE, zorder=3)
axa.semilogx(tol[:1], dip[:1], "o", color=VERM, zorder=4)
axa.axhline(conv, color=MUTED, lw=0.5, ls=":", zorder=1)
axa.invert_xaxis()
axa.set_xlabel(r"charge-solve tolerance $\tau$")
axa.set_ylabel(r"mean molecular dipole $\mu$  (D)")
axa.set_ylim(2.70, 2.83)
figstyle.panel_label(axa, "a", dx=-0.17)

axa.annotate("converged label,\ndipole 3.0% low",
             xy=(1.0e-1, dip[0]), xytext=(4.0e-3, 2.725),
             fontsize=6.5, color=VERM, ha="left", va="center",
             arrowprops=dict(arrowstyle="->", color=VERM, lw=0.7,
                             shrinkA=1, shrinkB=3))
axa.annotate(r"flat to $10^{-5}$ D from $10^{-3}$", xy=(1e-7, conv), xytext=(1e-6, 2.765),
             fontsize=6.5, color=INK, ha="center", va="center",
             arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.5,
                             shrinkA=1, shrinkB=2))

axb.semilogx(tol, mv, "s-", color=VERM, ms=3.4, zorder=3)
axb.invert_xaxis()
axb.set_xlabel(r"charge-solve tolerance $\tau$")
axb.set_ylabel("Krylov iterations")
axb.set_ylim(0, 34)
figstyle.panel_label(axb, "b", dx=-0.17)

for ax in (axa, axb):
    ax.axvline(DECK, color=INK, lw=0.6, ls="--", alpha=0.55, zorder=2)
axb.text(DECK, 33.0, "shipped decks", fontsize=6.5, color=INK,
         ha="center", va="top",
         bbox=dict(fc="white", ec="none", pad=1.0))

axb.annotate("", xy=(1.0e-11, 3.0), xytext=(1.0e-5, 3.0),
             arrowprops=dict(arrowstyle="<->", color=MUTED, lw=0.7))
axb.text(3.0e-8, 5.0, "+18 iterations,\nno change in $\\mu$",
         fontsize=6.5, color=MUTED, ha="center", va="bottom")

fig.savefig("fig_tolscan.pdf")
fig.savefig("fig_tolscan.png", dpi=200)
print("wrote fig_tolscan.pdf / .png")
print(f"converged reference mu = {conv:.9f} D")
print(f"max |dev| over tau <= 1e-4 : {max(abs(dip[i]-conv) for i in range(3, len(dip))):.2e} D")
print(f"deviation at 1e-1          : {abs(dip[0]-conv):.4f} D "
      f"({100*abs(dip[0]-conv)/conv:.1f}%)")
