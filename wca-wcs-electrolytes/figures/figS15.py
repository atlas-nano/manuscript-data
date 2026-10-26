"""Draws Figure S15 from data/experimental/figS15A.csv, figS15B.csv and figS15C.csv."""
import pandas as pd
from figstyle import *

REG = [("C1s", "C 1s"), ("O1s", "O 1s"), ("F1s", "F 1s"), ("Al2p", "Al 2p")]
ROWS = [("figS15A", "A", "Cu electrode,\nafter 50 Li\nplating/stripping\ncycles"),
        ("figS15B", "B", "Li metal,\nsoaked in\nelectrolyte\n(C–C/C–H\nat 284.30 eV)"),
        ("figS15C", "C", "Al electrode,\nafter CV\nto 5 V\n(C–C/C–H\nat 283.88 eV)")]
LAB = {"C-C_C-H": "C–C/C–H", "C-O_candidate": "C–O", "O-C=O_candidate": "O–C=O", "carbonate_overlap": "carbonate (overlap)",
       "CF3_candidate": "CF$_3$", "C_oxygenated_1": "oxygenated C", "C_oxygenated_2": "oxygenated C",
       "fluorinated_C_candidate": "fluorinated C", "O_main_unresolved": "main (unresolved)", "O_high_BE": "high BE",
       "O_low_BE_unassigned": "low BE (unassigned)", "LiF_candidate": "LiF", "covalent_C-F": "covalent C–F",
       "covalent_F_candidate": "covalent C–F", "Al_III_O_candidate": "Al(III)", "Al_III_candidate": "Al(III)",
       "Al_low_BE_unassigned": "low BE (unassigned)"}
XL = {"C1s": (296, 280), "O1s": (538, 524), "F1s": (695, 680), "Al2p": (81, 67)}
TK = {"C1s": [294, 290, 286, 282], "O1s": [536, 532, 528], "F1s": [692, 688, 684], "Al2p": [78, 74, 70]}
COMP = "#d9822b"; OBS = "#7f7f7f"; ENV = "#262626"; BG = "#b0b0b0"

fig = plt.figure(figsize=(W_FULL, 13.6 * CM))
x0, wcol, y0, hrow = 0.155, 0.195, 0.115, 0.255
used = {}
for ri, (sheet, L, desc) in enumerate(ROWS):
    d = pd.read_csv(data("experimental", sheet + ".csv"))
    ybase = y0 + (2 - ri) * (hrow + 0.03)
    letter(fig, 0.005, ybase + hrow + 0.005, L)
    fig.text(0.035, ybase + hrow / 2, desc, fontsize=6.5, va="center", ha="left", linespacing=1.15)
    for ci, (reg, title) in enumerate(REG):
        ax = fig.add_axes([x0 + ci * (wcol + 0.016), ybase, wcol, hrow])
        r = d[d.region == reg]
        tr = {t: g.sort_values("Binding_energy_eV") for t, g in r.groupby("trace")}
        bg = tr["Shirley_background_CPS"]
        ax.plot(bg.Binding_energy_eV, bg.intensity_CPS, color=BG, lw=0.5)
        ob = tr["Observed_CPS"]
        ax.plot(ob.Binding_energy_eV, ob.intensity_CPS, ls="none", marker="o", ms=1.0, color=OBS)
        comps = [t for t in tr if t.endswith("_plus_background_CPS")]
        ymax = ob.intensity_CPS.max(); ymin = bg.intensity_CPS.min()
        for k, t in enumerate(comps):
            g = tr[t]; name = t[:-len("_plus_background_CPS")]
            ax.plot(g.Binding_energy_eV, g.intensity_CPS, color=COMP, lw=0.6)
            bgi = np.interp(g.Binding_energy_eV, bg.Binding_energy_eV, bg.intensity_CPS)
            j = int(np.argmax(g.intensity_CPS.values - bgi)); xe = g.Binding_energy_eV.values[j]
            lab = LAB[name]; used[name] = lab
            ax.axvline(xe, color="#9a9a9a", lw=0.4, ls=(0, (2, 2)), zorder=0)
            ax.text(xe, 1.0, lab, rotation=90, fontsize=4.8, ha="center", va="top", transform=ax.get_xaxis_transform(),
                    bbox=dict(fc="white", ec="none", pad=0.3, alpha=0.8))
        ev = tr["Envelope_CPS"]
        ax.plot(ev.Binding_energy_eV, ev.intensity_CPS, color=ENV, lw=0.7)
        ax.set_xlim(*XL[reg]); ax.set_xticks(TK[reg]); ax.set_ylim(ymin - 0.03 * (ymax - ymin), ymax + 0.75 * (ymax - ymin))
        ax.set_yticks([])
        if ri == 0: ax.set_title(title, fontsize=8, fontweight="bold", pad=3)
        if ri == 2: ax.set_xlabel("Binding energy (eV)", fontsize=6.5)
        if ci == 0: ax.set_ylabel("Intensity (a.u.)", fontsize=6.5)
        ax.tick_params(labelsize=6)
from matplotlib.lines import Line2D
h = [Line2D([], [], ls="none", marker="o", ms=2, color=OBS), Line2D([], [], color=ENV, lw=0.8),
     Line2D([], [], color=COMP, lw=0.8), Line2D([], [], color=BG, lw=0.8)]
fig.legend(h, ["measured", "fit envelope", "fitted components", "Shirley background"], loc="lower center", ncol=4,
           bbox_to_anchor=(0.58, 0.0), fontsize=6.5, frameon=False)
save(fig, "figS15")
