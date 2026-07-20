#!/usr/bin/env python3
"""Add WATER to the universal cage-efficiency vs non-Markovianity curve.
Water is placed two ways:
 (A) per-CHANNEL, reference-FREE: the translational cage uses the measured universal
     translational deficit (0.333 k_B/DoF from gap-constancy on LJ+metals), so
     eff_trans = ΔS_cage^trans / 0.333 — needs NO water free-energy reference.
 (B) molecular TOTAL, with a reference band: eff = ΔS_cage^tot/(S_ref - S_2PTrig),
     S_ref swept over the firm-estimate range (FEP/FL-TI ~59-61.5 for TIP4P/2005,
     ~64-66 for SPC/Ew) → vertical error bar; horizontal = cage-weighted γ/Ω₀.
"""
import numpy as np, csv, glob, re, sys, matplotlib
from pathlib import Path
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
sys.path.insert(0, str(Path(__file__).resolve().parent))
from strengthen_nonmarkovian import kernel_metrics, load_lj, LJ8, LJ64, MET, sat
from cage_decompose_compare import diag
HERE = Path(__file__).resolve().parent
WDIR = HERE.parent / "data" / "cage_sens" / "water"

DEFICIT = 0.333  # k_B per translational DoF (gap-constancy, LJ+metals)

# universal fit from the pooled LJ+metals set
POOL = LJ8+LJ64+MET
X=np.array([p['x'] for p in POOL]); E=np.array([p['eff'] for p in POOL])
p0,_=curve_fit(sat, X, E, p0=[0.95,0.8], maxfev=10000)

# ---- water channels ----
WAT = {"TIP4P/2005": (str(WDIR/"cin_tip4p_grp1.npz"), str(WDIR/"cin_tip4p_grp1_rot.npz"), 54.578, (59.0,60.5,61.5)),
       "SPC/Ew":     (str(WDIR/"cin_spcew_grp1.npz"), str(WDIR/"cin_spcew_grp1_rot.npz"), 57.055, (64.0,65.33,66.0))}
R=8.314
print("=== Water on the universal curve ===")
wat_trans=[]; wat_tot=[]
for name,(tn,rn,base,refband) in WAT.items():
    dt=diag(tn); dr=diag(rn); kt=kernel_metrics(tn); kr=kernel_metrics(rn)
    # (A) trans channel reference-free
    eff_trans = dt['dS']/DEFICIT
    pred = sat(kt['nonmark'], *p0)
    wat_trans.append((kt['nonmark'], eff_trans, name))
    print(f" {name} TRANS: γ/Ω₀={kt['nonmark']:.2f}  cage={dt['dS']:.3f}  eff=cage/0.333={eff_trans:.2f}  curve-pred={pred:.2f}")
    print(f" {name} ROT  : γ/Ω₀={kr['nonmark']:.2f}  cage={dr['dS']:.3f}  (rot deficit≈{dr['dS']/sat(kr['nonmark'],*p0):.2f} k_B; saturated)")
    # (B) molecular total, ref band
    cage_tot=(dt['dS']+dr['dS'])  # k_B/mol (per-DoF prefactor already in each)
    gw=(dt['dS']*kt['nonmark']+dr['dS']*kr['nonmark'])/cage_tot  # cage-weighted γ/Ω₀
    effs=[cage_tot/((ref-base)/R) for ref in refband]
    wat_tot.append((gw, effs[1], min(effs),max(effs), name))
    print(f" {name} TOTAL: cage-wtd γ/Ω₀={gw:.2f}  cage_tot={cage_tot:.3f}k_B  eff={effs[1]:.2f} (band {min(effs):.2f}-{max(effs):.2f} over ref {refband[0]}-{refband[2]})")

# ---- plot ----
fig,ax=plt.subplots(figsize=(8,5.6))
xs=np.linspace(0.6,6.5,200); ax.plot(xs,sat(xs,*p0),'k-',lw=1.6,label=f"universal fit  0.97(1−e^(−0.70γ/Ω₀))")
ax.scatter([p['x'] for p in LJ64],[p['eff'] for p in LJ64],c='tab:orange',marker='s',s=28,alpha=.55,edgecolor='none',label="LJ 64fs (aliased)")
ax.scatter([p['x'] for p in LJ8],[p['eff'] for p in LJ8],c='tab:blue',marker='o',s=30,alpha=.7,edgecolor='k',lw=.2,label="LJ 8fs")
ax.scatter([p['x'] for p in MET],[p['eff'] for p in MET],c='crimson',marker='D',s=55,edgecolor='k',lw=.3,label="metals (TI)")
# water trans (reference-free)
for x,e,nm in wat_trans:
    ax.scatter(x,e,c='teal',marker='*',s=260,edgecolor='k',lw=.6,zorder=5)
    ax.annotate(nm+" trans",(x,e),fontsize=7.5,xytext=(6,-2),textcoords='offset points')
ax.scatter([],[],c='teal',marker='*',s=160,edgecolor='k',label="water trans (ref-free, deficit=0.333)")
# water total (ref band)
for gw,e,elo,ehi,nm in wat_tot:
    ax.errorbar(gw,e,yerr=[[e-elo],[ehi-e]],fmt='P',color='purple',ms=10,capsize=3,zorder=5,mec='k',mew=.5)
    ax.annotate(nm+" total",(gw,e),fontsize=7.5,xytext=(6,4),textcoords='offset points')
ax.errorbar([],[],fmt='P',color='purple',label="water total (FEP ref ± band)")
ax.axhline(1,color='0.6',ls=':',lw=1)
ax.set_xlabel(r"non-Markovianity  $\gamma/\Omega_0$"); ax.set_ylabel(r"cage efficiency  $\Delta S_{cage}/$deficit")
ax.set_title("Universal cage-efficiency law with WATER:\ntrans channels land on the curve reference-free; molecular totals sit at the saturated end")
ax.legend(fontsize=8,loc='lower right'); ax.grid(alpha=.3); ax.set_ylim(0.2,1.6)
fig.tight_layout(); fig.savefig(HERE/"fig_universal_nonmarkovian.png",dpi=140)
fig.savefig(HERE/"fig_universal_nonmarkovian.pdf")
print("wrote fig_universal_nonmarkovian.png/.pdf")
