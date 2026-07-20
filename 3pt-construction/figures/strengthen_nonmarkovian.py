#!/usr/bin/env python3
"""STRENGTHEN: cage efficiency is a UNIVERSAL function of kernel non-Markovianity γ/Ω₀.
Tests:
 (1) Collapse: LJ-64fs, LJ-8fs, and the 7 metals fall on ONE eff(γ/Ω₀) curve. The
     64fs aliasing just displaces LJ to higher γ/Ω₀ (fabricated memory) ALONG the curve.
 (2) Non-circular: γ/Ω₀ uses only the raw Volterra kernel K̃(0) and the DoS 2nd moment
     Ω₀ — NOT the cage clip/gate/filter or (Wg-Ws). Efficiency uses the cage.
 (3) Gap constancy: the target deficit (S_ref-S_rigHS) is ~flat, so efficiency is driven
     by the cage (numerator=memory), not the denominator.
 (4) A 2nd intrinsic memory metric (raw relative memory area ∫|F_K-F_M|/∫F_K) also predicts
     efficiency — robustness that it's memory, not a γ/Ω₀ coincidence.
"""
import numpy as np, csv, glob, re, sys, matplotlib
from pathlib import Path
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import czt
sys.path.insert(0, str(Path(__file__).resolve().parent))
from kernel_tools import invert_kernel_matrix, H, VLIGHT, KB, PI, NA
DATA = Path(__file__).resolve().parent.parent / "data"
HERE = Path(__file__).resolve().parent
c = VLIGHT*1e-10

def kernel_metrics(npz):
    """Return γ/Ω₀ (non-Markovianity) and the raw relative memory area, from the kernel+DoS
    ONLY (no clip/gate/filter)."""
    d = np.load(npz)
    nu = d['nu_cm'].astype(float); tot = d['dos_total'].astype(float)
    C = d['C_scalar'].astype(float); dt = float(d['dt']); dnu = nu[1]-nu[0]
    wa = 2*PI*nu*c
    Om0 = np.sqrt(np.trapezoid(wa**2*tot, dx=1.0)/np.trapezoid(tot, dx=1.0))/(2*PI*c)
    cn = C/C[0]; Cm = np.zeros((cn.size,3,3))
    for i in range(3): Cm[:,i,i] = cn
    K = np.einsum('tii->t', invert_kernel_matrix(dt, Cm))/3.0
    nz = np.nonzero(K)[0]
    if not nz.size: return None
    nK_full = int(nz[-1])+1
    # main-lobe safeguard (matches the production py-xPT engine / SI Sec. S2):
    # gamma = Ktilde(0) integrated over the full kernel support picks up a smooth,
    # spurious post-main-lobe tail (the nz[-1] artifact). Detect the main-lobe
    # cutoff (first |K|<0.02*K0 after the initial samples) and fall back to it when
    # the full-support friction inflates by more than 2x over the main-lobe value.
    def _czt(nK):
        Kv = K[:nK].copy(); Kv[0]*=0.5; Kv[-1]*=0.5
        return czt(Kv, m=wa.size, w=np.exp(-1j*2*PI*dnu*c*dt), a=1.0)*dt
    K0 = abs(K[0])
    _below = np.where(np.abs(K[3:]) < 0.02*K0)[0] if K0 > 0.0 else np.array([], dtype=int)
    nK_main = int(_below[0])+3 if _below.size else nK_full
    Ktil = _czt(nK_full); gamma = float(Ktil[0].real)
    if nK_main < nK_full:
        Ktil_m = _czt(nK_main); g_main = float(Ktil_m[0].real)
        if abs(g_main) > 0.0 and abs(gamma - g_main) > 1.0*abs(g_main):
            Ktil = Ktil_m; gamma = g_main
    F_K = np.real(1.0/(1j*wa+Ktil)); F_M = gamma/(gamma**2+wa**2)
    # raw relative memory area over the diffusive band (nu up to ~2*Om0)
    band = nu <= 2*Om0/(2*PI*c)
    memarea = float(np.trapezoid(np.abs(F_K-F_M)[band], dx=dnu)/max(np.trapezoid(np.abs(F_K)[band], dx=dnu),1e-30))
    g_cm = gamma/(2*PI*c)
    return dict(nonmark=g_cm/Om0, memarea=memarea, gam=g_cm, Om0=Om0)

def rt(nm):
    m=re.search(r"r(\d+)_T(\d+)",nm); r=m.group(1);T=m.group(2)
    return (f"{float(r[0]+'.'+r[1:]):.3f}", f"{float(T[0]+'.'+T[1:]):.3f}")

def load_lj(npglob, csvpath):
    cc={(r["rho_star"],r["T_star"]):r for r in csv.DictReader(open(csvpath))}
    pts=[]
    for npz in sorted(glob.glob(npglob)):
        nm=re.search(r"cin_(r\w+)_grp1",npz).group(1); k=rt(nm)
        row=cc.get(k)
        if not row or not (0.70<=float(k[0])<=1.02): continue
        km=kernel_metrics(npz)
        if km is None: continue
        targ=float(row["S_ref"])-float(row["rigHS"]); deliv=float(row["3PT"])-float(row["rigHS"])
        if targ<=0: continue
        pts.append(dict(x=km['nonmark'], mem=km['memarea'], eff=deliv/targ,
                        targ=targ, rho=float(k[0]), key=k))
    return pts

LJ8  = load_lj(str(DATA/"cage_sens"/"8fs"/"cin_r*_grp1.npz"),  str(DATA/"lj_grid_8fs.csv"))
LJ64 = load_lj(str(DATA/"cage_sens"/"64fs"/"cin_r*_grp1.npz"), str(DATA/"lj_grid_64fs.csv"))

METALS=["Al","Ag","Ni","Pd","Rh","Au","Ir"]
TI=dict(zip(METALS,[9.25,10.24,10.31,11.08,11.10,11.71,12.25]))
RIG=dict(zip(METALS,[9.089,9.897,10.023,10.725,10.757,11.365,11.848]))
S3=dict(zip(METALS,[9.314,10.148,10.294,11.051,11.078,11.674,12.160]))
MET=[]
for m in METALS:
    km=kernel_metrics(str(DATA/"cage_sens"/"metals"/f"cin_{m}_grp1.npz"))
    if km is None: continue
    targ=TI[m]-RIG[m]; deliv=S3[m]-RIG[m]
    MET.append(dict(x=km['nonmark'], mem=km['memarea'], eff=deliv/targ, targ=targ, name=m))

# ---- (3) gap constancy ----
allt=[p['targ'] for p in LJ8+MET]; allx=[p['x'] for p in LJ8+MET]
print(f"(3) target gap (ref-rigHS): mean={np.mean(allt):.3f} std={np.std(allt):.3f} "
      f"range[{min(allt):.2f},{max(allt):.2f}]  corr(gap,γ/Ω₀)={np.corrcoef(allt,allx)[0,1]:+.2f}")

# ---- (1)+(4) universal curve fit on pooled LJ8+LJ64+metals ----
POOL = LJ8+LJ64+MET
X=np.array([p['x'] for p in POOL]); E=np.array([p['eff'] for p in POOL]); M=np.array([p['mem'] for p in POOL])
def sat(x,emax,k): return emax*(1-np.exp(-k*x))
p0,_=curve_fit(sat, X, E, p0=[0.95,0.8], maxfev=10000)
Efit=sat(X,*p0); ss=1-np.sum((E-Efit)**2)/np.sum((E-E.mean())**2)
print(f"(1) universal fit eff={p0[0]:.2f}*(1-exp(-{p0[1]:.2f}·γ/Ω₀)) ; pooled R²={ss:.3f} (N={len(POOL)})")
print(f"    corr(eff, γ/Ω₀) pooled = {np.corrcoef(X,E)[0,1]:+.3f}")
print(f"(4) corr(eff, raw memory area) pooled = {np.corrcoef(M,E)[0,1]:+.3f} ; corr(γ/Ω₀, memarea)={np.corrcoef(X,M)[0,1]:+.3f}")

# ---- plot ----
fig,ax=plt.subplots(1,2,figsize=(13,5.3))
a=ax[0]
xs=np.linspace(0.6,5.2,200); a.plot(xs,sat(xs,*p0),'k-',lw=1.5,label=f"universal fit (R²={ss:.2f})")
a.scatter([p['x'] for p in LJ64],[p['eff'] for p in LJ64],c='tab:orange',marker='s',s=40,alpha=.8,edgecolor='k',lw=.3,label="LJ 64fs (aliased)")
a.scatter([p['x'] for p in LJ8],[p['eff'] for p in LJ8],c='tab:blue',marker='o',s=40,alpha=.8,edgecolor='k',lw=.3,label="LJ 8fs (resolved)")
a.scatter([p['x'] for p in MET],[p['eff'] for p in MET],c='crimson',marker='D',s=55,edgecolor='k',lw=.3,label="metals (TI)")
for p in MET: a.annotate(p['name'],(p['x'],p['eff']),fontsize=6.5,xytext=(3,2),textcoords='offset points')
# per-state 64fs->8fs shift arrows (aliasing slides each state ALONG the curve)
d8={p['key']:p for p in LJ8}
ndown=0
for p in LJ64:
    q=d8.get(p['key'])
    if not q: continue
    a.annotate("", xy=(q['x'],q['eff']), xytext=(p['x'],p['eff']),
               arrowprops=dict(arrowstyle="->", color="0.4", lw=0.6, alpha=0.6))
    if q['x']<p['x']: ndown+=1
print(f"    64fs->8fs shift arrows: {ndown}/{sum(1 for p in LJ64 if p['key'] in d8)} point to LOWER γ/Ω₀ (down the curve)")
a.axhline(1,color='0.5',ls=':',lw=1)
a.set_xlabel(r"non-Markovianity  $\gamma/\Omega_0$"); a.set_ylabel(r"cage efficiency $\Delta S_{cage}/(S_{ref}-S_{rigHS})$")
a.set_title("(1) universal collapse: 64fs aliasing slides LJ UP the same curve")
a.legend(fontsize=8); a.grid(alpha=.3)
# arrows linking 64fs->8fs for matched states
d64={ (round(p['x'],99)):p for p in LJ64}
b=ax[1]
b.scatter(M,E,c=['tab:blue']*len(LJ8)+['tab:orange']*len(LJ64)+['crimson']*len(MET),s=35,alpha=.7,edgecolor='k',lw=.3)
b.set_xlabel(r"raw relative memory area  $\int|F_K-F_M|/\int F_K$"); b.set_ylabel("cage efficiency")
b.set_title(f"(4) independent memory metric also predicts eff (r={np.corrcoef(M,E)[0,1]:.2f})")
b.grid(alpha=.3)
fig.suptitle("Cage efficiency is a universal function of kernel non-Markovianity (non-circular: γ/Ω₀ & memory-area use only the raw kernel)",fontsize=11)
fig.tight_layout(rect=[0,0,1,0.95]); fig.savefig(HERE/"strengthen_universal_nonmarkovian.png",dpi=140)
print("wrote strengthen_universal_nonmarkovian.png")
