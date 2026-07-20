#!/usr/bin/env python
"""3PT path-independence / state-function validation (ms2a).

Shows the point-by-point 3PT entropies lie on a single integrable surface:
(1) ΔS_3PT(A→B) = S_3PT(B)−S_3PT(A) matches the thermodynamic-integration ΔS from
    the MBWR EOS over all liquid pairs;
(2) the MBWR TI is path-independent (two L-paths agree to ~1e-9) — the reference is
    a consistent surface and ΔS_TI is well defined;
(3) the 23 independent 3PT values are fit by one smooth S(rho*,T*) surface as tightly
    as the MBWR EOS itself (integrability proxy).

Figure: (a) ΔS_3PT vs ΔS_TI scatter over all liquid pairs; (b) residual histogram.
"""
import csv, math, itertools
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys; sys.path.insert(0, str(Path(__file__).resolve().parent))
import mbwr_lj as mb

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"

# --- liquid 3PT values from the grid ---
liq = []
for r in csv.DictReader(open(DATA / "lj_fullgrid_cage.csv")):
    rho, T = float(r["rho_star"]), float(r["T_star"])
    if 0.70 <= rho <= 1.02:
        liq.append((rho, T, float(r["3PT"]), float(r["S_ref"])))   # ref = MBWR for liquids

# --- TI integrand from MBWR EOS (central differences) ---
def dSdT_rho(rho, T, h=1e-4): return (mb.s_total_reduced(rho, T+h) - mb.s_total_reduced(rho, T-h)) / (2*h)
def dSdrho_T(rho, T, h=1e-4): return (mb.s_total_reduced(rho+h, T) - mb.s_total_reduced(rho-h, T)) / (2*h)

def simpson(fn, a, b, n=400):
    n += n % 2
    x = np.linspace(a, b, n+1); y = np.array([fn(xi) for xi in x])
    return (b-a)/n/3 * (y[0] + y[-1] + 4*y[1:-1:2].sum() + 2*y[2:-1:2].sum())

def ti_two_paths(A, B):
    (rA, TA), (rB, TB) = A, B
    p1 = (simpson(lambda r: dSdrho_T(r, TA), rA, rB) +   # isotherm@TA then isochore@rB
          simpson(lambda T: dSdT_rho(rB, T),  TA, TB))
    p2 = (simpson(lambda T: dSdT_rho(rA, T),  TA, TB) +   # isochore@rA then isotherm@TB
          simpson(lambda r: dSdrho_T(r, TB), rA, rB))
    return p1, p2

# --- (1)+(2): all liquid pairs, ΔS_3PT vs ΔS_TI (=MBWR(B)-MBWR(A), the path integral) ---
dS3, dSti = [], []
for (rA,TA,t3A,_), (rB,TB,t3B,_) in itertools.combinations(liq, 2):
    dS3.append(t3B - t3A)
    dSti.append(mb.s_total_reduced(rB,TB) - mb.s_total_reduced(rA,TA))
dS3, dSti = np.array(dS3), np.array(dSti)
resid = dS3 - dSti
rms = math.sqrt((resid**2).mean())
print(f"(1) {len(dS3)} liquid pairs: ΔS_3PT−ΔS_TI  mean={resid.mean():+.4f}  RMS={rms:.4f}  max|{np.abs(resid).max():.3f}|  k_B")
_ptrms = math.sqrt(np.mean([(t-r)**2 for _,_,t,r in liq]))
print(f"    point |3PT−MBWR| RMS = {_ptrms:.4f};  sqrt2*that = {math.sqrt(2)*_ptrms:.4f}")

# representative distant pairs: two-path TI agreement
print("(2) two-path TI agreement (path1 vs path2; ΔS_3PT for reference):")
S3 = {(round(rho,3),round(T,3)): t3 for rho,T,t3,_ in liq}
for A,B in [((0.70,1.4),(1.00,1.4)), ((0.85,0.9),(0.85,1.4)), ((0.70,1.0),(0.85,1.4))]:
    p1,p2 = ti_two_paths(A,B)
    d3 = S3[(round(B[0],3),round(B[1],3))] - S3[(round(A[0],3),round(A[1],3))]
    print(f"    {A}->{B}: TI p1={p1:+.4f} p2={p2:+.4f} (Δpath={p1-p2:+.1e})  ΔS_3PT={d3:+.4f}  3PT−TI={d3-p1:+.4f}")

# --- (3): smooth-surface integrability proxy (quadratic fit residual) ---
def fit_resid(vals):
    rho = np.array([x[0] for x in liq]); T = np.array([x[1] for x in liq]); y = np.array(vals)
    M = np.column_stack([np.ones_like(rho), rho, T, rho**2, rho*T, T**2])
    c,_,_,_ = np.linalg.lstsq(M, y, rcond=None)
    return math.sqrt(np.mean((y - M@c)**2))
print(f"(3) quadratic-surface fit residual: 3PT={fit_resid([x[2] for x in liq]):.4f}  MBWR={fit_resid([x[3] for x in liq]):.4f} k_B")

# --- figure ---
plt.rcParams.update({"font.size": 10, "axes.linewidth": 0.8, "figure.dpi": 150, "savefig.bbox": "tight"})
fig, (axA, axB) = plt.subplots(1, 2, figsize=(8.0, 3.6))
lim = [min(dSti.min(),dS3.min())-0.1, max(dSti.max(),dS3.max())+0.1]
axA.plot(lim, lim, "-", color="0.5", lw=1.0, zorder=1)
axA.plot(dSti, dS3, "o", ms=3, color="#55a868", alpha=0.5, zorder=2)
axA.set_xlim(lim); axA.set_ylim(lim); axA.set_aspect("equal")
axA.set_xlabel(r"$\Delta S_{\rm TI}$ (MBWR path integral, $k_B$/atom)")
axA.set_ylabel(r"$\Delta S_{\rm 3PT}=S^*_{\rm 3PT}(B)-S^*_{\rm 3PT}(A)$")
axA.set_title(f"(a) all {len(dS3)} liquid pairs", fontsize=9.5)
axA.grid(ls=":", lw=0.5, alpha=0.6)
axB.hist(resid, bins=30, color="#55a868", alpha=0.85)
axB.axvline(0, color="k", lw=0.9)
axB.set_xlabel(r"$\Delta S_{\rm 3PT}-\Delta S_{\rm TI}$  ($k_B$/atom)")
axB.set_ylabel("pairs")
axB.set_title(f"(b) residual: mean {resid.mean():+.3f}, RMS {rms:.3f}", fontsize=9.5)
axB.grid(axis="y", ls=":", lw=0.5, alpha=0.6)
fig.tight_layout()
for ext in ("png","pdf"): fig.savefig(HERE / f"fig_path_independence.{ext}")
print("wrote fig_path_independence.png/.pdf")
