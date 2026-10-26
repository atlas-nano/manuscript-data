"""Assembles the ion-pairing thermodynamic cycle of Figures 3F and S3 into data/simulation/pairing_cycle.json."""
import csv, json, pathlib
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
SIM = ROOT / "data" / "simulation"
TABLES = ROOT / "data" / "tables"
COLS = ("c_ggla", "dH", "dS_Li", "dS_an", "dG_solv", "dG_mix", "band_lo", "band_hi")
C0 = 1.0 / 1660.5389


def per_box(path):
    out = {}
    with open(path, encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            out.setdefault(r["system"], []).append({c: float(r[c]) for c in COLS})
    return out


def diss_tds():
    with open(TABLES / "crystal_dissociation.csv", encoding="utf-8") as fh:
        diss = {r["salt"]: float(r["dG_diss"]) for r in csv.DictReader(fh)}
    with open(TABLES / "ideal_mixing.csv", encoding="utf-8") as fh:
        tds = {r["solvent"]: float(r["TdS_ideal"]) for r in csv.DictReader(fh)}
    return diss, tds


def dG_assoc(cen, W, rb=None):
    m = np.isfinite(W)
    r, w = cen[m], W[m]
    if rb is None:
        sel = (r >= 3.0) & (r <= 8.0)
        rb = float(r[sel][np.nanargmax(w[sel])]) if sel.any() else 5.5
    b = r <= rb + 1e-6
    if b.sum() < 3:
        return rb, np.nan, np.nan
    trap = getattr(np, "trapezoid", None) or np.trapz
    Ka = 4 * np.pi * trap(r[b] ** 2 * np.exp(-w[b]), r[b])
    return rb, Ka, -np.log(max(Ka * C0, 1e-300))


SALTS = ("LiTPFA", "LiPF6", "LiBF4"); SOLVS = ("FM", "ACN")
KB = 1.98720425e-3
KT_PMF = KB * 293.0
C_BOX = 1660.5389 / 38.0 ** 3
GGLA_CUT = -10.0
NBOOT, SEED = 20000, 20261005
NPZ = {("LiPF6", "FM"): "lipf6_fm.npz", ("LiBF4", "FM"): "libf4_fm.npz",
       ("LiTPFA", "ACN"): "litpfa_acn.npz", ("LiPF6", "ACN"): "lipf6_acn.npz",
       ("LiBF4", "ACN"): "libf4_acn.npz"}
FOOTING_MISMATCH = {("LiPF6", "FM"), ("LiTPFA", "FM")}


def assoc(salt, solv):
    if (salt, solv) == ("LiTPFA", "FM"):
        z = np.load(SIM / "opes_pmf" / "litpfa_fm.npz")
        _, _, dg = dG_assoc(z["cen"], z["pm_opes_analytic_mean"])
        return float(dg), float("nan")
    z = np.load(SIM / "opes_pmf" / NPZ[(salt, solv)])
    return float(z["dg"]), float(np.std(z["dgs"], ddof=1))


def g_pair(dg_kt, sd_kt):
    x = np.exp(-dg_kt) * C_BOX
    p = x / (1 + x)
    g = -KT_PMF * np.log1p(x)
    sig = p * (sd_kt if np.isfinite(sd_kt) else 0.0) * KT_PMF
    return g, sig, p


def ols(x, y):
    A = np.vstack([np.ones_like(x), x]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    return coef


def main():
    rows = per_box(SIM / "xpt_dgmix_anion_clamp_293K.csv")
    diss, tds = diss_tds()
    rng = np.random.default_rng(SEED)
    out = {}
    for solv in SOLVS:
        for salt in SALTS:
            b = rows[f"{salt.lower()}_{solv.lower()}"]
            U = np.array([r["c_ggla"] for r in b])
            G_tot0 = np.array([r["dG_mix"] for r in b])
            dga, sda = assoc(salt, solv)
            gp, sgp, p = g_pair(dga, sda)
            free = U > GGLA_CUT
            rec = dict(n=len(b), n_free=int(free.sum()), U_LA_median=float(np.median(U)),
                       U_LA_box=U.tolist(), dG_mix_box=G_tot0.tolist(), dG_assoc_kT=dga, sd_assoc_kT=sda,
                       p_paired=p, dG_pair=gp, s_pair=sgp, totals={})
            for tot, G in (("0", G_tot0), ("B", G_tot0 + U)):
                med = float(np.median(G))
                if free.sum() >= 3:
                    route = "free-box median"
                    ref_fn = lambda Gs, Us: np.median(Gs[Us > GGLA_CUT]) if (Us > GGLA_CUT).any() else np.nan
                elif (salt, solv) == ("LiPF6", "FM"):
                    route = "OLS intercept at U_LA=0"
                    ref_fn = lambda Gs, Us: ols(Us, Gs)[0] if np.ptp(Us) > 1e-9 else np.nan
                else:
                    route = None
                t = dict(med=med)
                bs_med = [np.median(G[rng.integers(0, len(G), len(G))]) for _ in range(4000)]
                t["s_med"] = float(np.std(bs_med, ddof=1))
                if route:
                    ref = float(ref_fn(G, U))
                    Rs, refs = [], []
                    for _ in range(NBOOT):
                        i = rng.integers(0, len(G), len(G))
                        r_ = ref_fn(G[i], U[i])
                        if np.isfinite(r_):
                            refs.append(r_)
                            Rs.append(np.median(G[i]) - r_ - gp)
                    R = med - ref - gp
                    t.update(route=route, ref=ref, s_ref=float(np.std(refs, ddof=1)), R=R,
                             s_R=float(np.hypot(np.std(Rs, ddof=1), sgp)),
                             dG_solv_free=ref - diss[salt] + tds[solv])
                else:
                    t.update(route="residual", dG_solv_free=med - diss[salt] + tds[solv] - gp,
                             s_ref=float(np.hypot(t["s_med"], sgp)))
                rec["totals"][tot] = t
            out[f"{salt}/{solv}"] = rec
    with open(SIM / "pairing_cycle.json", "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1, default=float)
    for k, r in out.items():
        print(f"{k:11s} median dG_mix {r['totals']['B']['med']:+8.2f} +- {r['totals']['B']['s_med']:.2f} kcal/mol")


if __name__ == "__main__":
    main()
