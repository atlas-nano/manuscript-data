"""Writes the ion-pair association table data/tables/opes_assoc.csv (contact minimum, barrier, association free energy at rb = 4.25 A and the LiTPFA - LiPF6 difference) from the OPES COLVAR files of each salt-solvent run."""
import argparse
import csv
import glob
import os
import pathlib
import warnings

import numpy as np

from analyse_opes_pmf import C0, KT, read_colvar

ROOT = pathlib.Path(__file__).resolve().parents[1]
LO, BW, NB = 1.5, 0.10, 110
CEN = LO + BW * (np.arange(NB) + 0.5)
DCUT = 12.6
EPS = 1e-6
RB = 4.25
DISCARD = 0.30
ZERO = (10.0, 12.0)
ALT_ZERO = [(9.0, 11.0), (9.0, 12.0), (10.0, 11.0), (11.0, 12.0), (9.5, 12.0), (11.5, 12.0)]
SALTS = ("litpfa", "lipf6", "libf4")
SOLVENTS = ("fm", "acn", "me2o", "dme", "dee")
SALT_NAME = {"litpfa": "LiTPFA", "lipf6": "LiPF6", "libf4": "LiBF4"}
SOLV_NAME = {"fm": "FM", "acn": "ACN", "me2o": "Me2O", "dme": "DME", "dee": "DEE"}


def load_run(run_dir, groups=None):
    out = {}
    gdirs = sorted(glob.glob(os.path.join(run_dir, "g*")), key=lambda s: int(os.path.basename(s)[1:]))
    for g in gdirs:
        name = os.path.basename(g)
        if groups and name not in groups:
            continue
        ws = []
        for f in sorted(glob.glob(os.path.join(g, "COLVAR.*")), key=lambda s: int(s.rsplit(".", 1)[1])):
            got = read_colvar(f)
            if got is not None:
                ws.append(got)
        out[name] = ws
    return out


def frames(ws, part=None):
    D, L = [], []
    for d, lw in ws:
        n = len(d)
        idx = np.arange(int(DISCARD * n), n)
        if part is not None:
            k, K = part
            m = len(idx)
            idx = idx[k * m // K:(k + 1) * m // K]
        idx = idx[d[idx] <= DCUT]
        D.append(d[idx])
        L.append(lw[idx])
    return np.concatenate(D), np.concatenate(L)


def pmf(d, lw, zero=ZERO):
    wt = np.exp(lw - lw.max())
    i = np.floor((d - LO) / BW + 1e-12).astype(int)
    ok = (i >= 0) & (i < NB)
    P = np.bincount(i[ok], weights=wt[ok], minlength=NB)
    with np.errstate(divide="ignore", invalid="ignore"):
        W = np.where(P > 0, -np.log(P), np.nan) + 2.0 * np.log(CEN)
    zw = (CEN >= zero[0] - EPS) & (CEN <= zero[1] + EPS) & np.isfinite(W)
    if zw.sum() < 5:
        return np.full(NB, np.nan)
    return W - W[zw].mean()


def dG(W, rb=RB):
    m = np.isfinite(W) & (CEN <= rb + EPS)
    if m.sum() < 3:
        return np.nan
    r = CEN[m]
    return -np.log(4 * np.pi * np.trapezoid(r ** 2 * np.exp(-W[m]), r) * C0)


def profiles(run, part=None, zero=ZERO):
    Ws = []
    for ws in run.values():
        W = pmf(*frames(ws, part), zero=zero)
        if np.isfinite(W).any():
            Ws.append(W)
    return np.vstack(Ws)


def pooled(Ws):
    with np.errstate(invalid="ignore"):
        return np.nanmean(Ws, axis=0)


def band(W, lo, hi):
    return (CEN >= lo - EPS) & (CEN <= hi + EPS) & np.isfinite(W)


def analyse(run, salt):
    Ws = profiles(run)
    mean = pooled(Ws)
    base = dG(mean)
    per = np.array([dG(W) for W in Ws])
    gsd = float(np.std(per[np.isfinite(per)], ddof=1))
    q = [dG(pooled(profiles(run, part=(k, 4)))) for k in range(4)]
    zmeans = {z: pooled(profiles(run, zero=z)) for z in ALT_ZERO}
    zterm = float(np.sqrt(np.mean([(dG(zmeans[z]) - base) ** 2 for z in ALT_ZERO])))
    lo = 1.5 if salt == "libf4" else 2.0
    m = band(mean, lo, 3.5)
    j = np.flatnonzero(m)[np.argmin(mean[m])]
    contact_edge = j == np.flatnonzero(band(np.zeros(NB), lo, 3.5))[0] or j == np.flatnonzero(m)[-1]
    b = band(mean, 3.0, 8.0)
    k = np.flatnonzero(b)[np.argmax(mean[b])]
    edge = k in (np.flatnonzero(b)[0], np.flatnonzero(b)[-1])
    return dict(n=len(Ws), mean=mean, base=base, gsd=gsd, qsd=float(np.std(q, ddof=1)), zterm=zterm, zmeans=zmeans,
                cmin=float(mean[j]), cmin_r=float(CEN[j]), cmin_sd=float(np.nanstd(Ws[:, j], ddof=1)),
                contact_edge=bool(contact_edge), cwin=(float(CEN[m][0]), float(CEN[m][-1]), float(np.min(mean[m]))),
                bar=None if edge else (float(mean[k]), float(np.nanstd(Ws[:, k], ddof=1)), float(CEN[k])))


def fmt(x, s=None, nd=1):
    t = f"{x:+.{nd}f}" if s is None else f"{x:+.{nd}f} ± {s:.{nd}f}"
    return t.replace("-", "−")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run", action="append", required=True, metavar="SYSTEM=DIR[:g1,g2,...]",
                    help="salt_solvent label, its run directory of bias groups and optionally the groups to use")
    ap.add_argument("--out", default=str(ROOT / "data" / "tables" / "opes_assoc.csv"))
    a = ap.parse_args()
    warnings.simplefilter("ignore", RuntimeWarning)
    res = {}
    for spec in a.run:
        name, rest = spec.split("=", 1)
        path, _, gl = rest.partition(":")
        run = load_run(path, gl.split(",") if gl else None)
        res[name] = analyse(run, name.split("_")[0])
    rows = []
    for solv in SOLVENTS:
        for salt in SALTS:
            key = f"{salt}_{solv}"
            if key not in res:
                continue
            r = res[key]
            unc = r["gsd"] if salt == "libf4" else float(np.sqrt(r["gsd"] ** 2 + r["qsd"] ** 2 + r["zterm"] ** 2))
            if r["contact_edge"] and r["cmin"] > 0:
                contact = f"none sampled; W ≥ {fmt(r['cwin'][2], nd=0)} over {r['cwin'][0]:.2f}–{r['cwin'][1]:.2f} Å"
            else:
                word = "minimum" if r["cmin"] < 0 else "local minimum"
                contact = f"{word} {fmt(r['cmin'], r['cmin_sd'])} at {r['cmin_r']:.2f} Å"
            barrier = "none" if r["bar"] is None else f"{fmt(r['bar'][0], r['bar'][1])} at {r['bar'][2]:.2f} Å"
            ddg, ddg_v, ddg_u = "", "", ""
            if salt == "lipf6" and f"litpfa_{solv}" in res:
                t = res[f"litpfa_{solv}"]
                dd = t["base"] - r["base"]
                gsd = float(np.hypot(t["gsd"], r["gsd"]))
                qsd = float(np.hypot(t["qsd"], r["qsd"]))
                zt = float(np.sqrt(np.mean([(dG(t["zmeans"][z]) - dG(r["zmeans"][z]) - dd) ** 2 for z in ALT_ZERO])))
                tot = float(np.sqrt(gsd ** 2 + qsd ** 2 + zt ** 2))
                ddg, ddg_v, ddg_u = fmt(dd, tot), f"{dd:.4f}", f"{tot:.4f}"
            rows.append(dict(salt=SALT_NAME[salt], solvent=SOLV_NAME[solv], contact_region=contact, barrier=barrier,
                             dG_assoc=fmt(r["base"], unc), ddG_assoc=ddg, n_groups=r["n"], dG_assoc_kT=f"{r['base']:.4f}",
                             dG_group_sd=f"{r['gsd']:.4f}", dG_quarter_sd=f"{r['qsd']:.4f}", dG_zero_window=f"{r['zterm']:.4f}",
                             dG_uncertainty=f"{unc:.4f}", ddG_assoc_kT=ddg_v, ddG_uncertainty=ddg_u))
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0]))
        wr.writeheader()
        wr.writerows(rows)
    for row in rows:
        print(" | ".join(str(row[k]) for k in ("salt", "solvent", "contact_region", "barrier", "dG_assoc", "ddG_assoc")))
    print(f"-> {a.out}")


if __name__ == "__main__":
    main()
