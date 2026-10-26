"""Writes the reweighted Li-anion radial potentials of mean force data/simulation/opes_pmf/<salt>_<solvent>.npz from the OPES COLVAR files of one salt-solvent run."""
import argparse
import glob
import os
import pathlib
import warnings

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
KB = 8.617333262e-5
T = 293.0
KT = KB * T
EQ_FRAC = 0.30
EDGES = np.arange(1.5, 12.55, 0.10)
CEN = 0.5 * (EDGES[1:] + EDGES[:-1])
PLATEAU = (CEN >= 10.0) & (CEN <= 12.0)
DMAX = 12.6
C0 = 1.0 / 1660.5389


def read_colvar(path):
    fields = None
    d, b = [], []
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("#! FIELDS"):
                fields = ln.split()[2:]
                continue
            if ln.startswith("#") or not ln.strip():
                continue
            p = ln.split()
            if fields is None or len(p) != len(fields):
                continue
            try:
                v = [float(x) for x in p]
            except ValueError:
                continue
            r = dict(zip(fields, v))
            if "opes.bias" not in r:
                return None
            d.append(r["d"])
            b.append(r["opes.bias"] - r.get("opes.rct", 0.0))
    if not d:
        return None
    return np.asarray(d), np.asarray(b) / KT


def hist_pmf(d, lw):
    lw = lw - lw.max()
    w = np.exp(lw)
    idx = np.digitize(d, EDGES) - 1
    ok = (idx >= 0) & (idx < len(EDGES) - 1)
    P = np.bincount(idx[ok], weights=w[ok], minlength=len(EDGES) - 1).astype(float)
    n_eff = w.sum() ** 2 / np.sum(w ** 2)
    with np.errstate(divide="ignore", invalid="ignore"):
        P = P / P.sum()
        F = -np.log(np.where(P > 0, P, np.nan))
        W = F + 2.0 * np.log(CEN)
    if np.isfinite(W[PLATEAU]).sum() >= 5:
        W = W - np.nanmean(W[PLATEAU])
    else:
        W = np.full_like(W, np.nan)
    return W, n_eff


def group_frames(gdir, dmax):
    D, LW, WID = [], [], []
    for i, f in enumerate(sorted(glob.glob(os.path.join(gdir, "COLVAR.*")))):
        got = read_colvar(f)
        if got is None:
            continue
        d, lw = got
        n0 = int(EQ_FRAC * len(d))
        D.append(d[n0:])
        LW.append(lw[n0:])
        WID.append(np.full(len(d) - n0, i, dtype=int))
    if not D:
        return None
    d = np.concatenate(D)
    lw = np.concatenate(LW)
    wid = np.concatenate(WID)
    keep = d <= dmax
    return d[keep], lw[keep], len(d), int((~keep).sum()), wid[keep]


def dG_assoc(cen, W, rb=None):
    m = np.isfinite(W)
    r, w = cen[m], W[m]
    if rb is None:
        sel = (r >= 3.0) & (r <= 8.0)
        rb = float(r[sel][np.nanargmax(w[sel])]) if sel.any() else 5.5
    b = r <= rb + 1e-6
    if b.sum() < 3:
        return rb, np.nan, np.nan
    Ka = 4 * np.pi * np.trapezoid(r[b] ** 2 * np.exp(-w[b]), r[b])
    return rb, Ka, -np.log(max(Ka * C0, 1e-300))


def profile(run_dir, groups=None, dmax=DMAX):
    gd = sorted(glob.glob(os.path.join(run_dir, "g*")), key=lambda s: int(os.path.basename(s)[1:]))
    if groups:
        gd = [g for g in gd if os.path.basename(g) in groups]
    Ws, neffs, kept, dgs, cuts, Ds, LWs, WIDs = [], [], [], [], [], [], [], []
    for g in gd:
        got = group_frames(g, dmax)
        if got is None:
            continue
        d, lw, n_raw, n_cut, wid = got
        W, ne = hist_pmf(d, lw)
        if not np.isfinite(W).any():
            continue
        Ws.append(W)
        neffs.append(ne)
        kept.append(os.path.basename(g))
        cuts.append(n_cut)
        Ds.append(d)
        LWs.append(lw)
        WIDs.append(wid)
        dgs.append(dG_assoc(CEN, W)[2])
    Ws = np.vstack(Ws)
    with np.errstate(invalid="ignore"):
        mean = np.nanmean(Ws, axis=0)
        sd = np.nanstd(Ws, axis=0, ddof=1) if len(Ws) > 1 else np.zeros_like(mean)
    cov = np.isfinite(Ws).sum(axis=0)
    rb, ka, dg = dG_assoc(CEN, mean)
    qdg = []
    for q in range(4):
        Wq = []
        for d, lw, wid in zip(Ds, LWs, WIDs):
            sel = []
            for w in np.unique(wid):
                idx = np.flatnonzero(wid == w)
                m = len(idx)
                sel.append(idx[q * m // 4:(q + 1) * m // 4])
            sel = np.concatenate(sel)
            Wg, _ = hist_pmf(d[sel], lw[sel])
            if np.isfinite(Wg).any():
                Wq.append(Wg)
        if Wq:
            with np.errstate(invalid="ignore"):
                mq = np.nanmean(np.vstack(Wq), axis=0)
            qdg.append(dG_assoc(CEN, mq, rb)[2])
    return dict(cen=CEN, mean=mean, sd=sd, cov=cov, Ws=Ws, kept=np.array(kept), neff=np.array(neffs),
                dgs=np.array(dgs), rb=rb, dg=dg, ka=ka, cuts=np.array(cuts), qdg=np.array(qdg))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("system", help="salt_solvent label of the output file, e.g. lipf6_acn")
    ap.add_argument("run_dir", help="directory holding the bias groups g1, g2, ... each with COLVAR.0-3")
    ap.add_argument("--groups", nargs="*", help="bias groups to use (default: every group with a usable profile)")
    ap.add_argument("--dmax", type=float, default=DMAX, help="frames with d above this (A) are dropped")
    ap.add_argument("--out", default=str(ROOT / "data" / "simulation" / "opes_pmf"))
    a = ap.parse_args()
    warnings.simplefilter("ignore", RuntimeWarning)
    res = profile(a.run_dir, a.groups, a.dmax)
    os.makedirs(a.out, exist_ok=True)
    out = os.path.join(a.out, f"{a.system}.npz")
    np.savez(out, **res)
    i = int(np.nanargmin(np.where((CEN >= 1.5) & (CEN <= 3.5), res["mean"], np.inf)))
    print(f"{a.system}: {len(res['kept'])} groups ({', '.join(res['kept'])}); frames dropped at d > {a.dmax} A: "
          f"{int(res['cuts'].sum())}; contact W {res['mean'][i]:+.3f} kT at {CEN[i]:.2f} A; "
          f"dG_assoc(rb {res['rb']:.2f} A) {res['dg']:+.3f} kT -> {out}")


if __name__ == "__main__":
    main()
