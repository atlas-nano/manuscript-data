"""Reweights the two-dimensional OPES run at the Li(100) interface (collective variables h and cn) into the free-energy surfaces opes2d_surfaces.npz, the region levels and the contact cost in opes2d_results.json from the COLVAR, KERNELS, STATE and PLUMED files of its consecutive segment directories."""
import argparse
import glob
import json
import math
import os
import pathlib
import re
import warnings

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
KT = 0.0252488
EQ_FRAC = 0.30
DELTA1 = 0.7179
GRID = np.arange(0.0, 8.0001, 0.05)
GRID2 = np.arange(0.0, 8.0001, 0.1)
HF = np.round(np.arange(-1.0, 16.0 + 1e-9, 0.1), 10)
H2 = np.round(np.arange(-1.0, 16.0 + 1e-9, 0.2), 10)
HFC = 0.5 * (HF[1:] + HF[:-1])
H2C = 0.5 * (H2[1:] + H2[:-1])
C_HI, L1, RR, H_QHI = 3.0, (3.0, 6.0), (11.0, 14.0), 14.5
NAN = float("nan")
KEYS = ("t", "h", "cn", "cn12", "nc", "bias", "rct", "lw", "dli0")


def read_colvar(path):
    names, rows = None, []
    for l in open(path, errors="ignore"):
        if l.startswith("#! FIELDS"):
            names = l.split()[2:]
            continue
        if l.startswith("#"):
            continue
        s = l.split()
        if names and len(s) == len(names):
            try:
                rows.append([float(x) for x in s])
            except ValueError:
                pass
    return names, np.array(rows)


def cidx(x, grid=GRID):
    return np.clip(np.round((x - grid[0]) / (grid[1] - grid[0])).astype(int), 0, len(grid) - 1)


def hist_w(x, w, grid=GRID):
    return np.bincount(cidx(x, grid), weights=w, minlength=len(grid))


def hidx(h, edges):
    dh = edges[1] - edges[0]
    i = np.floor((h - edges[0]) / dh + 1e-9).astype(int)
    return np.clip(i, 0, len(edges) - 2)


def nlog(x):
    with np.errstate(divide="ignore", invalid="ignore"):
        y = -np.log(x)
    y = np.asarray(y, float)
    y[~np.isfinite(y)] = np.nan
    return y


def kish(w):
    s2 = (w ** 2).sum()
    return float(w.sum() ** 2 / s2) if s2 > 0 else 0.0


def nwalkers_per_bin(idx, wk, nbins):
    if len(idx) == 0:
        return np.zeros(nbins, int)
    key = np.unique(wk.astype(np.int64) * nbins + idx)
    return np.bincount(key % nbins, minlength=nbins)


def tau_int(x):
    d = np.asarray(x, float) - np.mean(x)
    tau = 1.0
    for k in range(1, max(2, len(d) // 4)):
        rho = (d[:-k] * d[k:]).mean() / d.var() if d.var() > 0 else 0
        if rho <= 0:
            break
        tau += 2 * rho
    return tau


def read_barrier(seg_dirs):
    files = sorted(set(f for d in seg_dirs for f in glob.glob(os.path.join(d, "plumed.w*.s*.dat"))))
    if not files:
        raise SystemExit("no plumed.w*.s*.dat file in the segment directories")
    vals = set()
    for f in files:
        m = re.search(r"BARRIER=([0-9.eE+-]+)", open(f, errors="ignore").read())
        vals.add(float(m.group(1)) if m else None)
    if len(vals) != 1 or None in vals:
        raise SystemExit("the PLUMED files do not carry one common OPES BARRIER")
    return vals.pop()


def load_segments(seg_dirs, nwalk, nseg, dt):
    seg_len = nseg * dt
    W, names_ref = [], None
    for k in range(1, nwalk + 1):
        parts = []
        for s, d in enumerate(seg_dirs, 1):
            fn = os.path.join(d, "COLVAR.%d" % (k - 1))
            names, A = read_colvar(fn)
            if names is None or not len(A):
                raise SystemExit(f"{fn} has no readable rows")
            names_ref = names_ref or names
            if names != names_ref:
                raise SystemExit(f"FIELDS differ in {fn}")
            ti = names.index("time")
            t0 = (s - 1) * seg_len
            t = A[:, ti]
            if s > 1:
                if t.min() >= t0 - 1e-6:
                    pass
                elif abs(t.min()) < 1e-6 and t.max() <= seg_len + 1e-6:
                    A = A.copy()
                    A[:, ti] += t0
                elif t.max() > t0 - 1e-6:
                    A = A[A[:, ti] >= t0 - 1e-6]
                else:
                    raise SystemExit(f"{fn} times {t.min():g}..{t.max():g} do not fit segment {s}")
            parts.append((s, A))
        A = np.concatenate([p[1] for p in parts])
        segid = np.concatenate([np.full(len(p[1]), p[0]) for p in parts])
        t = A[:, names_ref.index("time")]
        _, idx = np.unique(t[::-1], return_index=True)
        keep = np.sort(len(t) - 1 - idx)
        A = A[keep]
        segid = segid[keep]
        o = np.argsort(A[:, names_ref.index("time")], kind="stable")
        A = A[o]
        F = {n: A[:, i] for i, n in enumerate(names_ref)}
        F["_seg"] = segid[o]
        W.append(F)
    return W


def read_plumed_table(path):
    names, sets, rows = None, {}, []
    for l in open(path, errors="ignore"):
        if l.startswith("#! FIELDS"):
            names = l.split()[2:]
            continue
        if l.startswith("#! SET"):
            p = l.split()
            sets[p[2]] = p[3] if len(p) > 3 else ""
            continue
        if l.startswith("#"):
            continue
        s = l.split()
        if names and len(s) == len(names):
            try:
                rows.append([float(x) for x in s])
            except ValueError:
                pass
    return names, sets, np.array(rows)


class OpesBias:

    def __init__(self, seg_dirs, barrier, nseg, dt):
        rows = []
        prev_max = -np.inf
        for s, d in enumerate(seg_dirs, 1):
            names, sets, K = read_plumed_table(os.path.join(d, "KERNELS"))
            col = [names.index(n) for n in ("time", "h", "cn", "sigma_h", "sigma_cn", "height", "logweight")]
            K = K[:, col]
            t0 = (s - 1) * nseg * dt
            if s > 1 and K[:, 0].min() < t0 - 1e-6 and K[:, 0].max() <= nseg * dt + 1e-6 and abs(K[:, 0].min()) < 1e-6:
                K = K.copy()
                K[:, 0] += t0
            K = K[K[:, 0] > prev_max + 1e-9]
            if len(K):
                prev_max = K[:, 0].max()
            rows.append(K)
        self.K = np.concatenate(rows)
        names, sets, S = read_plumed_table(os.path.join(seg_dirs[-1], "STATE"))
        col = [names.index(n) for n in ("time", "h", "cn", "sigma_h", "sigma_cn", "height")]
        self.S = S[:, col]
        self.gamma = float(sets["biasfactor"])
        self.eps = float(sets["epsilon"])
        self.cut = float(sets["kernel_cutoff"])
        self.thr2 = float(sets["compression_threshold"]) ** 2
        self.zed_state = float(sets["zed"])
        self.sw_state = float(sets["sum_weights"])
        self.kbt = barrier / self.gamma
        self.vac = math.exp(-0.5 * self.cut ** 2)

    def kde(self, c, s, h, X, chunk=256):
        out = np.zeros(len(X))
        c2 = self.cut ** 2
        for i0 in range(0, len(X), chunk):
            x = X[i0:i0 + chunk]
            d2 = (((x[:, None, :] - c[None]) / s[None]) ** 2).sum(2)
            v = np.where(d2 < c2, np.exp(-0.5 * d2) - self.vac, 0.0)
            out[i0:i0 + chunk] = v @ h
        return out

    def zed(self, c, s, h, norm):
        return self.kde(c, s, h, c).sum() / norm / len(h)

    def replay(self, cuts):
        cuts = sorted(cuts)
        K = self.K
        n, cap = 0, len(K) + 1
        C = np.empty((cap, 2))
        Sg = np.empty((cap, 2))
        H = np.empty(cap)
        out, ci = {}, 0
        norm = 0.0

        def merge(a, b_h, b_c, b_s):
            hh = H[a] + b_h
            cc = (H[a] * C[a] + b_h * b_c) / hh
            ss = (H[a] * (Sg[a] ** 2 + C[a] ** 2) + b_h * (b_s ** 2 + b_c ** 2)) / hh - cc ** 2
            H[a] = hh
            C[a] = cc
            Sg[a] = np.sqrt(ss)

        for row in K:
            while ci < len(cuts) and row[0] > cuts[ci] + 1e-9:
                out[cuts[ci]] = (C[:n].copy(), Sg[:n].copy(), H[:n].copy(), norm)
                ci += 1
            c = row[1:3]
            s = row[3:5]
            h = row[5]
            norm += math.exp(row[6])
            if n:
                d2 = (((c - C[:n]) / Sg[:n]) ** 2).sum(1)
                k = int(np.argmin(d2))
                if d2[k] < self.thr2:
                    merge(k, h, c.copy(), s.copy())
                    giver = k
                    while True:
                        d2 = (((C[giver] - C[:n]) / Sg[:n]) ** 2).sum(1)
                        d2[giver] = np.inf
                        t = int(np.argmin(d2))
                        if d2[t] >= self.thr2:
                            break
                        taker, g = (t, giver) if t < giver else (giver, t)
                        merge(taker, H[g], C[g].copy(), Sg[g].copy())
                        C[g:n - 1] = C[g + 1:n].copy()
                        Sg[g:n - 1] = Sg[g + 1:n].copy()
                        H[g:n - 1] = H[g + 1:n].copy()
                        n -= 1
                        giver = taker
                    continue
            C[n] = c
            Sg[n] = s
            H[n] = h
            n += 1
        while ci < len(cuts):
            out[cuts[ci]] = (C[:n].copy(), Sg[:n].copy(), H[:n].copy(), norm)
            ci += 1
        return out

    def F2d(self, snap, zed=None):
        c, s, h, norm = snap
        z = self.zed(c, s, h, norm) if zed is None else zed
        HH, CC = np.meshgrid(HFC, GRID, indexing="ij")
        X = np.column_stack([HH.ravel(), CC.ravel()])
        return (-np.log(self.kde(c, s, h, X) / norm / z + self.eps)).reshape(HH.shape)

    @staticmethod
    def Fh_from_2d(F2):
        m = np.nanmin(F2)
        return nlog(np.exp(-(F2 - m)).sum(1) * (GRID[1] - GRID[0])) + m


def make_rec(F, w, face, sl):
    b, r = F["opes.bias"][sl], F["opes.rct"][sl]
    return dict(w=w, face=face, t=F["time"][sl], h=F["h"][sl], cn=F["cn"][sl], cn12=F["cn12"][sl], nc=F["ncount"][sl],
                bias=b, rct=r, lw=(b - r) / KT, dli0=F["dli0.lowest"][sl])


def pool(recs):
    recs = [r for r in recs if len(r["h"])]
    if not recs:
        return None
    P = {k: np.concatenate([r[k] for r in recs]) for k in KEYS}
    P["wk"] = np.concatenate([np.full(len(r["h"]), r["w"]) for r in recs])
    P["w"] = np.exp(P["lw"] - P["lw"].max())
    return P


def h_profile(P):
    return np.bincount(hidx(P["h"], HF), weights=P["w"], minlength=len(HF) - 1)


def levels_from_mass(Pf, h_lo):
    tot = Pf.sum()

    def lev(lo, hi):
        if not np.isfinite(lo) or hi - lo <= 1e-9:
            return NAN
        m = (HF[:-1] >= lo - 1e-9) & (HF[1:] <= hi + 1e-9)
        mass = Pf[m].sum()
        return float(-np.log(mass / tot / (hi - lo))) if mass > 0 else NAN

    FR, FC, FL1 = lev(*RR), lev(h_lo, C_HI), lev(*L1)
    P2 = Pf.reshape(-1, 2).sum(1)
    F2 = nlog(P2 / tot / 0.2)
    mC = (H2[:-1] >= h_lo - 1e-9) & (H2[1:] <= C_HI + 1e-9) & np.isfinite(F2)
    if mC.any():
        j = np.where(mC)[0][np.argmin(F2[mC])]
        minC, locC = float(F2[j] - FR), float(H2C[j])
    else:
        minC, locC = NAN, NAN
    return dict(FC=FC, FL1=FL1, FR=FR, dGC=FC - FR, dGL1=FL1 - FR, dGminC=minC, locC=locC,
                F1112=lev(11.0, 12.0), F1314=lev(13.0, 14.0))


def cn_cond(P, mask):
    if P is None or not mask.any():
        return None
    acc = hist_w(P["cn"][mask], P["w"][mask])
    F = nlog(acc / acc.sum())
    F = F - np.nanmin(F)
    nw = nwalkers_per_bin(cidx(P["cn"][mask]), P["wk"][mask], len(GRID))
    return dict(F=F, nw=nw, nstar=float(GRID[int(np.argmax(acc))]), acc=acc)


def dg1_read(F, nw, c, dshift, cov_min):
    b = (GRID >= c - 0.5 - 1e-9) & (GRID <= c + 0.5 + 1e-9)
    if not np.isfinite(F[b]).any():
        return NAN
    Fa = F - np.nanmin(F[b])
    ok = np.isfinite(Fa) & (nw >= cov_min)
    need = (GRID >= c - dshift - 0.05 - 1e-9) & (GRID <= c + 1e-9)
    if cov_min > 0 and not ok[need].all():
        return NAN
    return float(np.interp(c - dshift, GRID[ok], Fa[ok])) if ok.sum() > 2 else NAN


def quantities(P, h_lo, dshift):
    if P is None:
        return dict(dGC=NAN, dGL1=NAN, dGminC=NAN, locC=NAN, dG1R=NAN, nR=NAN, F1112=NAN, F1314=NAN, FR=NAN)
    q = levels_from_mass(h_profile(P), h_lo)
    cc = cn_cond(P, (P["h"] >= RR[0]) & (P["h"] <= RR[1]))
    if cc is None:
        q.update(dG1R=NAN, nR=NAN)
    else:
        q.update(nR=cc["nstar"], dG1R=dg1_read(cc["F"], cc["nw"], cc["nstar"], dshift, 0))
    return q


def spread(vals):
    v = np.array(vals, float)
    v = v[np.isfinite(v)]
    return float(v.max() - v.min()) if len(v) >= 2 else NAN


def analyse(W, faces, barrier, dshift, bias):
    res = dict(nwalk=len(W))
    t_start = min(F["time"][0] for F in W)
    t_end = min(F["time"][-1] for F in W)
    recs = []
    for k, F in enumerate(W):
        n0 = int(EQ_FRAC * len(F["time"]))
        recs.append(make_rec(F, k + 1, faces[k], slice(n0, None)))
    t_a = max(r["t"][0] for r in recs)
    P = pool(recs)
    Pf = h_profile(P)
    res.update(t_start=t_start, t_end=t_end, t_analysed_from=t_a, rows_analysed=len(P["h"]), kish_neff=kish(P["w"]))

    i2 = hidx(P["h"], H2)
    nw_h2 = nwalkers_per_bin(i2, P["wk"], len(H2) - 1)
    occ = np.where(nw_h2 >= 4)[0]
    h_lo = float(H2[occ[0]]) if len(occ) else NAN
    mR = (P["h"] >= RR[0]) & (P["h"] <= RR[1])
    ccR = cn_cond(P, mR)
    nR = ccR["nstar"] if ccR else NAN
    cap = nR + 1.5
    res.update(h_lo=h_lo, nstar_R=nR, cn_cap=cap)

    tot = P["w"].sum()
    ic2 = cidx(P["cn"], GRID2)
    flat = i2 * len(GRID2) + ic2
    nb2 = (len(H2) - 1) * len(GRID2)
    P2 = np.bincount(flat, weights=P["w"], minlength=nb2).reshape(len(H2) - 1, len(GRID2))
    P2sq = np.bincount(flat, weights=P["w"] ** 2, minlength=nb2).reshape(P2.shape)
    nw2 = nwalkers_per_bin(flat, P["wk"], nb2).reshape(P2.shape)
    with np.errstate(invalid="ignore", divide="ignore"):
        neff2 = np.where(P2sq > 0, P2 ** 2 / P2sq, 0.0)
    F2 = nlog(P2 / tot / (0.2 * 0.1))
    F2 -= np.nanmin(F2)
    Fh01 = nlog(Pf / tot / 0.1)
    Fh02 = nlog(Pf.reshape(-1, 2).sum(1) / tot / 0.2)
    med2 = np.full(nb2, np.nan)
    o = np.argsort(flat, kind="stable")
    fs = flat[o]
    bs = P["bias"][o]
    u, st = np.unique(fs, return_index=True)
    for a, b0, b1 in zip(u, st, list(st[1:]) + [len(fs)]):
        med2[a] = np.median(bs[b0:b1])
    med2 = med2.reshape(P2.shape)
    arrays = dict(h_edges_02=H2, cn_centres_01=GRID2, F_h_cn=F2, nwalk_h_cn=nw2, neff_h_cn=neff2,
                  median_bias_eV_h_cn=med2, h_edges_01=HF, F_h_01=Fh01, F_h_02=Fh02, cn_grid_005=GRID)
    hin = ((H2[:-1] >= h_lo - 1e-9) & (H2[1:] <= H_QHI + 1e-9))[:, None]
    cin = (GRID2 <= cap + 1e-9)[None, :]
    arrays["reliable_mask_h_cn"] = (nw2 >= 4) & (med2 >= -barrier + 2 * KT) & hin & cin

    q = levels_from_mass(Pf, h_lo)
    res["levels"] = q
    sel = P["cn"] <= cap + 1e-9
    res["levels_cn_capped"] = levels_from_mass(h_profile(dict(P, h=P["h"][sel], w=P["w"][sel])), h_lo)
    regs = {"C": (P["h"] >= h_lo - 1e-9) & (P["h"] <= C_HI) if np.isfinite(h_lo) else np.zeros(len(P["h"]), bool),
            "L1": (P["h"] > L1[0]) & (P["h"] <= L1[1]), "R": mR, "all": np.ones(len(P["h"]), bool)}
    Fcn = {k: cn_cond(P, m) for k, m in regs.items()}
    arrays.update({f"F_cn_{k}": (v["F"] if v else np.full(len(GRID), np.nan)) for k, v in Fcn.items()})
    arrays.update({f"nwalk_cn_{k}": (v["nw"] if v else np.zeros(len(GRID), int)) for k, v in Fcn.items()})
    res["dG1R"] = dg1_read(ccR["F"], ccR["nw"], nR, dshift, 8) if ccR else NAN
    res["dG1R_alt"] = {f"{c:.2f}": dg1_read(ccR["F"], ccR["nw"], c, dshift, 8) for c in (nR - 0.05, nR, nR + 0.05)} if ccR else {}
    res["nstar_C"] = Fcn["C"]["nstar"] if Fcn["C"] else NAN

    path = []
    if np.isfinite(h_lo):
        e = np.arange(math.floor(h_lo / 0.5) * 0.5, H_QHI + 1e-9, 0.5)
        e = np.unique(np.concatenate([[h_lo], e[e > h_lo + 1e-9]]))
        PR_n = (ccR["acc"][int(round(nR / 0.05))] / tot) if ccR else NAN
        FpathR = float(-np.log(PR_n / ((RR[1] - RR[0]) * 0.05)) - q["FR"]) if ccR and PR_n > 0 else NAN
        res["F_path_R"] = FpathR
        for lo, hi in zip(e[:-1], e[1:]):
            m = (P["h"] >= lo - 1e-9) & (P["h"] < hi - 1e-9)
            cc = cn_cond(P, m)
            if cc is None:
                path.append(dict(lo=float(lo), hi=float(hi), n=0))
                continue
            j = int(round(cc["nstar"] / 0.05))
            Fp = float(-np.log(cc["acc"][j] / tot / ((hi - lo) * 0.05)) - q["FR"])
            ww = P["w"][m]
            path.append(dict(lo=float(lo), hi=float(hi), n=int(m.sum()), nw=int(len(np.unique(P["wk"][m]))),
                             nstar=cc["nstar"], dF_path=Fp - FpathR,
                             cn_mean=float(np.average(P["cn"][m], weights=ww)),
                             nc_mean=float(np.average(P["nc"][m], weights=ww)),
                             cn12_mean=float(np.average(P["cn12"][m], weights=ww)),
                             dnstar=cc["nstar"] - nR, dG1=dg1_read(cc["F"], cc["nw"], cc["nstar"], dshift, 8), neff=kish(ww)))
    res["path"] = path
    res["per_walker"] = [dict(w=r["w"], face=r["face"], **{k: quantities(pool([r]), h_lo, dshift)[k] for k in ("dGC", "dGL1")})
                         for r in recs]

    tcs = list(np.arange(t_start + 50.0, t_end + 1e-6, 50.0))
    if not tcs or abs(tcs[-1] - t_end) > 1e-6:
        tcs.append(t_end)
    bias_ts = {}
    if bias is not None:
        snaps = bias.replay(tcs)
        for tc in tcs:
            bias_ts[tc] = OpesBias.Fh_from_2d(bias.F2d(snaps[tc]))
    ts_rows = []
    for tc in tcs:
        rr = []
        for k, F in enumerate(W):
            m = np.where(F["time"] <= tc + 1e-9)[0]
            n0 = int(EQ_FRAC * len(m))
            rr.append(make_rec(F, k + 1, faces[k], m[n0:]))
        qq = quantities(pool(rr), h_lo, dshift)
        row = dict(t=float(tc), dGC=qq["dGC"], dGL1=qq["dGL1"], dGminC=qq["dGminC"])
        if tc in bias_ts:
            Pb = np.exp(-(bias_ts[tc] - np.nanmin(bias_ts[tc])))
            Pb[~np.isfinite(Pb)] = 0.0
            qb = levels_from_mass(Pb * 0.1, h_lo)
            row.update(bias_dGC=qb["dGC"], bias_dGL1=qb["dGL1"])
        ts_rows.append(row)
    res["timeseries"] = ts_rows

    quarts = []
    for qi in range(4):
        rr = []
        for r in recs:
            n = len(r["h"]) // 4
            sl = slice(qi * n, (qi + 1) * n)
            rr.append({k: (r[k][sl] if k in KEYS else r[k]) for k in r})
        quarts.append(quantities(pool(rr), h_lo, dshift))
    halves = {nm: quantities(pool([r for r in recs if r["w"] % 2 == par]), h_lo, dshift) for nm, par in (("odd", 1), ("even", 0))}
    blocks = []
    tb = t_a
    while tb + 50.0 <= t_end + 1e-6:
        rr = [{k: (r[k][(r["t"] >= tb - 1e-9) & (r["t"] < tb + 50 - 1e-9)] if k in KEYS else r[k]) for k in r} for r in recs]
        blocks.append((tb, quantities(pool(rr), h_lo, dshift)))
        tb += 50.0

    def block_se(key):
        v = np.array([b[1][key] for b in blocks], float)
        v = v[np.isfinite(v)]
        if len(v) < 3:
            return NAN, len(v)
        return float(v.std(ddof=1) / math.sqrt(len(v)) * math.sqrt(tau_int(v))), len(v)

    U = {}
    for key in ("dGC", "dGL1", "dG1R", "dGminC"):
        qs = [x[key] for x in quarts]
        hd = abs(halves["odd"][key] - halves["even"][key])
        se, nbk = block_se(key)
        fin = [c for c in (hd / 2, spread(qs) / 2, se) if np.isfinite(c)]
        U[key] = dict(quarters=qs, quarter_spread=spread(qs), odd=halves["odd"][key], even=halves["even"][key], half_diff=hd,
                      block_se=se, n_blocks=nbk, u=(max(fin) if fin else NAN))
    res["u"] = U
    jk = {}
    for key in ("dGC", "dGminC", "dGL1"):
        th = np.array([quantities(pool([r for r in recs if r["w"] != wi]), h_lo, dshift)[key]
                       for wi in sorted(set(r["w"] for r in recs))], float)
        fin = th[np.isfinite(th)]
        n_ = len(fin)
        jk[key] = dict(se=(float(math.sqrt((n_ - 1) / n_ * ((fin - fin.mean()) ** 2).sum())) if n_ > 1 else NAN), n=n_)
    res["jackknife"] = jk

    if bias is not None:
        snap = (bias.S[:, 1:3], bias.S[:, 3:5], bias.S[:, 5], bias.sw_state)
        F2b = bias.F2d(snap, zed=bias.zed_state)
        Fhb = OpesBias.Fh_from_2d(F2b)
        Pb = np.exp(-(Fhb - np.nanmin(Fhb)))
        Pb[~np.isfinite(Pb)] = 0.0
        Fb02 = nlog(Pb.reshape(-1, 2).sum(1) / Pb.sum() / 0.2)
        res["levels_bias"] = levels_from_mass(Pb * 0.1, h_lo)
        arrays.update(F_bias_h_01=nlog(Pb / Pb.sum() / 0.1), F_bias_h_02=Fb02, F_bias_h_cn_fine=(F2b - np.nanmin(F2b)))

    if set(faces) == {"top", "bottom"}:
        res["faces"] = {}
        for fc in ("bottom", "top"):
            Pp = pool([r for r in recs if r["face"] == fc])
            Pq = h_profile(Pp)
            arrays[f"F_h_02_{fc}"] = nlog(Pq.reshape(-1, 2).sum(1) / Pq.sum() / 0.2)
            qf = quantities(Pp, h_lo, dshift)
            res["faces"][fc] = dict(n_walkers=sum(f == fc for f in faces), dGC=qf["dGC"], dGL1=qf["dGL1"])
    return res, arrays


def clean(o):
    if isinstance(o, dict):
        return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, (np.floating, float)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


def fmt(x, n=2):
    return "nan" if x is None or not np.isfinite(x) else f"{x:.{n}f}"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seg-dirs", nargs="+", required=True, help="segment directories in time order, each with "
                    "COLVAR.0-15, KERNELS, STATE and plumed.w*.s*.dat")
    ap.add_argument("--walker-map", help="json with a 'walkers' list giving each walker's 'face' (top or bottom)")
    ap.add_argument("--out", default=str(ROOT / "data" / "simulation" / "opes2d"))
    ap.add_argument("--nwalk", type=int, default=16)
    ap.add_argument("--nseg", type=int, default=1000000, help="MD steps per segment")
    ap.add_argument("--dt", type=float, default=0.0005, help="MD time step (ps)")
    ap.add_argument("--delta1", type=float, default=DELTA1, help="bulk first-shell coordination step (cn units)")
    ap.add_argument("--no-bias", action="store_true", help="skip the bias-derived surfaces")
    a = ap.parse_args()
    warnings.simplefilter("ignore", RuntimeWarning)
    seg_dirs = [os.path.abspath(d) for d in a.seg_dirs]
    barrier = read_barrier(seg_dirs)
    W = load_segments(seg_dirs, a.nwalk, a.nseg, a.dt)
    faces = ["?"] * a.nwalk
    if a.walker_map:
        faces = [w["face"] for w in json.load(open(a.walker_map))["walkers"]][:a.nwalk]
    bias = None if a.no_bias else OpesBias(seg_dirs, barrier, a.nseg, a.dt)
    res, arrays = analyse(W, faces, barrier, a.delta1, bias)
    res.update(barrier_eV=barrier, kT_eV=KT, delta1=a.delta1)
    os.makedirs(a.out, exist_ok=True)
    np.savez(os.path.join(a.out, "opes2d_surfaces.npz"), **{k: np.asarray(v) for k, v in arrays.items()})
    json.dump(clean(res), open(os.path.join(a.out, "opes2d_results.json"), "w"), indent=1)
    q, U = res["levels"], res["u"]
    print(f"h_lo {fmt(res['h_lo'], 1)} A; n*_R {fmt(res['nstar_R'])}; n*_C {fmt(res['nstar_C'])}")
    print(f"dG_C = {fmt(q['dGC'], 3)} kT = {fmt(q['dGC'] * KT, 4)} eV +/- {fmt(U['dGC']['u'])} kT "
          f"(jackknife SE {fmt(res['jackknife']['dGC']['se'])}, block SE {fmt(U['dGC']['block_se'])})")
    print(f"dG_min,C = {fmt(q['dGminC'], 3)} kT at h = {fmt(q['locC'])} A +/- {fmt(U['dGminC']['u'])} kT")
    print(f"dG_L1 = {fmt(q['dGL1'], 3)} kT +/- {fmt(U['dGL1']['u'])} kT; dG1,R = {fmt(res['dG1R'], 3)} kT")
    if "levels_bias" in res:
        b = res["levels_bias"]
        print(f"bias-derived: dG_C {fmt(b['dGC'], 3)}  dG_L1 {fmt(b['dGL1'], 3)}  dG_min,C {fmt(b['dGminC'], 3)} kT")
    print(f"-> {a.out}/opes2d_surfaces.npz, opes2d_results.json")


if __name__ == "__main__":
    main()
