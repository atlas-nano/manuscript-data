"""Reweights the bulk Li+ desolvation OPES runs (collective variable cn) into the per-group free-energy profiles F(cn) and the desolvation free energies dG1 and dG2, written to data/simulation/desolv_pmf/<solvent>.npz and <solvent>.json."""
import argparse
import glob
import json
import os
import pathlib
import warnings

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
KT = 0.0252488
EQ_FRAC = 0.30
GRID = np.arange(0.0, 8.0001, 0.05)


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


def hist_w(x, w, grid=GRID):
    idx = np.clip(np.round((x - grid[0]) / (grid[1] - grid[0])).astype(int), 0, len(grid) - 1)
    return np.bincount(idx, weights=w, minlength=len(grid))


def profile(cn, wts, cn_min):
    acc = hist_w(cn, wts)
    with np.errstate(divide="ignore"):
        F = -np.log(acc / acc.sum())
    F[~np.isfinite(F)] = np.nan
    b = (GRID >= cn_min - 0.5) & (GRID <= cn_min + 0.5)
    return F - np.nanmin(F[b]), acc


def load_groups(root, cn_min):
    groups = {}
    for gd in sorted(glob.glob(os.path.join(root, "g*")), key=lambda p: int(os.path.basename(p)[1:])):
        recs = []
        for cf in sorted(glob.glob(os.path.join(gd, "COLVAR.*"))):
            names, cv = read_colvar(cf)
            if not len(cv):
                continue
            i = {n: j for j, n in enumerate(names)}
            cv = cv[int(EQ_FRAC * len(cv)):]
            recs.append(dict(cn=cv[:, i["cn"]], nc=cv[:, i["ncount"]],
                             lw=(cv[:, i["opes.bias"]] - cv[:, i["opes.rct"]]) / KT))
        if not recs:
            continue
        lwall = np.concatenate([r["lw"] for r in recs])
        wts = np.exp(lwall - lwall.max())
        F, acc = profile(np.concatenate([r["cn"] for r in recs]), wts, cn_min)
        ncr = np.round(np.concatenate([r["nc"] for r in recs])).astype(int)
        pn = {int(n): float(wts[ncr == n].sum() / wts.sum()) for n in np.unique(ncr)}
        neff_n = {}
        for n in np.unique(ncr):
            m = ncr == n
            neff_n[int(n)] = float(wts[m].sum() ** 2 / (wts[m] ** 2).sum()) if wts[m].sum() else 0.0
        groups[os.path.basename(gd)] = dict(F=F, cov=acc > 0, pn=pn, neff_n=neff_n, recs=recs)
    if not groups:
        raise SystemExit(f"no COLVAR data under {root}")
    return groups


def stats(v):
    v = np.asarray(v, float)
    return dict(mean=float(np.nanmean(v)), median=float(np.nanmedian(v)), sd=float(np.nanstd(v, ddof=1)),
                sem=float(np.nanstd(v, ddof=1) / len(v) ** 0.5))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("solvent", help="label of the output files, e.g. fm")
    ap.add_argument("run_dir", help="directory holding the bias groups g1, g2, ... each with COLVAR.0-3")
    ap.add_argument("--cn-min", type=float, required=True, help="bulk coordination-number minimum")
    ap.add_argument("--delta1", type=float, required=True, help="one contact-equivalent in cn units")
    ap.add_argument("--cn-min-alt", type=float, nargs="*", default=[], help="alternative minima for the cn_min sensitivity")
    ap.add_argument("--out", default=str(ROOT / "data" / "simulation" / "desolv_pmf"))
    a = ap.parse_args()
    warnings.simplefilter("ignore", RuntimeWarning)
    cn_min, D = a.cn_min, a.delta1
    groups = load_groups(a.run_dir, cn_min)
    gl = list(groups)
    Fs = np.array([groups[g]["F"] for g in gl])
    cov_all = np.all([groups[g]["cov"] for g in gl], axis=0)
    Fm, Fmed, Fsd = np.nanmean(Fs, 0), np.nanmedian(Fs, 0), np.nanstd(Fs, 0, ddof=1)

    def interp(F, x):
        ok = np.isfinite(F) & cov_all
        return float(np.interp(x, GRID[ok], F[ok])) if ok.sum() > 2 else float("nan")

    dg1 = {g: interp(groups[g]["F"], cn_min - D) for g in gl}
    dg2 = {g: interp(groups[g]["F"], cn_min - 2 * D) for g in gl}
    v1, v2 = np.array(list(dg1.values())), np.array(list(dg2.values()))

    def quarter(recs_by_q, q, pooled_max):
        cq, wq = [], []
        for r in recs_by_q:
            n = len(r["cn"]) // 4
            cq.append(r["cn"][q * n:(q + 1) * n])
            wq.append(r["lw"][q * n:(q + 1) * n])
        cq = np.concatenate(cq)
        if pooled_max:
            w = np.concatenate([np.exp(l - r["lw"].max()) for l, r in zip(wq, recs_by_q)])
        else:
            lq = np.concatenate(wq)
            w = np.exp(lq - lq.max())
        Fq, _ = profile(cq, w, cn_min)
        return interp(Fq, cn_min - D), interp(Fq, cn_min - 2 * D)

    per_group_q = {}
    for g in gl:
        qs = []
        for q in range(4):
            qs.append(quarter(groups[g]["recs"], q, False))
        q1, q2 = [x[0] for x in qs], [x[1] for x in qs]
        per_group_q[g] = dict(dG1_quarters=q1, dG1_spread=float(np.nanmax(q1) - np.nanmin(q1)),
                              dG2_quarters=q2, dG2_spread=float(np.nanmax(q2) - np.nanmin(q2)))
    allrecs = [r for g in gl for r in groups[g]["recs"]]
    pq = []
    for q in range(4):
        pq.append(quarter(allrecs, q, True))
    pq1, pq2 = [x[0] for x in pq], [x[1] for x in pq]

    alt = {}
    for c in sorted(set(a.cn_min_alt) | ({cn_min} if a.cn_min_alt else set())):
        vals = []
        for g in gl:
            F = groups[g]["F"]
            b = (GRID >= c - 0.5) & (GRID <= c + 0.5)
            Fa = F - np.nanmin(F[b])
            ok = np.isfinite(Fa) & cov_all
            vals.append(float(np.interp(c - D, GRID[ok], Fa[ok])) if ok.sum() > 2 else float("nan"))
        alt[f"{c:.2f}"] = dict(per_group=dict(zip(gl, vals)), **{k: v for k, v in stats(vals).items() if k != "sem"})

    os.makedirs(a.out, exist_ok=True)
    np.savez(os.path.join(a.out, f"{a.solvent}.npz"), grid=GRID, F_groups=Fs, F_mean=Fm, F_median=Fmed, F_sd=Fsd,
             coverage=cov_all, dG1=v1, dG2=v2, cn_min=cn_min, Delta1=D,
             pn=np.array([[groups[g]["pn"].get(n, 0.0) for n in range(0, 12)] for g in gl]),
             neff_n=np.array([[groups[g]["neff_n"].get(n, 0.0) for n in range(0, 12)] for g in gl]))
    res = dict(solvent=a.solvent, cn_min=cn_min, Delta1=D, groups=gl,
               dG1=dict(stats(v1), per_group=dg1), dG2=dict(stats(v2), per_group=dg2),
               quarters_per_group=per_group_q,
               quarters_pooled=dict(dG1_quarters=pq1, dG1_spread=float(np.nanmax(pq1) - np.nanmin(pq1)),
                                    dG2_quarters=pq2, dG2_spread=float(np.nanmax(pq2) - np.nanmin(pq2))),
               cn_min_sensitivity=alt)
    json.dump(res, open(os.path.join(a.out, f"{a.solvent}.json"), "w"), indent=1)
    s1, s2 = res["dG1"], res["dG2"]
    print(f"{a.solvent}: {len(gl)} groups; cn_min {cn_min:.2f}, Delta1 {D:.4f}; dG1 = {s1['mean']:.3f} kT "
          f"(median {s1['median']:.3f}, sd {s1['sd']:.3f}); dG2 = {s2['mean']:.3f} kT (sd {s2['sd']:.3f}) -> {a.out}")


if __name__ == "__main__":
    main()
