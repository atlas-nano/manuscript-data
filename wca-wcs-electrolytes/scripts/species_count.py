"""Counts the ionic species (connected Li-anion contact components) in every production frame of the fluoromethane salt trajectories and writes data/simulation/species_count/species_count.json with one ions-only snapshot per salt."""
import argparse
import json
import os
import pathlib

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
CIP_CUT = 3.0
DROP = 0.5
CUT_SWEEP = (2.5, 3.0, 3.5, 4.0, 4.5)
TYPES = {"lipf6": dict(anion=range(6, 8), li=8), "litpfa": dict(anion=range(6, 11), li=11)}
ELEMENTS = {"litpfa": {1: "F", 2: "C", 3: "H", 4: "H", 5: "H", 6: "C", 7: "C", 8: "O", 9: "F", 10: "Al", 11: "Li"},
            "lipf6": {1: "F", 2: "C", 3: "H", 4: "H", 5: "H", 6: "P", 7: "F", 8: "Li"}}


def frames(path):
    f = open(path)
    while True:
        line = f.readline()
        if not line:
            return
        if not line.startswith("ITEM: TIMESTEP"):
            continue
        ts = int(f.readline())
        f.readline()
        n = int(f.readline())
        f.readline()
        L = []
        for _ in range(3):
            lo, hi = (float(x) for x in f.readline().split()[:2])
            L.append(hi - lo)
        hdr = f.readline().split()[2:]
        col = {k: i for i, k in enumerate(hdr)}
        ids = np.empty(n, int)
        mol = np.empty(n, int)
        typ = np.empty(n, int)
        r = np.empty((n, 3))
        for k in range(n):
            p = f.readline().split()
            ids[k] = int(p[col["id"]])
            mol[k] = int(p[col["mol"]])
            typ[k] = int(p[col["type"]])
            r[k] = [float(p[col[c]]) for c in ("xu", "yu", "zu")]
        yield ts, (np.array(L), ids, mol, typ, r)


def dedup_frames(path):
    keep = {}
    order = []
    for ts, fr in frames(path):
        if ts not in keep:
            order.append(ts)
        keep[ts] = fr
    return [(ts, keep[ts]) for ts in sorted(order)]


def ion_blocks(typ, T):
    li = np.where(typ == T["li"])[0]
    an = np.where(np.isin(typ, list(T["anion"])))[0]
    nat = len(an) // max(1, len(li))
    if len(an) != nat * len(li):
        raise SystemExit("anion atoms do not divide evenly among the Li+")
    return li, an, np.arange(len(an)) // nat, nat


def min_dist_matrix(L, r, li, an, parent, n_ion):
    d = r[li][:, None, :] - r[an][None, :, :]
    d -= L * np.round(d / L)
    dist = np.linalg.norm(d, axis=-1)
    out = np.full((len(li), n_ion), np.inf)
    for p in range(n_ion):
        out[:, p] = dist[:, parent == p].min(axis=1)
    return out


def n_species_from(md, cut):
    near = md < cut
    n_li, n_an = near.shape
    par = list(range(n_li + n_an))

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x

    for i, j in zip(*np.where(near)):
        a, b = find(int(i)), find(n_li + int(j))
        if a != b:
            par[a] = b
    return len({find(k) for k in range(n_li + n_an)}), near.any(axis=0)


def scan(path, salt, drop):
    dd = dedup_frames(path)
    n_all = len(dd)
    keep = dd[int(drop * n_all):]
    T = TYPES[salt]
    ts, mds, frs = [], [], []
    for t, (L, ids, mol, typ, r) in keep:
        li, an, parent, nat = ion_blocks(typ, T)
        ts.append(t)
        mds.append(min_dist_matrix(L, r, li, an, parent, len(li)))
        frs.append((L, typ, r, li, an, parent))
    mds = np.array(mds)
    nsp, paired = {}, {}
    for c in CUT_SWEEP:
        rows = [n_species_from(m, c) for m in mds]
        nsp[c] = np.array([x[0] for x in rows])
        paired[c] = np.array([x[1] for x in rows])
    return dict(ts=np.array(ts), nsp=nsp, paired=paired, md=mds, frames=frs, n_all=n_all, n_kept=len(keep))


def quartiles(a):
    a = np.asarray(a, float)
    return tuple(float(np.percentile(a, p)) for p in (25, 50, 75))


def choose(per_rep):
    meds = {k: float(np.median(v["nsp"][CIP_CUT])) for k, v in per_rep.items()}
    across = float(np.median(list(meds.values())))
    rep = sorted(meds, key=lambda k: (abs(meds[k] - across), k))[0]
    s = per_rep[rep]
    nsp = s["nsp"][CIP_CUT]
    med = float(np.median(nsp))
    idx = int(np.argmin(np.abs(np.asarray(nsp, float) - med)))
    q1, q2, q3 = quartiles(nsp)
    return rep, idx, dict(median=med, q1=q1, q3=q3, iqr=q3 - q1, n_considered=int(s["n_kept"]), n_raw=int(s["n_all"]),
                          ts=int(s["ts"][idx]), nsp=int(nsp[idx]), across_replica_median=across, replica_medians=meds,
                          ts_first=int(s["ts"][0]), ts_last=int(s["ts"][-1]))


def stability(s, idx):
    cls = {c: s["paired"][c][idx] for c in CUT_SWEEP}
    ref = cls[CIP_CUT]
    ever = np.zeros_like(ref)
    for c in CUT_SWEEP:
        ever |= (cls[c] != ref)
    return dict(per_cut={f"{c:.1f}": int(cls[c].sum()) for c in CUT_SWEEP},
                flips_vs_render_cut={f"{c:.1f}": int((cls[c] != ref).sum()) for c in CUT_SWEEP},
                n_anions_unstable=int(ever.sum()), n_anions=int(len(ref)))


def wrap_and_center(L, r, li, an, parent, n_ion):
    out = r.copy()
    for idx in [np.array([i]) for i in li] + [an[parent == p] for p in range(n_ion)]:
        c = out[idx].mean(axis=0)
        out[idx] -= np.floor(c / L) * L
        if not (out[idx].max(axis=0) - out[idx].min(axis=0) < L / 2).all():
            raise SystemExit("an ion is not intact in the dump; coordinates must be unwrapped")
    return out - L / 2.0


def write_xyz(path, L, r, li, an, parent, n_ion, paired):
    lines = [("LI", r[i]) for i in li]
    for p in range(n_ion):
        nm = "AP" if paired[p] else "AF"
        lines.extend((nm, r[k]) for k in an[parent == p])
    with open(path, "w") as f:
        f.write(f"{len(lines)}\n")
        f.write(f"cell {L[0]:.4f} {L[1]:.4f} {L[2]:.4f} centred on origin; "
                f"LI=Li+ AP=Li-contacted anion AF=free anion\n")
        for nm, xyz in lines:
            f.write(f"{nm} {xyz[0]:10.4f} {xyz[1]:10.4f} {xyz[2]:10.4f}\n")
    return len(lines)


def nspecies_summary(per_rep):
    out = {}
    for rep, s in sorted(per_rep.items()):
        md = s["md"]
        row = dict(n_kept=int(s["n_kept"]), n_raw=int(s["n_all"]), ts_first=int(s["ts"][0]), ts_last=int(s["ts"][-1]),
                   closest=float(md.min()), closest_median=float(np.median(md.min(axis=(1, 2)))))
        for c in CUT_SWEEP:
            v = s["nsp"][c]
            n_obs, cnt = np.unique(v, return_counts=True)
            row[f"{c:.1f}"] = dict(mean=float(v.mean()), median=float(np.median(v)), min=int(v.min()), max=int(v.max()),
                                   frac_ceiling=float((v == 16).mean()),
                                   hist={str(int(k)): int(x) for k, x in zip(n_obs, cnt)})
        out[rep] = row
    agg = {}
    for c in CUT_SWEEP:
        m = np.array([out[r][f"{c:.1f}"]["mean"] for r in out])
        agg[f"{c:.1f}"] = dict(mean=float(m.mean()), sd=float(m.std(ddof=1)) if len(m) > 1 else None,
                               lo=float(m.min()), hi=float(m.max()), n_replicas=int(len(m)))
    return dict(per_replica=out, across_replicas=agg,
                closest=min(out[r]["closest"] for r in out), closest_hi=max(out[r]["closest"] for r in out),
                closest_median_lo=min(out[r]["closest_median"] for r in out),
                closest_median_hi=max(out[r]["closest_median"] for r in out),
                ts_first=min(out[r]["ts_first"] for r in out), ts_last=max(out[r]["ts_last"] for r in out),
                frames_total=sum(out[r]["n_kept"] for r in out))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--traj", action="append", required=True, metavar="SALT:REPLICA=PATH",
                    help="LAMMPS ions dump (id mol type xu yu zu) of one replica, salt litpfa or lipf6")
    ap.add_argument("--drop", type=float, default=DROP, help="fraction of de-duplicated frames discarded as equilibration")
    ap.add_argument("--outdir", default=str(ROOT / "data" / "simulation" / "species_count"))
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    runs = {}
    for spec in a.traj:
        key, path = spec.split("=", 1)
        salt, rep = key.split(":", 1)
        runs.setdefault(salt, {})[rep] = path
    out = {"contact_cut": CIP_CUT, "cut_sweep": list(CUT_SWEEP), "drop": a.drop, "systems": {}}
    for salt, reps in runs.items():
        per_rep = {rep: scan(path, salt, a.drop) for rep, path in reps.items()}
        rep, idx, st = choose(per_rep)
        s = per_rep[rep]
        L, typ, r, li, an, parent = s["frames"][idx]
        n_ion = len(li)
        rw = wrap_and_center(L, r, li, an, parent, n_ion)
        paired = s["paired"][CIP_CUT][idx]
        xyz = os.path.join(a.outdir, f"{salt}_fm_snapshot.xyz")
        natoms = write_xyz(xyz, L, rw, li, an, parent, n_ion, paired)
        conc = n_ion / (6.02214076e23 * float(np.prod(L)) * 1e-27)
        rows = {}
        for k, v in per_rep.items():
            q1, q2, q3 = quartiles(v["nsp"][CIP_CUT])
            rows[k] = dict(n_raw=int(v["n_all"]), n_kept=int(v["n_kept"]), ts_first=int(v["ts"][0]), ts_last=int(v["ts"][-1]),
                           median=q2, q1=q1, q3=q3, mean=float(v["nsp"][CIP_CUT].mean()),
                           min=int(v["nsp"][CIP_CUT].min()), max=int(v["nsp"][CIP_CUT].max()),
                           by_cut={f"{c:.1f}": dict(median=float(np.median(v["nsp"][c])), mean=float(v["nsp"][c].mean()))
                                   for c in CUT_SWEEP})
        T = TYPES[salt]
        out["systems"][salt] = dict(
            replica=rep, frame_index=idx, xyz=os.path.basename(xyz), n_atoms=natoms, box=[float(x) for x in L], conc_M=conc,
            n_ions=2 * n_ion, n_pairs=n_ion, n_paired_anions=int(paired.sum()), paired_mask=[bool(x) for x in paired],
            anion_types=sorted(T["anion"]), li_type=T["li"], elements={str(k): v for k, v in sorted(ELEMENTS[salt].items())},
            selection=st, stability=stability(s, idx), replicas=rows, nspecies=nspecies_summary(per_rep))
        g = out["systems"][salt]["nspecies"]["across_replicas"][f"{CIP_CUT:.1f}"]
        print(f"{salt}: n_species at {CIP_CUT:.1f} A = {g['mean']:.2f} +/- {g['sd']:.2f} over {g['n_replicas']} replicas; "
              f"snapshot {rep} frame {idx} (timestep {st['ts']}, n_species {st['nsp']}) -> {xyz}")
    jf = os.path.join(a.outdir, "species_count.json")
    json.dump(out, open(jf, "w"), indent=1)
    print(f"-> {jf}")


if __name__ == "__main__":
    main()
