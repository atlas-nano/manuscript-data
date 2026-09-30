#!/usr/bin/env python3
"""Tabulate experiment I (runs/exp1) -> exp1.csv and a printed table.

Per arm and eta: lambda_min (Lanczos, lmin/ run), CG iterations, total water charge (= charge moved between
water and metal), max |net charge| of one molecule, max |q| of any atom, mean water dipole (D, about the
molecule's O; for a charged molecule this is origin-dependent and flagged), slab layer charges.
"""
import csv, math, re
from pathlib import Path

HERE = Path(__file__).resolve().parent
EA2D = 4.80320471          # e*A -> Debye
LX, LY = 17.65065805, 20.38122435
LAYERS = [10.0, 12.40195032, 14.80390065, 17.20585097, 19.60780130]


def read_dump(p):
    rows, go = [], False
    for l in Path(p).read_text().splitlines():
        if l.startswith("ITEM: ATOMS"): go = True; continue
        if go and l.strip():
            a = l.split()
            rows.append((int(a[0]), int(a[1]), int(a[2]), float(a[3]), *map(float, a[4:7])))
    return rows


def original_mol():
    """mol id from data.exp1 (the global arm overwrites mol ids with 1)."""
    m, go = {}, False
    for l in (HERE / "data.exp1").read_text().splitlines():
        s = l.strip()
        if s.startswith("Atoms"): go = True; continue
        if go and s.startswith("Bonds"): break
        if go and s:
            a = s.split(); m[int(a[0])] = (int(a[1]), int(a[2]))
    return m


def lmin(d):
    f = d / "lmin/log.lammps"
    if not f.exists(): return float("nan")
    m = re.search(r"Lanczos lambda_min=([-+0-9.eE]+)", f.read_text())
    return float(m.group(1)) if m else float("nan")


def cg_iters(d):
    m = re.search(r"SOLVER-DIAG step 0: \S+ (\d+) iters", (d / "log.lammps").read_text())
    return int(m.group(1)) if m else -1


def analyse(d, mol0):
    rows = read_dump(d / "charges.dump")
    wat = {}
    lay = [0.0] * 5
    qmax = 0.0
    for (i, _, t, q, x, y, z) in rows:
        qmax = max(qmax, abs(q))
        if t == 1:
            k = min(range(5), key=lambda j: abs(z - LAYERS[j])); lay[k] += q
        else:
            wat.setdefault(mol0[i][0], []).append((t, q, x, y, z))
    net = [sum(a[1] for a in at) for at in wat.values()]
    dips = []
    for at in wat.values():
        o = [a for a in at if a[0] == 2][0]
        # minimum image in x, y: the dump is wrapped, so a molecule can straddle the lateral boundary
        mx = lambda d: d - LX * round(d / LX); my = lambda d: d - LY * round(d / LY)
        px = sum(a[1] * mx(a[2] - o[2]) for a in at); py = sum(a[1] * my(a[3] - o[3]) for a in at)
        pz = sum(a[1] * (a[4] - o[4]) for a in at)
        dips.append(math.sqrt(px * px + py * py + pz * pz) * EA2D)
    return dict(qwater=sum(net), maxnet=max(abs(v) for v in net), qmax=qmax,
                mu=sum(dips) / len(dips), layers=lay, nwat=len(wat))


if __name__ == "__main__":
    mol0 = original_mol()
    out = []
    ref = analyse(HERE / "runs/exp1/samqeq_eta5.172_nometal", mol0)
    print(f"no-metal reference film: {ref['nwat']} waters, mean dipole {ref['mu']:.4f} D, max|q| {ref['qmax']:.3f}")
    print(f"{'arm':<10}{'eta':>7}{'lmin':>9}{'CG':>5}{'Q_water':>10}{'max|Qm|':>9}{'max|q|':>8}"
          f"{'mu/mu0':>8}   layers L1..L5 (bottom..top)")
    for arm in ["qeq", "qeq_gself", "mol", "samqeq"]:
        for eta in ["5.172", "2.5", "1.0", "0.3", "0.05"]:
            d = HERE / f"runs/exp1/{arm}_eta{eta}"
            r = analyse(d, mol0); lm = lmin(d); it = cg_iters(d)
            print(f"{arm:<10}{eta:>7}{lm:>9.3f}{it:>5}{r['qwater']:>10.4f}{r['maxnet']:>9.4f}{r['qmax']:>8.3f}"
                  f"{r['mu']/ref['mu']:>8.3f}   " + " ".join(f"{v:+.3f}" for v in r["layers"]))
            out.append(dict(arm=arm, eta=eta, lmin=lm, cg=it, qwater=r["qwater"], maxnet=r["maxnet"],
                            qmax=r["qmax"], mu=r["mu"], mu_ratio=r["mu"] / ref["mu"],
                            **{f"L{k+1}": v for k, v in enumerate(r["layers"])}))
    with open(HERE / "exp1.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
