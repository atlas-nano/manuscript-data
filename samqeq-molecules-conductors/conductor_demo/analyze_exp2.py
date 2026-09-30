#!/usr/bin/env python3
"""Experiment II (runs/exp2): image law. For each arm: charge induced on the upper half of the slab by the +1
charge above it (a conductor gives -1), and the energy E(z) fitted to  E = C - 2 k_e / (4 (z - z0))
(two charges, each with its own image; k_e = 14.399645 eV A)
plus the uniform-sheet term of the lateral image lattice, 2 pi k_e q^2 (z - z0)/area per charge. Writes exp2.csv."""
import csv, math, re
from pathlib import Path

HERE = Path(__file__).resolve().parent
KE = 14.3996454784
AREA = 3 * 17.65065805 * 3 * 20.38122435
ZMID = (10.0 + 19.60780130) / 2 + 15.0 + 0.5   # upper half: top two layers (the middle layer excluded)


def dump_upper(p):
    s, go = 0.0, False
    for l in Path(p).read_text().splitlines():
        if l.startswith("ITEM: ATOMS"): go = True; continue
        if go and l.strip():
            a = l.split()
            if a[2] == "1" and float(a[6]) > ZMID: s += float(a[3])
    return s


def fit(zs, es):
    """least squares for C, z0 in E = C - A/(z - z0) + B (z - z0): A = 2 k_e/4 (two charges, each with its own
    image), B = 2 * 2 pi k_e / area (the uniform-sheet part of the lateral image lattice). Both fixed; grid on z0."""
    A = 2 * KE / 4
    B = 2 * 2 * math.pi * KE / AREA
    best = None
    for k in range(-3000, 1500):
        z0 = k * 1e-3
        if z0 >= min(zs) - 0.2: break
        f = [-A / (z - z0) + B * (z - z0) for z in zs]
        C = sum(e - g for e, g in zip(es, f)) / len(es)
        r = math.sqrt(sum((e - C - g) ** 2 for e, g in zip(es, f)) / len(es))
        if best is None or r < best[2]: best = (z0, C, r)
    return best


if __name__ == "__main__":
    rows = []
    for arm in ["gself", "nogself"]:
        for eta in ["5.172", "0.3"]:
            base = HERE / f"runs/exp2/{arm}_eta{eta}"
            pts = []
            for d in sorted(base.glob("z*"), key=lambda p: float(p.name[1:])):
                if not (d / "pe.txt").exists(): continue
                z = float(d.name[1:]); e = float((d / "pe.txt").read_text())
                it = re.search(r"SOLVER-DIAG step 0: \S+ (\d+) iters", (d / "log.lammps").read_text())
                pts.append((z, e, dump_upper(d / "charges.dump"), int(it.group(1)) if it else -1))
            if not pts: continue
            far = [p for p in pts if p[0] >= 3.0]
            z0, C, r = fit([p[0] for p in far], [p[1] for p in far])
            print(f"{arm:8s} eta {eta:>5}: image fit (z>=3) z0 = {z0:+.3f} A, rms {r*1000:.1f} meV")
            for (z, e, qu, it) in pts:
                law = C - 2 * KE / (4 * (z - z0)) + 2 * 2 * math.pi * KE / AREA * (z - z0)
                print(f"   z {z:5.1f}  E {e:12.5f}  law {law:12.5f}  dev {1000*(e-law):+8.1f} meV"
                      f"  Q_upper {qu:+.4f}  CG {it}")
                rows.append(dict(arm=arm, eta=eta, z=z, E=e, law=law, Q_upper=qu, cg=it, z0=z0, rms=r))
    with open(HERE / "exp2.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
