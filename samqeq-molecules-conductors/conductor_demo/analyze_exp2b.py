#!/usr/bin/env python3
"""Experiment II, done properly (referee P3).

(1) Screening from charge DENSITIES, not site sums: each solved site carries a Gaussian density of the pair kernel's
    width sigma_i = (2 alpha_i)^{-1/2} (2.38 A for Au); the charge on the near side of the slab's mid-plane is
    sum_i q_i Phi((z_i - z_mid)/sigma_i). With the +1 charge above, a conductor holds exactly -1 above the mid-plane
    (zero field inside). Also reported with the self-energy width w = 0.5 A, and as bare site sums.
(2) Image energy against the EXACT lateral-lattice image sum: for a rectangular lateral cell of area A, a lattice of
    charges q at height z above an image plane z0 has, per charge, E = (k_e q^2 / 2)(2 pi / A)[d - sum_{G != 0} e^{-G d}/G]
    + const, d = 2 (z - z0)  (the image lattice of -q, 2D Poisson summation). Two charges (one per face) -> x2.
    Fit C and z0 on z >= ZFIT.
"""
import csv, math, re
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
KE = 14.3996454784
LX, LY = 3 * 17.65065805, 3 * 20.38122435
AREA = LX * LY
ZOFF = 15.0
ZBOT, ZTOP = 10.0 + ZOFF, 19.60780130 + ZOFF
ZMID = 0.5 * (ZBOT + ZTOP)
LAM, RC_AU, W = 0.462770, 1.618, 0.5
SIG_PAIR = 1.0 / math.sqrt(2 * LAM / (2 * RC_AU ** 2))      # (2 alpha)^-1/2 = 2.378 A
ZFIT = 4.0
PHI = lambda x: 0.5 * (1 + math.erf(x / math.sqrt(2)))

# reciprocal vectors of the lateral cell for the image sum
nG = 60
gx = 2 * math.pi / LX * np.arange(-nG, nG + 1)
gy = 2 * math.pi / LY * np.arange(-nG, nG + 1)
GX, GY = np.meshgrid(gx, gy)
G = np.sqrt(GX ** 2 + GY ** 2).ravel(); G = G[G > 0]


def e_image(z, z0):
    d = 2 * (z - z0)
    return 2 * (KE / 2) * (2 * math.pi / AREA) * (d - np.sum(np.exp(-G * d) / G))


def fit(zs, es):
    best = None
    for k in range(-2000, 3000):
        z0 = k * 1e-3
        if z0 > min(zs) - 1.0: break
        f = np.array([e_image(z, z0) for z in zs]); C = np.mean(np.array(es) - f)
        r = math.sqrt(np.mean((np.array(es) - C - f) ** 2))
        if best is None or r < best[2]: best = (z0, C, r)
    return best


def densities(dump):
    s_pair = s_w = s_site = 0.0; go = False
    for l in Path(dump).read_text().splitlines():
        if l.startswith("ITEM: ATOMS"): go = True; continue
        if go and l.strip():
            a = l.split()
            if a[2] != "1": continue
            q, z = float(a[3]), float(a[6])
            s_pair += q * PHI((z - ZMID) / SIG_PAIR)
            s_w += q * PHI((z - ZMID) / W)
            s_site += q if z > ZMID + 0.5 else (0.5 * q if abs(z - ZMID) < 0.5 else 0.0)
    return s_pair, s_w, s_site


if __name__ == "__main__":
    rows = []
    for arm in ["gself", "nogself"]:
        for eta in ["5.172", "0.3"]:
            pts = []
            for d in sorted((HERE / f"runs/exp2/{arm}_eta{eta}").glob("z*"), key=lambda p: float(p.name[1:])):
                if not (d / "pe.txt").exists(): continue
                z = float(d.name[1:]); e = float((d / "pe.txt").read_text())
                pts.append((z, e, *densities(d / "charges.dump")))
            fz = [p for p in pts if p[0] >= ZFIT]
            z0, C, r = fit([p[0] for p in fz], [p[1] for p in fz])
            print(f"{arm:8s} eta {eta:>5}: exact-image fit (z>={ZFIT}) z0 = {z0:+.3f} A (from the top plane: "
                  f"{z0 - 0:+.3f}), rms {1000*r:.1f} meV")
            for (z, e, sp, sw, ss) in pts:
                law = C + e_image(z, z0)
                print(f"   z {z:5.1f}  dev {1000*(e-law):+8.1f} meV   near-side charge: density(pair) {sp:+.4f}"
                      f"  density(w) {sw:+.4f}  sites {ss:+.4f}")
                rows.append(dict(arm=arm, eta=eta, z=z, E=e, law=law, dev_meV=1000 * (e - law), z0=z0, rms_meV=1000 * r,
                                 Q_density_pair=sp, Q_density_w=sw, Q_sites=ss))
    with open(HERE / "exp2b.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
