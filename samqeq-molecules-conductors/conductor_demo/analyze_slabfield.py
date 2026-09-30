#!/usr/bin/env python3
"""The section-2.4 field test (ms_alphasigma_overleaf/conductor_slab) on the density definition used for the
point-charge test: charge above the slab mid-plane from Gaussian densities of the pair-kernel width (2.378 A),
relative to the Gauss limit eps0 E A = 0.9940 e; also the top-half site sum (layers 4+5 plus half of layer 3)."""
import math
from pathlib import Path
D = Path.home() / "manuscripts_todo/ms_alphasigma_overleaf/conductor_slab"
ZMID, SIG, GAUSS = 14.80390065, 2.378, 0.9940
PHI = lambda x: 0.5 * (1 + math.erf(x / math.sqrt(2)))
z = {}; go = False
for l in (D / "data.au_ortho").read_text().splitlines():
    s = l.split()
    if l.startswith("Atoms"): go = True; continue
    if go and len(s) >= 7: z[int(s[0])] = float(s[6])
print("eta  gself    density(pair)/Gauss  sites/Gauss")
out = []
for d in sorted(D.glob("runs/eta*_gself*")):
    q = {}; go = False
    for l in (d / "charges.dump").read_text().splitlines():
        if l.startswith("ITEM: ATOMS"): go = True; continue
        if go and l.strip(): a = l.split(); q[int(a[0])] = float(a[2])
    dens = sum(q[i] * PHI((z[i] - ZMID) / SIG) for i in q)
    site = sum(q[i] * (1.0 if z[i] > ZMID + 0.5 else 0.5 if abs(z[i] - ZMID) < 0.5 else 0.0) for i in q)
    eta, gs = d.name.split("_")
    print(f"{eta[3:]:>6} {gs[5:]:<8} {dens/GAUSS:+10.3f}           {site/GAUSS:+8.3f}")
