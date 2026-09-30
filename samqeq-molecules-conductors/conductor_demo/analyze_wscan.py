#!/usr/bin/env python3
"""Width scan table -> wscan.csv. Field: lambda_min (Lanczos), density face ratio (pair width, mid-plane),
sub-surface layer sign (staggered = L4 opposite to L5 face... reported as L5, L4 site charges). exp1: lambda_min,
water dipole ratio to the metal-free film (min image). exp2: near-side density charge at z = 6."""
import csv, math, re
from pathlib import Path
import analyze_exp1 as a1
import analyze_exp2b as a2
HERE = Path(__file__).resolve().parent
KE, LAM = 14.3996454784, 0.462770
PHI = lambda x: 0.5 * (1 + math.erf(x / math.sqrt(2)))
SLAB = Path.home() / "manuscripts_todo/ms_alphasigma_overleaf/conductor_slab"
zs = {}; go = False
for l in (SLAB / "data.au_ortho").read_text().splitlines():
    s = l.split()
    if l.startswith("Atoms"): go = True; continue
    if go and len(s) >= 7: zs[int(s[0])] = float(s[6])
LAY = [10.0, 12.40195032, 14.80390065, 17.20585097, 19.60780130]
def lmin(d):
    m = re.search(r"Lanczos lambda_min=([-+0-9.eE]+)", (d / "lmin/log.lammps").read_text()); return float(m.group(1))
mol0 = a1.original_mol()
mu0 = a1.analyse(HERE / "runs/exp1/samqeq_eta5.172_nometal", mol0)["mu"]
rows = []
print(f"{'eta':>6} {'w':>6} {'bracket':>8} | field: {'lmin':>7} {'face/G':>7} {'L4':>7} {'L5':>7} | exp1: {'lmin':>6} {'mu/mu0':>7} | exp2 Q(z=6)")
for eta in ["5.172", "0.3"]:
    for w in ["0.3", "0.5", "0.75", "1.0", "1.5", "2.0", "2.378"]:
        fd = HERE / f"runs/wscan/field_eta{eta}_w{w}"
        q = {}; go = False
        for l in (fd / "charges.dump").read_text().splitlines():
            if l.startswith("ITEM: ATOMS"): go = True; continue
            if go and l.strip(): a = l.split(); q[int(a[0])] = float(a[2])
        dens = sum(q[i] * PHI((zs[i] - 14.80390065) / 2.378) for i in q) / 0.9940
        lay = [0.0] * 5
        for i in q: lay[min(range(5), key=lambda j: abs(zs[i] - LAY[j]))] += q[i]
        e1 = HERE / f"runs/wscan/exp1_eta{eta}_w{w}"
        r1 = a1.analyse(e1, mol0)
        _, _, _ = None, None, None
        sp, sw, ss = a2.densities(HERE / f"runs/wscan/exp2_eta{eta}_w{w}/charges.dump")
        br = float(eta) + KE / (math.sqrt(math.pi) * float(w)) - KE * math.sqrt(2 * (LAM / (2 * 1.618 ** 2)) / math.pi)
        row = dict(eta=eta, w=w, bracket=br, field_lmin=lmin(fd), face_gauss=dens, L4=lay[3], L5=lay[4],
                   exp1_lmin=lmin(e1), mu_ratio=r1["mu"] / mu0, exp2_Qnear=sp)
        rows.append(row)
        print(f"{eta:>6} {w:>6} {br:8.3f} |        {row['field_lmin']:7.3f} {dens:7.3f} {lay[3]:+7.3f} {lay[4]:+7.3f} |"
              f"       {row['exp1_lmin']:6.3f} {row['mu_ratio']:7.4f} | {sp:+.4f}")
with open(HERE / "wscan.csv", "w", newline="") as f:
    wr = csv.DictWriter(f, fieldnames=list(rows[0])); wr.writeheader(); wr.writerows(rows)
