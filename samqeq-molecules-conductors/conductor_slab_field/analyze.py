#!/usr/bin/env python3
"""Layer-resolved induced charge on the Au(111) slab, eta sweep, gself on/off.

Layers are assigned by z from data.au_ortho (5 layers). The Gauss limit is the
charge a perfect conductor places on each face, eps0*E*A, with E = 0.5 V/A.
Writes layers.csv next to this script and prints a table.
"""
import os, glob, csv

HERE = os.path.dirname(os.path.abspath(__file__))
EPS0 = 0.00552635  # e / (V A)
EFIELD = 0.5       # V / A


def read_data(path):
    lines = open(path).read().splitlines()
    box = {}
    atoms = {}
    i = 0
    while i < len(lines):
        l = lines[i]
        for ax in ('x', 'y', 'z'):
            if l.strip().endswith(f'{ax}lo {ax}hi'):
                lo, hi = map(float, l.split()[:2])
                box[ax] = hi - lo
        if l.strip().startswith('Atoms'):
            i += 2
            while i < len(lines) and lines[i].strip():
                p = lines[i].split()
                atoms[int(p[0])] = float(p[6])  # full: id mol type q x y z
                i += 1
            continue
        i += 1
    return box, atoms


def read_q(path):
    q = {}
    on = False
    for l in open(path):
        if l.startswith('ITEM: ATOMS'):
            on = True
            continue
        if on:
            p = l.split()
            q[int(p[0])] = float(p[2])
    return q


box, zs = read_data(os.path.join(HERE, 'data.au_ortho'))
zvals = sorted(set(round(z, 2) for z in zs.values()))
assert len(zvals) == 5, zvals
layer = {i: zvals.index(round(z, 2)) + 1 for i, z in zs.items()}
gauss = EPS0 * EFIELD * box['x'] * box['y']

rows = []
for d in sorted(glob.glob(os.path.join(HERE, 'runs', 'eta*_gself*'))):
    name = os.path.basename(d)
    eta = float(name[3:9])
    g = name.split('gself')[1]
    q = read_q(os.path.join(d, 'charges.dump'))
    lay = [sum(v for i, v in q.items() if layer[i] == k) for k in range(1, 6)]
    pe = float(open(os.path.join(d, 'pe.txt')).read().split()[0])
    face = 0.5 * (abs(lay[0]) + abs(lay[4]))
    rows.append(dict(eta=eta, gself=g, L1=lay[0], L2=lay[1], L3=lay[2], L4=lay[3],
                     L5=lay[4], total=sum(lay), face_over_gauss=face / gauss, pe=pe))

rows.sort(key=lambda r: (r['gself'], -r['eta']))
with open(os.path.join(HERE, 'layers.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)
print(f"Gauss limit per face = {gauss:.4f} e  (A = {box['x']*box['y']:.3f} A^2)")
print(f"{'gself':>5} {'eta':>6} {'L1':>8} {'L2':>8} {'L3':>8} {'L4':>8} {'L5':>8} {'sum':>9} {'face/G':>7}")
for r in rows:
    print(f"{r['gself']:>5} {r['eta']:6.3f} {r['L1']:8.4f} {r['L2']:8.4f} {r['L3']:8.4f} "
          f"{r['L4']:8.4f} {r['L5']:8.4f} {r['total']:9.2e} {r['face_over_gauss']:7.3f}")
