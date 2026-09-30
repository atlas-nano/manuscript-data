#!/usr/bin/env python3
"""Regenerate the samQEq parameter files used for the ReaxFF port from the ReaxFF ffields.

ReaxFF: eta = 2*etaEEM (reaxff_ffield.cpp:200); fix qeq/reaxff uses chi, eta in eV with a 14.4 eV*A
Coulomb constant. samQEq (units real) reads chi, eta in kcal/mol and uses qqrd2e = 332.06371.
  opmatch: x 332.06371/14.4 -> operator identical to qeq/reaxff's
  phys:    x 23.060549      -> samQEq's own eV->kcal/mol (fix_qeq_base_sam.cpp:108)
Writes to ./regen/ so the files actually used are never overwritten.
"""
import os

def eem(ffield, elements):
    lines = open(ffield).read().splitlines()
    i = next(k for k, l in enumerate(lines) if 'Nr of atoms' in l)
    nat = int(lines[i].split()[0]); i += 4
    out = {}
    for _ in range(nat):
        l1, l2 = lines[i].split(), lines[i + 1].split()
        out[l1[0]] = (float(l2[5]), 2.0 * float(l2[6]), float(l1[6]))  # chi, eta, gamma
        i += 4
    return [out[e] for e in elements]

SCALES = {'opmatch': 332.06371 / 14.4, 'phys': 23.060549}
os.makedirs('regen', exist_ok=True)
for name, ff, els in (('sam', 'qeq_ff.water', 'OH'), ('rdx_sam', 'ffield.reax', 'CHON')):
    for tag, s in SCALES.items():
        with open(f'regen/{name}_{tag}.param', 'w') as f:
            f.write('# gamma_align kappa_bond r_ov r_loc lr_alpha lr_ewald  (units real)\n0.30  0.0  2.0  2.5  0.0  2\n')
            for t, (el, (chi, eta, g)) in enumerate(zip(els, eem(ff, els)), 1):
                f.write(f'{t}  {chi*s:.12f}  {eta*s:.12f}  {g:.4f}  0.0  {-7*s:.6f}  {3*s:.6f}   # {el}\n')
