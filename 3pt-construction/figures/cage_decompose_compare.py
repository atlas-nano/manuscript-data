#!/usr/bin/env python3
"""Decompose the 3PT cage construction across material classes at matched
(fine-dump) resolution. For each state compute:
  f          fluidicity (gas DoF fraction)
  Scage      dS_cage delivered [k_B/DoF]  (= p*gate*int cage(1-w)(Wg-Ws))
  intcage    int cage d_nu       [cage DoF, per atom]
  cagepk     cage-band peak frequency [cm^-1]
  Om0,gam    Einstein freq & friction Ktilde(0) [cm^-1]; gam/Om0 = non-Markovianity
  Wg         gas weight per DoF
  WgWs       mean (Wg-Ws) over the cage band
and (where a reference is available) the TARGET gap = ref - rigHS the cage must close.
"""
import numpy as np
from kernel_tools import (cage_memory_entropy, packing_from_f_dgen,
                          hs_excess_entropy_dgen, invert_kernel_matrix,
                          H, VLIGHT, KB, PI, NA)
c = VLIGHT * 1e-10
hc_k = 100.0 * H * VLIGHT / KB


def diag(npz):
    d = np.load(npz)
    nu = d['nu_cm'].astype(float); tot = d['dos_total'].astype(float)
    gas = d['dos_gas'].astype(float); C = d['C_scalar'].astype(float)
    T = float(d['T_K']); D = int(d['dimension']); dnu = nu[1] - nu[0]
    wg_ov = float(d['wsr']) if 'wsr' in d.files else None
    co = []
    dS = cage_memory_entropy(float(d['dt']), C, nu, tot, gas, T,
                             float(d['mass_amu']) if 'mass_amu' in d.files else 18.0,
                             float(d['vol_per_atom']) if 'vol_per_atom' in d.files else 30.0,
                             dimension=D, prefactor=1.0 / D, gate_f0=0.01, nuc_scale=1.0,
                             Wg_override=wg_ov, cage_out=co)
    if dS is None or not co:
        return None
    cage = co[0]
    u = np.where(nu > 0, hc_k * nu / T, 1e-9); Ws = np.where(nu > 0, 1.0 - np.log(u), 0.0)
    f = float(np.trapezoid(gas, dx=dnu) / D)
    if wg_ov is not None:
        Wg = wg_ov
    else:
        y = packing_from_f_dgen(f, D); m = float(d['mass_amu']) * 1e-3 / NA
        lam = (2 * PI * m * KB * T / H**2)**(D / 2.0); V = float(d['vol_per_atom']) * (1e-10)**D
        Wg = ((D / 2.0 + 1.0 + np.log(lam * V / f)) / D + hs_excess_entropy_dgen(y, D) / D)
    intcage = float(np.trapezoid(cage, dx=dnu))
    cagepk = float(nu[np.argmax(cage)]) if cage.max() > 0 else 0.0
    wa = 2 * PI * nu * c
    Om0 = np.sqrt(np.trapezoid(wa**2 * tot, dx=1.0) / np.trapezoid(tot, dx=1.0)) / (2 * PI * c)
    cn = C / C[0]; Cm = np.zeros((cn.size, 3, 3))
    for i in range(3):
        Cm[:, i, i] = cn
    K = np.einsum('tii->t', invert_kernel_matrix(float(d['dt']), Cm)) / 3.0
    gam = float(np.trapezoid(K, dx=float(d['dt']))) / (2 * PI * c)
    wgws = float(np.trapezoid(cage * (Wg - Ws), dx=dnu) / max(intcage, 1e-12))
    return dict(f=f, dS=dS, intcage=intcage, cagepk=cagepk, Om0=Om0, gam=gam,
                nonmark=gam / Om0, Wg=Wg, WgWs=wgws, D=D)
