"""Explicit coul/shield/intra coefficients for EVERY type pair (i <= j), gamma_ij = sqrt(gamma_i gamma_j).
Under pair_style hybrid/overlay LAMMPS does not mix a sub-style across types once `* * lj/cut` and `* * coul/long`
cover the cross pairs; a per-type-only list silently drops the shielding correction from every cross pair's FORCES
while the charge solve keeps it (RESULT_conductor_demo_20260927.md, ledger 2026-09-28)."""
import math


def shield_pairs(g):
    """g: {type: gamma}; returns the pair_coeff lines"""
    ts = sorted(g)
    return "\n".join(f"pair_coeff      {i} {j} coul/shield/intra {math.sqrt(g[i] * g[j]):.6f}"
                     for a, i in enumerate(ts) for j in ts[a:])
