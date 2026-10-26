"""Write the pair term phi(r) of a Li EAM file in fs format as a LAMMPS pair_style table file on stdout, for example: python3 li_phi_table.py Li_v2.eam.fs > li_phi.table."""
import sys

import numpy as np
from scipy.interpolate import CubicSpline

N, RMIN, RMAX = 8000, 0.5, 4.563

lines = open(sys.argv[1]).readlines()
nrho, drho, nr, dr, rc = lines[4].split()
nrho, nr, dr = int(nrho), int(nr), float(dr)
values = np.array(" ".join(lines[6:]).split(), dtype=float)
rphi = values[nrho + nr:nrho + 2 * nr]
spline = CubicSpline(np.arange(nr) * dr, rphi)
r = np.linspace(RMIN, RMAX, N)
phi = spline(r) / r
force = -(spline(r, 1) / r - spline(r) / r ** 2)
print("LI_PHI")
print(f"N {N} R {RMIN} {RMAX}")
print()
for i in range(N):
    print(f"{i + 1} {r[i]:.8f} {phi[i]:.10e} {force[i]:.10e}")
