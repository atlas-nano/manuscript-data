#!/usr/bin/env python3
"""Experiment II of DESIGN.md: does the slab obey the image law?

The Au(111) slab replicated 3x3 laterally (2160 atoms). A fixed +1 point charge sits a height z above the top
atomic plane and a fixed -1 the same height below the bottom plane (neutral cell), both at the lateral centre
and both OUTSIDE the solve group, so they enter only through the fixed-charge field. Arms: self-energy on/off
x eta_Au in {5.172, 0.3}; z in HEIGHTS. Writes runs/exp2/<arm>_eta<eta>/z<z>/ and runs them.
"""
import os, subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
LMP = os.environ.get("LMP", str(Path.home() / "codes/lammps/lammps-release/build/lmp"))
AU = Path.home() / "manuscripts_todo/ms_alphasigma_overleaf/conductor_slab/data.au_ortho"
LX1, LY1, LZ = 17.65065805, 20.38122435, 60.0
NREP = 3
ZOFF = 15.0                    # slab shifted up so the lower charge stays inside the box at z = 12
Z_BOT, Z_TOP = 10.0 + ZOFF, 19.60780130 + ZOFF
HEIGHTS = ["2.0", "2.5", "3.0", "4.0", "5.0", "6.0", "8.0", "10.0", "12.0"]
ETAS = ["5.172", "0.3"]
ARMS = {"nogself": False, "gself": True}

PARAM = """1  2.0  2.0  12.0  0.0  2  0.0
1  5.3100  {eta}  1.6180  0.0  -7.0  3.0
2  0.0000  10.000  0.5000  0.0  -7.0  3.0
"""

DECK = """# conductor demonstration, experiment II: image law. arm {arm}, eta_Au {eta}, z {z} A (run_exp2.py)
units           metal
atom_style      full
boundary        p p f
read_data       {data}
group           au   type 1
group           ion  type 2
pair_style      hybrid/overlay lj/cut 12.0 coul/long 12.0 coul/shield/intra 12.0 all gaussian
pair_coeff      * * lj/cut 0.0 1.0
pair_coeff      * * coul/long
pair_coeff      1 1 coul/shield/intra 1.618000
pair_coeff      1 2 coul/shield/intra 0.899444
pair_coeff      2 2 coul/shield/intra 0.500000
kspace_style    pppm/samqeq 1.0e-6
kspace_modify   slab 3.0 mesh 72 84 80 gewald 0.35
neighbor        2.0 bin
fix             chg au qeq/sam 1 0.0 10.0 1.0e-8 param
fix_modify      chg shield gaussian
fix_modify      chg energy yes
{gself}
run             0
write_dump      all custom charges.dump id mol type q x y z modify sort id format line "%d %d %d %.12g %.8f %.8f %.8f"
variable        peval equal pe
print           "${{peval}}" file pe.txt screen no
"""


def write_data(path, z):
    au = []
    for l in AU.read_text().splitlines():
        p = l.split()
        if len(p) >= 7 and p[2] == "1" and p[0].isdigit():
            au.append(tuple(map(float, p[4:7])))
    atoms = []
    for ix in range(NREP):
        for iy in range(NREP):
            for (x, y, zz) in au:
                atoms.append((1, 1, 0.0, x + ix * LX1, y + iy * LY1, zz + ZOFF))
    cx, cy = NREP * LX1 / 2, NREP * LY1 / 2
    atoms.append((2, 2, +1.0, cx, cy, Z_TOP + z))
    atoms.append((3, 2, -1.0, cx, cy, Z_BOT - z))
    with open(path, "w") as f:
        f.write(f"# Au(111) 5L slab 3x3 + fixed +-1 point charges at {z} A (run_exp2.py)\n\n")
        f.write(f"{len(atoms)} atoms\n2 atom types\n\n")
        f.write(f"0.0 {NREP*LX1:.8f} xlo xhi\n0.0 {NREP*LY1:.8f} ylo yhi\n0.0 {LZ:.8f} zlo zhi\n\n")
        f.write("Masses\n\n1 196.966569\n2 22.99\n\nAtoms # full\n\n")
        for k, (m, t, q, x, y, zz) in enumerate(atoms, 1):
            f.write("%d %d %d %.6f %.8f %.8f %.8f\n" % (k, m, t, q, x, y, zz))


def make(arm, eta, z):
    d = HERE / "runs/exp2" / f"{arm}_eta{eta}" / f"z{z}"
    d.mkdir(parents=True, exist_ok=True)
    write_data(d / "data", float(z))
    (d / "param").write_text(PARAM.format(eta=eta))
    (d / "in.lammps").write_text(DECK.format(arm=arm, eta=eta, z=z, data="data",
        gself="fix_modify      chg gself on width 0.50" if ARMS[arm] else ""))
    return d


def run(d):
    if (d / "pe.txt").exists():
        return d, 0
    r = subprocess.run([LMP, "-in", "in.lammps", "-log", "log.lammps", "-screen", "none"], cwd=d,
                       capture_output=True, text=True)
    (d / "run.out").write_text(r.stdout + r.stderr)
    return d, r.returncode


if __name__ == "__main__":
    dirs = [make(a, e, z) for a in ARMS for e in ETAS for z in HEIGHTS]
    with ThreadPoolExecutor(int(os.environ.get("NJ", "6"))) as ex:
        for d, rc in ex.map(run, dirs):
            if rc: print(f"FAIL rc={rc} {d.relative_to(HERE)}")
    print("done", len(dirs))
