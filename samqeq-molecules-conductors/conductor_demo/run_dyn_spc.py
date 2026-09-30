#!/usr/bin/env python3
"""Referee item 3, dynamics, with the paper's benchmark water (SPC-FQ, cube-root kernel gamma 1.11, spcfq_expt.param).

Why not TIP4P-FQ: its M site sits at the positivity margin in every kernel (Slater 2s self-Coulomb 16.1 eV vs
eta_M 16.11 eV; Gaussian port 17.53 vs 16.11), and a metal-free film of the Gaussian port over-polarizes in dynamics
(median q_M -1.37 e vs -1.15 in the bulk liquid). See RESULT_conductor_demo_20260927.md.

System: the section-2.4 Au(111) slab (cube-root kernel, gamma_Au = 0.2372 1/A, i.e. the same contact value
3.416 eV as its Gaussian) + 180 SPC-FQ molecules in the gap between the slab and its periodic image (fully periodic, no vacuum:
a vacuum-bounded film of the same water heats and collapses from its grid start, with or without gold), random
orientations; the metal-free control is the same 180 molecules as bulk water. Solve and shielding cutoffs
matched at 12 A. Rigid water (fix rigid/nvt/small on a per-atom body label, so the global-constraint arms can put
every atom in one fragment), frozen gold, 298 K, 1 fs, NEQ + NPROD steps. Illustrative Au-O LJ (sigma 3.0 A,
epsilon 0.05 eV), not a fitted gold-water potential. Arms: per molecule + self-energy (eta 0.3, 5.172), global
without self-energy (eta 5.172 definite, 0.3 indefinite), per molecule without self-energy (eta 0.3), metal-free film.
"""
import math, os, random, subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
LMP = os.environ.get("LMP", str(Path.home() / "codes/lammps/lammps-release/build/lmp"))
AU = Path.home() / "manuscripts_todo/ms_alphasigma_overleaf/conductor_slab/data.au_ortho"
LX, LY, ZTOP = 17.65065805, 20.38122435, 19.60780130
NZL, DZL = 6, 2.9                 # six water layers, 2.9 A apart
LZ_CELL = 30.3                    # the slab's periodic image confines the water (gap 20.7 A, ~0.95 g/cc)
LZ_BULK = 15.9                    # metal-free control: the same 180 molecules as bulk water, ~0.94 g/cc
NEQ, NPROD = int(os.environ.get("NEQ", "5000")), int(os.environ.get("NPROD", "15000"))
GEOM = [(0.0, 0.0, 0.0), (0.81650, 0.57735, 0.0), (-0.81650, 0.57735, 0.0)]    # spcfq.mol, O at origin

CASES = {  # name: (constraint, self-energy, eta, metal)
    "samqeq_eta0.3":   ("mol", True, "0.3", True),
    "samqeq_eta5.172": ("mol", True, "5.172", True),
    "qeq_eta5.172":    ("global", False, "5.172", True),
    "qeq_eta0.3":      ("global", False, "0.3", True),
    "mol_eta0.3":      ("mol", False, "0.3", True),
    "nometal":         ("mol", False, "5.172", False),
}


def rot(v, a, b, c):
    ca, sa, cb, sb, cc, sc = math.cos(a), math.sin(a), math.cos(b), math.sin(b), math.cos(c), math.sin(c)
    x, y, z = v
    y, z = y * ca - z * sa, y * sa + z * ca
    x, z = x * cb + z * sb, -x * sb + z * cb
    x, y = x * cc - y * sc, x * sc + y * cc
    return x, y, z


WATER = HERE / "runs/spc_bulk180/water180.data"     # equilibrated SPC-FQ, 176 molecules, 17.65 x 20.38 x 14.97 A
LZW = 14.97
GAPW = 3.0


def read_water():
    lines = WATER.read_text().splitlines()
    i0 = next(k for k, l in enumerate(lines) if l.startswith("Atoms"))
    mols = {}
    for l in lines[i0 + 2:]:
        p = l.split()
        if not p or not p[0].isdigit(): break
        aid, m, t = int(p[0]), int(p[1]), int(p[2])
        x, y, z = (float(v) for v in p[4:7]); ix, iy, iz = (int(v) for v in p[7:10])
        mols.setdefault(m, []).append((t, x + ix * LX, y + iy * LY, z + iz * LZW))
    out = []
    for m in sorted(mols):
        at = sorted(mols[m], key=lambda a: a[0])            # O first
        o = at[0]
        sx, sy, sz = (math.floor(o[1] / LX) * LX, math.floor(o[2] / LY) * LY, math.floor(o[3] / LZW) * LZW)
        out.append([(t, x - sx, y - sy, z - sz) for (t, x, y, z) in at])   # whole molecule, O inside the box
    return out


def build(path, metal=True):
    water = read_water()
    atoms, bodies, aid = [], [], 0
    zoff = ZTOP + GAPW if metal else 0.0
    for mid, mol in enumerate(water, start=1):
        for (t, x, y, z) in mol:
            aid += 1
            atoms.append((aid, mid, t + 1, 0.0, x, y, z + zoff)); bodies.append((aid, mid))
    nw = len(water)
    if metal:
        for l in AU.read_text().splitlines():
            p = l.split()
            if len(p) >= 7 and p[0].isdigit() and p[2] == "1":
                aid += 1
                atoms.append((aid, nw + 1, 1, 0.0, float(p[4]), float(p[5]), float(p[6]))); bodies.append((aid, 0))
    lz = (LZW + 2 * GAPW + (ZTOP - 10.0)) if metal else LZW   # the slab's periodic image confines the water
    with open(path, "w") as f:
        # NO bonds/angles: as in the shipped SPC-FQ decks, every intramolecular pair is a full shielded pair
        f.write(f"# Au(111) + {nw} SPC-FQ (run_dyn_spc.py)\n\n{len(atoms)} atoms\n\n3 atom types\n\n")
        f.write(f"0.0 {LX:.8f} xlo xhi\n0.0 {LY:.8f} ylo yhi\n0.0 {lz:.8f} zlo zhi\n\n")
        f.write("Masses\n\n1 196.966569\n2 15.9994\n3 1.008\n\nAtoms # full\n\n")
        # wrapped coordinates WITH image flags: without them read_data wraps an H that sits outside the box to the far
        # side, and fix rigid/small builds a split molecule (found 2026-09-28; it broke every earlier dynamics deck)
        for (aid_, m_, t_, q_, x, y, z) in atoms:
            ix, iy, iz = math.floor(x / LX), math.floor(y / LY), math.floor(z / lz)
            f.write("%d %d %d %.6f %.8f %.8f %.8f %d %d %d\n" % (aid_, m_, t_, q_, x - ix * LX, y - iy * LY, z - iz * lz, ix, iy, iz))
        f.write("\nRigidBodies\n\n" + "".join(f"{a} {b}\n" for a, b in bodies))
    return nw


DECK = """# Au(111) + SPC-FQ film dynamics, arm {name} (run_dyn_spc.py)
units           metal
atom_style      full
boundary        p p p
fix             pb all property/atom i_body ghost yes
read_data       data fix pb NULL RigidBodies
group           water type 2 3
group           fq    type {fqtypes}
pair_style      hybrid/overlay lj/cut 12.0 coul/long 12.0 coul/shield/intra 12.0 all cbrt
pair_coeff      * * lj/cut 0.0 1.0
pair_coeff      2 2 lj/cut 0.01275 3.144
pair_coeff      1 2 lj/cut 0.05 3.0            # illustrative Au-O term (not a fitted gold-water potential)
pair_coeff      * * coul/long
pair_coeff      1 1 coul/shield/intra 0.237200
pair_coeff      1 2 coul/shield/intra 0.513120
pair_coeff      1 3 coul/shield/intra 0.513120
pair_coeff      2 2 coul/shield/intra 1.110000
pair_coeff      2 3 coul/shield/intra 1.110000
pair_coeff      3 3 coul/shield/intra 1.110000
kspace_style    pppm/samqeq 1.0e-6
kspace_modify   gewald 0.35 mesh {mesh}   # pinned: all-zero start charges make the automatic mesh far too coarse
neighbor        2.0 bin
neigh_modify    every 1 delay 0 check yes
{metal_del}
{cons}
fix             chg fq qeq/sam 1 0.0 12.0 1.0e-7 param
fix_modify      chg shield cbrt
fix_modify      chg energy yes
{gself}
velocity        water create 298.0 4928459 mom yes rot yes dist gaussian
fix             rig water rigid/nvt/small custom i_body temp 298.0 298.0 0.02
compute         kew water ke
variable        Tw equal c_kew*2.0/((6*{nw}-3)*8.617333e-5)
variable        aq atom abs(q)
compute         qmax fq reduce max v_aq
thermo_style    custom step v_Tw pe c_qmax
thermo          100
timestep        0.00025
run             2000            # gentle start from the grid (tight coupling, 0.25 fs)
unfix           rig
fix             rig water rigid/nvt/small custom i_body temp 298.0 298.0 0.1
timestep        0.001
run             {neq}
dump            dq all custom 500 charges.lammpstrj id mol type q x y z
dump_modify     dq sort id
run             {nprod}
print           "DYN_DONE" file done.txt screen no
"""


def make(name):
    cons, gs, eta, metal = CASES[name]
    d = HERE / "runs/dyn_spc" / name; d.mkdir(parents=True, exist_ok=True)
    nw = build(d / "data", metal=metal)
    (d / "param").write_text(f"0.30  10.0  2.0  3.5  0.0  2\n1  5.3100  {eta}  0.2372  0.0  -7.0  3.0\n"
                             "2  6.7721  15.1145  1.11  0.0  -7.0  3.0\n3  4.50  16.1595  1.11  0.0  -7.0  3.0\n")
    (d / "in.lammps").write_text(DECK.format(
        name=name, nw=nw, neq=NEQ, nprod=NPROD, mesh='24 28 42' if metal else '24 28 20', fqtypes="1 2 3" if metal else "2 3",
        metal_del="" if metal else "",
        cons="set             group all mol 1   # one global constraint (bodies keep i_body)" if cons == "global"
             else "# per-molecule constraints",
        gself="fix_modify      chg gself on width 0.50 types 1" if gs else ""))
    return d


def run(d):
    if (d / "done.txt").exists(): return d, 0
    r = subprocess.run([LMP, "-in", "in.lammps", "-log", "log.lammps", "-screen", "none"], cwd=d,
                       capture_output=True, text=True)
    (d / "run.out").write_text(r.stdout + r.stderr)
    return d, r.returncode


if __name__ == "__main__":
    ds = [make(n) for n in CASES]
    with ThreadPoolExecutor(len(ds)) as ex:
        for d, rc in ex.map(run, ds):
            print(("FAIL rc=%d " % rc if rc else "ok ") + d.name)
