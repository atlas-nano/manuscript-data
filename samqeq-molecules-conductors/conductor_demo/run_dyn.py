#!/usr/bin/env python3
"""Referee item 3 (dynamics): the Au(111) + TIP4P-FQ film of experiment I under thermostatted dynamics.

Rigid water (fix rigid/nvt/small keyed on a per-atom body label, so that the global-constraint arms can set one
fragment for the solve without merging the rigid bodies), frozen gold, 298 K, 0.2 fs (the shipped TIP4P-FQ timestep; at 1 fs the metal-free film itself
goes near-critical), NSTEP steps, charges solved every
step. An ILLUSTRATIVE Au-O Lennard-Jones term (sigma 3.0 A, epsilon 0.05 eV) keeps the film on the slab; it is not a
fitted gold-water potential: the test concerns the charge solve only. Output: thermo every 100 steps (max |q|,
CG iterations from SOLVER-DIAG), charges dumped every 1000 steps. Analysis: analyze_dyn.py.
"""
import os, subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import run_exp1

HERE = Path(__file__).resolve().parent
LMP = os.environ.get("LMP", str(Path.home() / "codes/lammps/lammps-release/build/lmp"))
NSTEP = int(os.environ.get("NSTEP", "25000"))
CASES = [("samqeq", "0.3"), ("samqeq", "5.172"), ("qeq", "5.172"), ("qeq", "0.3"), ("mol", "0.3"), ("nometal", "-")]

DECK = """# conductor demonstration, dynamics: arm {arm}, eta_Au {eta} (run_dyn.py)
units           metal
atom_style      full
boundary        p p f
fix             pb all property/atom i_body ghost yes
read_data       {data} fix pb NULL RigidBodies
group           au   type 1
group           water type 2 3 4
group           fq   type 1 3 4
pair_style      hybrid/overlay lj/cut 12.0 coul/long 12.0 coul/shield/intra 12.0 all gaussian
pair_coeff      * * lj/cut 0.0 1.0
pair_coeff      2 2 lj/cut 0.012412 3.159
pair_coeff      1 2 lj/cut 0.05 3.0            # illustrative Au-O term (not a fitted gold-water potential)
pair_coeff      * * coul/long
pair_coeff      1 1 coul/shield/intra 1.618000
pair_coeff      1 2 coul/shield/intra 0.714283
pair_coeff      1 3 coul/shield/intra 0.808855
pair_coeff      1 4 coul/shield/intra 0.714283
pair_coeff      2 2 coul/shield/intra 0.315328
pair_coeff      2 3 coul/shield/intra 0.357078
pair_coeff      2 4 coul/shield/intra 0.315328
pair_coeff      3 3 coul/shield/intra 0.404355
pair_coeff      3 4 coul/shield/intra 0.357078
pair_coeff      4 4 coul/shield/intra 0.315328
special_bonds   lj/coul 1.0e-8 1.0e-8 1.0
kspace_style    pppm/samqeq 1.0e-6
kspace_modify   slab 3.0 mesh 24 28 80 gewald 0.35
neighbor        2.0 bin
neigh_modify    every 1 delay 0 check yes
{cons}
fix             chg fq qeq/sam 1 0.0 10.0 1.0e-8 param
fix_modify      chg shield gaussian
fix_modify      chg energy yes
{gself}
velocity        water create 298.0 4928459 mom yes rot yes dist gaussian
fix             rig water rigid/nvt/small custom i_body temp 298.0 298.0 0.1
variable        aq atom abs(q)
compute         qmax fq reduce max v_aq
thermo_style    custom step pe c_qmax
thermo          500
timestep        0.0002
dump            dq all custom 2500 charges.lammpstrj id mol type q x y z
dump_modify     dq sort id
run             {nstep}
print           "DYN_DONE" file done.txt screen no
"""


def write_data(path):
    """data.exp1 with explicit image flags (molecules whole: without flags read_data wraps an H/M outside the box to the
    far side and fix rigid/small builds a split body) plus a RigidBodies section (body id = original molecule id)"""
    import math
    LX, LY, LZ = 17.65065805, 20.38122435, 60.0
    lines = (HERE / "data.exp1").read_text().splitlines()
    ia = next(k for k, l in enumerate(lines) if l.startswith("Atoms"))
    for k in range(ia + 2, len(lines)):
        p = lines[k].split()
        if not p: break
        x, y, z = (float(v) for v in p[4:7])
        ix, iy, iz = math.floor(x / LX), math.floor(y / LY), math.floor(z / LZ)
        lines[k] = "%s %s %s %s %.8f %.8f %.8f %d %d %d" % (p[0], p[1], p[2], p[3], x - ix * LX, y - iy * LY, z - iz * LZ, ix, iy, iz)
    s = "\n".join(lines) + "\n"
    i0 = next(k for k, l in enumerate(lines) if l.startswith("Atoms"))
    i1 = next(k for k, l in enumerate(lines) if l.startswith("Bonds"))
    body = ["RigidBodies", ""]
    for l in lines[i0 + 2:i1]:
        if l.strip():
            a = l.split(); body.append(f"{a[0]} {a[1]}")
    Path(path).write_text(s.rstrip("\n") + "\n\n" + "\n".join(body) + "\n")


def make(arm, eta):
    if arm == "nometal":          # control: the same film, gold deleted, per-molecule constraints
        d = make("samqeq", "5.172")
        c = HERE / "runs/dyn/nometal"; c.mkdir(parents=True, exist_ok=True)
        for f in ("data", "param"): (c / f).write_text((d / f).read_text())
        s = (d / "in.lammps").read_text()
        s = s.replace("group           fq   type 1 3 4", "delete_atoms    group au bond yes\ngroup           fq   type 3 4")
        s = s.replace("fix_modify      chg gself on width 0.50 types 1\n", "")
        (c / "in.lammps").write_text(s)
        return c
    cons, gs = run_exp1.ARMS[arm]
    d = HERE / f"runs/dyn/{arm}_eta{eta}"; d.mkdir(parents=True, exist_ok=True)
    write_data(d / "data")
    (d / "param").write_text(run_exp1.PARAM.format(eta=eta))
    (d / "in.lammps").write_text(DECK.format(arm=arm, eta=eta, data="data", nstep=NSTEP,
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
    ds = [make(a, e) for a, e in CASES]
    with ThreadPoolExecutor(len(ds)) as ex:
        for d, rc in ex.map(run, ds):
            print(("FAIL rc=%d " % rc if rc else "ok ") + d.name)
