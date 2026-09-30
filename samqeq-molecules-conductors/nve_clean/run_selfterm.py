#!/usr/bin/env python3
"""Referee P10: how much the grid self-term moves the liquid dipole, at fixed configurations.

One trajectory of the 64-molecule SPC-FQ benchmark (the shipped settings: shielding cutoff 5 A, solve cutoff
10 A, pinned recip_self 0.37933), 20 ps NVT to equilibrate, then 50 ps NVT with positions dumped every 1 ps.
The 50 frames are re-solved (rerun) with three self-terms: the pin, the probe calibration (default) and the exact
per-atom route, at the fixed Ewald split and mesh the pin was measured at (PPPM otherwise re-tunes g_ewald
from the current charges at every run). The same frames in each arm, so the difference in the mean dipole carries no sampling noise.
Analysis: python3 run_selfterm.py --analyze.
"""
import os, re, subprocess, sys
from pathlib import Path
import numpy as np
import run_nve

HERE = Path(__file__).resolve().parent
LMP = os.environ.get("LMP", str(Path.home() / "codes/lammps/lammps-release/build/lmp"))
D = HERE / "runs_selfterm"
HEAD = run_nve.DECK.split("velocity")[0].format(name="selfterm", shield="5.0", shift="", tol="1.0e-7",
                                                extra="{extra}", dt="0.001", thermo_nvt=1000, thermo_nve=100)
TRAJ = HEAD.format(extra="fix_modify   chg recip_self 0.37933") + """velocity     all create 298 123456 mom yes rot yes dist gaussian
fix          integ all rigid/nvt/small molecule temp 298 298 0.1
thermo_style custom step temp pe c_dip
thermo       1000
run          20000
dump         d all custom 1000 traj.lammpstrj id type x y z
dump_modify  d sort id
run          50000
print        "TRAJ_DONE" file done.txt screen no
"""
ARMS = {"pin": "fix_modify   chg recip_self 0.37933", "probe": "", "peratom": "fix_modify   chg recip_self peratom"}
RERUN = """{head}variable     f loop 51
label        frame
variable     s equal 20000+1000*(v_f-1)
read_dump    ../traj/traj.lammpstrj ${{s}} x y z box yes
reset_timestep 0   # the fix re-solves at run setup only from step 0
run          0
write_dump   fq custom q.${{s}}.dump id type q x y z modify sort id
next         f
jump         SELF frame
print        "RERUN_DONE" file done.txt screen no
"""


def run(d):
    r = subprocess.run([LMP, "-in", "in.lammps", "-log", "log.lammps", "-screen", "none"], cwd=d, capture_output=True, text=True)
    (d / "run.out").write_text(r.stdout + r.stderr); return r.returncode


def mean_dipole(f):
    L = Path(f).read_text().splitlines()
    lo = [float(L[k].split()[0]) for k in (5, 6, 7)]; hi = [float(L[k].split()[1]) for k in (5, 6, 7)]
    box = np.array(hi) - np.array(lo)
    a = np.array([[float(x) for x in l.split()] for l in L[9:]])
    a = a[np.argsort(a[:, 0])]
    q, r = a[:, 2], a[:, 3:6]
    mu = []
    for k in range(0, len(a), 3):                      # O, H, H of one molecule (ids 3k+1..3k+3)
        d = r[k:k + 3] - r[k]; d -= box * np.round(d / box)
        mu.append(np.linalg.norm((q[k:k + 3, None] * d).sum(0)))
    return 4.80320 * np.mean(mu)


if "--analyze" in sys.argv:
    mu = {a: np.array([mean_dipole(D / a / f"q.{20000 + 1000 * i}.dump") for i in range(51)]) for a in ARMS}
    for a in ARMS: print(f"{a:<8} frames {len(mu[a])}  <mu> {mu[a].mean():.4f} D")
    for a in ("probe", "peratom"):
        d = mu[a] - mu["pin"]
        print(f"{a} - pin: {d.mean():+.4f} D ({100*d.mean()/mu['pin'].mean():+.2f}%), frame spread {d.std():.4f} D")
    sys.exit()

if __name__ == "__main__":
  for sub, deck in [("traj", TRAJ)] + [(a, RERUN.format(head=HEAD.format(extra=x).replace("kspace_style pppm/samqeq 1.0e-6", "kspace_style pppm/samqeq 1.0e-6\nkspace_modify gewald 0.33957131 mesh 18 18 18   # the split the pin was measured at"))) for a, x in ARMS.items()]:
    d = D / sub; d.mkdir(parents=True, exist_ok=True)
    for f in ("spcfq.mol", "spcfq_expt.param"): (d / f).write_text((HERE / f).read_text())
    (d / "in.lammps").write_text(deck)
    if (d / "done.txt").exists(): continue
    print(sub, "rc", run(d), flush=True)
