#!/usr/bin/env python3
"""Referee P11/P12: an uncertainty for the ported TIP4P-FQ liquid dipole (paper: 2.6404 D from one 3.75 ps run of
examples/tip4pfq_liquid/in.tip4pfq_xl_slater_lr1) and the effect of the thermostat coupling.

The shipped deck unchanged except: velocity seed, 20 ps equilibration at tau_T 0.02 ps (the shipped coupling),
then 100 ps production at tau_T = 0.02 or 0.1 ps; mean molecular dipole every 0.1 ps. Four seeds per coupling.
Run on Perlmutter (one interactive CPU node, 8 replicas x 16 ranks): `--stage` writes the decks, pm_drive.sh runs them. Analysis: python3 run_port.py --analyze (ten 10 ps blocks per run)."""
import os, re, subprocess, sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
LMP = os.environ.get("LMP", str(Path.home() / "codes/lammps/lammps-release/build/lmp"))
RUNS = {f"tau{t}_s{s}": (t, s) for t in ("0.02", "0.1") for s in (65491, 20260928, 4928459, 777)}


def deck(tau, seed):
    s = (HERE / "in.tip4pfq_xl_slater_lr1").read_text()
    s = s.split("\ntimestep ")[0] + "\n"
    return s + f"""timestep        0.00025
velocity        all create 298.0 {seed} mom yes rot yes dist gaussian
fix             rig all rigid/nvt/small molecule temp 298.0 298.0 0.02
thermo_style    custom step temp v_mu
thermo          400
run             80000
unfix           rig
fix             rig all rigid/nvt/small molecule temp 298.0 298.0 {tau}
reset_timestep  0
run             400000
print           "PORT_DONE" file done.txt screen no
"""


def run(name):
    d = HERE / "runs" / name; d.mkdir(parents=True, exist_ok=True)
    for f in ("tip4pfq_liquid.data", "water4.mol", "water4_liquid_slater_lr1.param"):
        (d / f).write_text((HERE / f).read_text())
    (d / "in.lammps").write_text(deck(*RUNS[name]))
    if "--stage" in sys.argv or (d / "done.txt").exists(): return name, 0
    r = subprocess.run(["mpirun", "-np", "4", LMP, "-in", "in.lammps", "-log", "log.lammps", "-screen", "none"], cwd=d)
    return name, r.returncode


if __name__ == "__main__":
    if "--analyze" in sys.argv:
        allm = {}
        for name in RUNS:
            L = (HERE / "runs" / name / "log.lammps").read_text().split("reset_timestep")[-1].splitlines()
            v = np.array([float(l.split()[2]) for l in L if re.match(r"^\s+\d+\s+[0-9.]+\s+[0-9.]+\s*$", l)])
            b = np.array([x.mean() for x in np.array_split(v[1:], 10)])
            allm.setdefault(RUNS[name][0], []).append(b)
            print(f"{name:<18} {len(v)} pts  <mu> {v[1:].mean():.4f} D  block sd {b.std(ddof=1):.4f}  se {b.std(ddof=1)/np.sqrt(10):.4f}")
        for t, bs in allm.items():
            b = np.concatenate(bs); print(f"tau {t}: <mu> {b.mean():.4f} +/- {b.std(ddof=1)/np.sqrt(len(b)):.4f} D ({len(b)} blocks)")
        sys.exit()
    with ThreadPoolExecutor(4) as ex:
        for n, rc in ex.map(run, RUNS): print(n, "rc", rc, flush=True)
