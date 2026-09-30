#!/usr/bin/env python3
"""Referee item 4: energy conservation of the full solve on a correctly configured deck.

Base: the student S4 deck (64 SPC-FQ waters, spcfq_expt.param, rigid molecules, 1 fs), which carried coul/shield/intra
5.0 against a solve cutoff of 10.0. Here the shielding cutoff matches the solve cutoff (10.0), the LJ is shifted to zero
at its cutoff, the density is fixed (the 64-molecule sc 3.1 box, 1.00 g/cc), and the protocol is 100 ps NVT then
200 ps NVE. Arms: tolerance 1e-5 / 1e-7 / 1e-9 at 1 fs, 1e-9 at 0.5 fs, the original mismatched deck (control), and
the predictor-corrector (aspc 2) on the matched deck. Release binary, one rank each. Analysis: analyze_nve.py.
"""
import os, subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
LMP = os.environ.get("LMP", str(Path.home() / "codes/lammps/lammps-release/build/lmp"))
NVT_PS, NVE_PS = float(os.environ.get("NVT_PS", "100")), float(os.environ.get("NVE_PS", "200"))

ARMS = {  # name: (shield cutoff, lj shift, tolerance, dt ps, extra fix_modify)
    "matched_tol1e-5":        ("10.0", True,  "1.0e-5", 0.001,  ""),
    "matched_tol1e-7":        ("10.0", True,  "1.0e-7", 0.001,  ""),
    "matched_tol1e-9":        ("10.0", True,  "1.0e-9", 0.001,  ""),
    "matched_tol1e-9_dt0.5":  ("10.0", True,  "1.0e-9", 0.0005, ""),
    "mismatched_tol1e-5":     ("5.0",  False, "1.0e-5", 0.001,  ""),
    "matched_aspc_tol1e-7":   ("10.0", True,  "1.0e-7", 0.001,  "fix_modify   chg aspc 2"),
    "matched_warm_tol1e-7":   ("10.0", True,  "1.0e-7", 0.001,  "fix_modify   chg solver warmstart on"),
    "matched_xl_tol1e-7":     ("10.0", True,  "1.0e-7", 0.001,  "fix_modify   chg xl 5.0e-5 0.5"),
}

DECK = """# referee item 4: NVE of SPC-FQ water, arm {name} (run_nve.py)
units        metal
atom_style   full
boundary     p p p
lattice      sc 3.1
region       box block 0 4 0 4 0 4
create_box   2 box
mass         1 15.9994
mass         2 1.008
molecule     h2o spcfq.mol
create_atoms 0 box mol h2o 12345
group        fq type 1 2
set          type 1 charge -0.85
set          type 2 charge  0.425
pair_style   hybrid/overlay lj/cut 10.0 coul/long 10.0 coul/shield/intra {shield} all
pair_coeff   1 1 lj/cut 0.01275 3.176
pair_coeff   1 2 lj/cut 0.0 1.0
pair_coeff   2 2 lj/cut 0.0 1.0
pair_coeff   * * coul/long
pair_coeff   * * coul/shield/intra 1.11
{shift}
kspace_style pppm/samqeq 1.0e-6
neighbor     2.0 bin
neigh_modify every 1 delay 0 check yes
fix          chg fq qeq/sam 1 0.0 10.0 {tol} spcfq_expt.param
fix_modify   chg energy yes
{extra}
compute      dip fq dipole/samqeq
timestep     {dt}
velocity     all create 298 123456 mom yes rot yes dist gaussian
fix          integ all rigid/nvt/small molecule temp 298 298 0.1
thermo_style custom step temp pe etotal c_dip
thermo       {thermo_nvt}
run          {nvt}
unfix        integ
reset_timestep 0
fix          integ all rigid/nve/small molecule
compute      rt all temp/rigid ...
thermo       {thermo_nve}
run          {nve}
print        "NVE_DONE" file done.txt screen no
"""


def make(name):
    shield, shift, tol, dt, extra = ARMS[name]
    d = HERE / "runs" / name; d.mkdir(parents=True, exist_ok=True)
    for f in ("spcfq.mol", "spcfq_expt.param"):
        (d / f).write_text((HERE / f).read_text())
    deck = DECK.format(name=name, shield=shield, shift="pair_modify  pair lj/cut shift yes" if shift else "",
                       tol=tol, dt=dt, extra=extra, nvt=int(round(NVT_PS / dt)), nve=int(round(NVE_PS / dt)),
                       thermo_nvt=int(round(1.0 / dt)), thermo_nve=int(round(0.1 / dt)))
    deck = deck.replace("compute      rt all temp/rigid ...\n", "")
    (d / "in.lammps").write_text(deck)
    return d


def run(d):
    if (d / "done.txt").exists(): return d, 0
    r = subprocess.run([LMP, "-in", "in.lammps", "-log", "log.lammps", "-screen", "none"], cwd=d,
                       capture_output=True, text=True)
    (d / "run.out").write_text(r.stdout + r.stderr)
    return d, r.returncode


if __name__ == "__main__":
    ds = [make(n) for n in ARMS]
    with ThreadPoolExecutor(len(ds)) as ex:
        for d, rc in ex.map(run, ds):
            print(("FAIL rc=%d " % rc if rc else "ok ") + d.name)
