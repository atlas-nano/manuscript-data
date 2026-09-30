#!/usr/bin/env python3
"""Referee P1: TIP4P-FQ electronegativities are defined only up to a per-molecule constant. Shift the water's
chi (H and M) by DCHI and re-solve: under a global constraint the water->metal transfer follows the offset;
under per-molecule constraints the charges must not change. Arms qeq, qeq_gself (eta 5.172, 0.05 where
definite), samqeq (5.172). Writes runs/offset/...; analysis inline (offset.csv)."""
import csv, os, subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import run_exp1, analyze_exp1 as a1
HERE = Path(__file__).resolve().parent
LMP = os.environ.get("LMP", str(Path.home() / "codes/lammps/lammps-release/build/lmp"))
DCHI = ["-2.0", "-1.0", "-0.5", "0.0", "0.5", "1.0", "2.0"]
CASES = [("qeq", "5.172"), ("qeq_gself", "5.172"), ("qeq_gself", "0.05"), ("samqeq", "5.172")]
PARAM = """0.30  10.0  2.0  3.5  0.0  2
1  5.3100  {eta}  1.6180  0.0  -7.0  3.0
2  0.00  100.00  0.315328  0.0  -7.0  3.0
3  {chH:.4f}   15.31  0.404355  0.0  -7.0  3.0
4  {chM:.4f}   16.11  0.315328  0.0  -7.0  3.0
"""
def make(arm, eta, dc):
    cons, gs = run_exp1.ARMS[arm]
    d = HERE / f"runs/offset/{arm}_eta{eta}_d{dc}"; d.mkdir(parents=True, exist_ok=True)
    (d / "param").write_text(PARAM.format(eta=eta, chH=4.50 + float(dc), chM=7.47 + float(dc)))
    deck = run_exp1.DECK.format(arm=arm, eta=eta, nometal="",
        cons="set             group all mol 1" if cons == "global" else "# per-molecule",
        gself="fix_modify      chg gself on width 0.50 types 1" if gs else "", lmin="")
    (d / "in.lammps").write_text(deck.replace("../../../data.exp1", str(HERE / "data.exp1")))
    return d
def run(d):
    if (d / "pe.txt").exists(): return d, 0
    r = subprocess.run([LMP, "-in", "in.lammps", "-log", "log.lammps", "-screen", "none"], cwd=d, capture_output=True, text=True)
    return d, r.returncode
if __name__ == "__main__":
    ds = [make(a, e, dc) for a, e in CASES for dc in DCHI]
    with ThreadPoolExecutor(8) as ex:
        for d, rc in ex.map(run, ds):
            if rc: print("FAIL", d)
    mol0 = a1.original_mol(); rows = []
    for a, e in CASES:
        for dc in DCHI:
            r = a1.analyse(HERE / f"runs/offset/{a}_eta{e}_d{dc}", mol0)
            rows.append(dict(arm=a, eta=e, dchi=float(dc), qwater=r["qwater"], maxnet=r["maxnet"]))
            print(f"{a:<10}{e:>6}  dchi {dc:>5}  Q_water {r['qwater']:+8.4f}  max|Q_mol| {r['maxnet']:.4f}")
    with open(HERE / "offset.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
