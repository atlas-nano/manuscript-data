#!/usr/bin/env python3
"""Dynamics of the conductor demonstration (run_dyn.py): per arm, max |q| over the trajectory (thermo), ridge
engagements (last '[N engaged solves so far]' counter the log printed; LAMMPS stops printing warnings after 100),
and from the charge dumps: net water charge Q_w, max |Q_mol| over water molecules, gold charge spread.
Writes dyn.csv."""
import csv, re
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent


def frames(p):
    L = Path(p).read_text().splitlines(); k = 0
    while k < len(L):
        n = int(L[k + 3]); a = np.array([[float(x) for x in l.split()] for l in L[k + 9:k + 9 + n]])
        yield int(L[k + 1]), a; k += 9 + n


rows = []
for d in sorted((HERE / "runs/dyn").iterdir()):
    if not (d / "done.txt").exists(): print(d.name, "not finished"); continue
    log = (d / "log.lammps").read_text()
    th = [l.split() for l in log.splitlines() if re.match(r"^\s+\d+\s+-?\d", l) and len(l.split()) == 3]
    q = np.array([float(t[2]) for t in th])
    eng = re.findall(r"at step (\d+) .*?\[(\d+) engaged", log)
    ncg = len(re.findall(r"CG did not converge", log))
    Qw, Qm, sAu = [], [], []
    for step, a in frames(d / "charges.lammpstrj"):
        au, w = a[a[:, 2] == 1], a[a[:, 2] != 1]
        Qw.append(w[:, 3].sum())
        mol = np.bincount(((w[:, 0] - 1) // 4).astype(int), weights=w[:, 3])   # water ids 4k+1..4k+4 (the global arms reset mol)
        Qm.append(np.abs(mol).max())
        sAu.append(au[:, 3].std() if len(au) else 0.0)
    r = dict(arm=d.name, qmax_max=q.max(), qmax_last=q[-1], ridge_last_step=int(eng[-1][0]) if eng else 0,
             ridge_count=int(eng[-1][1]) if eng else 0, cg_fail_reported=ncg, Qw_mean=np.mean(Qw[1:]),
             Qw_absmax=np.abs(Qw).max(), Qmol_max=max(Qm), sdAu=np.mean(sAu[1:]))
    rows.append(r)
    print(f"{d.name:<16} max|q| {r['qmax_max']:.3f} (last {r['qmax_last']:.3f})  ridge {r['ridge_count']:>5} by step "
          f"{r['ridge_last_step']:>5}  CGfail {ncg}  Q_w {r['Qw_mean']:+.3f} (|max| {r['Qw_absmax']:.3f})  "
          f"max|Q_mol| {r['Qmol_max']:.2e}  sd(q_Au) {r['sdAu']:.4f}")
with open(HERE / "dyn.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
