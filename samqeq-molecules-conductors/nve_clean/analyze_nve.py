#!/usr/bin/env python3
"""Energy drift of the NVE leg (run_nve.py). The second thermo block of each log (after reset_timestep) is the NVE
leg; drift = slope of etotal vs time, converted to K/ns with the rigid-body degrees of freedom of 64 molecules
(6N - 3 = 381); uncertainty = standard error of the slopes of ten contiguous segments (their standard deviation over sqrt(10));
the standard deviation itself is also written (seg_sd).
Also the mean temperature and dipole of the leg. Writes nve.csv."""
import csv, math, re
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
KB = 8.617333262e-5
DOF = 6 * 64 - 3


def blocks(log):
    out, cur, hdr = [], None, None
    for l in Path(log).read_text().splitlines():
        s = l.split()
        if s[:1] == ["Step"]:
            hdr = s; cur = []; out.append((hdr, cur)); continue
        if cur is not None:
            if s and s[0].lstrip("-").isdigit() and len(s) == len(hdr):
                try: cur.append([float(x) for x in s])
                except ValueError: cur = None
            elif s[:1] == ["Loop"]:
                cur = None
    return out


rows = []
for d in sorted((HERE / "runs").iterdir()):
    log = d / "log.lammps"
    partial = not (d / "done.txt").exists()
    if partial and "reset_timestep" not in log.read_text():
        print(f"{d.name:<26} not finished"); continue
    dt = float([l.split()[1] for l in (d / "in.lammps").read_text().splitlines() if l.startswith("timestep")][0])
    hdr, data = blocks(log)[-1]
    a = np.array(data); c = {h: i for i, h in enumerate(hdr)}
    t = a[:, c["Step"]] * dt / 1000.0                      # ns
    e = a[:, c["TotEng"]]
    slope = np.polyfit(t, e, 1)[0]                          # eV/ns
    segs = np.array_split(np.arange(len(t)), 10)
    ss = [np.polyfit(t[s], e[s], 1)[0] for s in segs]
    k = 2.0 / (DOF * KB)
    row = dict(arm=d.name + (' (partial)' if partial else ''), ns=t[-1] - t[0], drift_K_per_ns=slope * k, seg_sd=np.std(ss, ddof=1) * k, seg_se=np.std(ss, ddof=1) * k / math.sqrt(len(ss)),
               T_mean=a[:, c["Temp"]].mean(), mu_first=a[:10, c["c_dip"]].mean(), mu_last=a[-10:, c["c_dip"]].mean(),
               E_sd_meV=1000 * np.std(e - np.polyval(np.polyfit(t, e, 1), t)))
    txt = (d / "log.lammps").read_text(); nve = txt[txt.rfind("reset_timestep"):]
    it = [int(m) for m in re.findall(r"SOLVER-DIAG step \d+: CG (\d+) iters", nve)]
    rb = [float(m) for m in re.findall(r"resid/b=([0-9.e+-]+)", nve)]
    lt = re.findall(r"Loop time of ([0-9.]+) on \d+ procs for (\d+) steps", nve)
    row.update(matvec=np.mean(it) if it else float("nan"), resid_b=np.median(rb) if rb else float("nan"),
               ms_step=1000 * float(lt[-1][0]) / int(lt[-1][1]) if lt else float("nan"))
    rows.append(row)
    print(f"{d.name:<26} {row['ns']*1000:6.0f} ps  drift {row['drift_K_per_ns']:+8.2f} ± {row['seg_se']:5.1f} (sd {row['seg_sd']:4.0f}) K/ns"
          f"  <T> {row['T_mean']:6.1f}  mu {row['mu_first']:.3f} -> {row['mu_last']:.3f} D  sd(E) {row['E_sd_meV']:.1f} meV  mv {row['matvec']:.2f}  r/b {row['resid_b']:.1e}  {row['ms_step']:.2f} ms")
if rows:
    with open(HERE / "nve.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
