#!/usr/bin/env python3
"""Harvest the Perlmutter CPC-1 scaling runs into one CSV.

Emits: tag,nodes,gpus,natoms,steps_per_s,katom_step_per_s,pair,bond,kspace,neigh,comm,output,modify,other
The percentages are the %total column of the LAMMPS MPI task timing breakdown.
'modify' is where the fix qeq/sam charge solve is accounted.
"""
import glob, os, re, sys

SEC = ["Pair", "Bond", "Kspace", "Neigh", "Comm", "Output", "Modify", "Other"]
rows = []
for log in sorted(glob.glob(os.path.expanduser("~/alphasigma_pm_bench/log.rep*_n*"))):
    tag = os.path.basename(log)[4:]
    m = re.match(r"rep(\d+)_n(\d+)", tag)
    if not m:
        continue
    rep, nodes = int(m.group(1)), int(m.group(2))
    txt = open(log, errors="replace").read()
    perf = re.search(r"Performance:.*?([\d.]+) timesteps/s, ([\d.]+) ([kM])atom-step/s", txt)
    nat = re.search(r"([\d,]+) atoms", txt)
    if not perf:
        print(f"# {tag}: no Performance line", file=sys.stderr)
        continue
    pct = {}
    for s in SEC:
        mm = re.search(r"^%s\s*\|.*\|\s+([\d.]+)\s*\r?$" % s, txt, re.M)
        pct[s] = mm.group(1) if mm else ""
    out = re.search(r"natoms=(\d+)", open(log.replace("log.", "out."), errors="replace").read()) \
        if os.path.exists(log.replace("log.", "out.")) else None
    natoms = out.group(1) if out else (nat.group(1).replace(",", "") if nat else "")
    ks = float(perf.group(2)) * (1000.0 if perf.group(3) == "M" else 1.0)
    it = re.findall(r"CG (\d+) iters", txt)
    cg = it[-1] if it else ""
    rows.append([tag, nodes, nodes * 4, natoms, perf.group(1), f"{ks:.1f}", cg] +
                [pct[s] for s in SEC])

print("tag,nodes,gpus,natoms,steps_per_s,katom_step_per_s,cg_iters," + ",".join(s.lower() for s in SEC))
for r in sorted(rows, key=lambda r: (r[3], r[1])):
    print(",".join(str(x) for x in r))
