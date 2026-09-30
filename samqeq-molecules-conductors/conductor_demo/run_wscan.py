#!/usr/bin/env python3
"""Self-energy width scan (referee P4): w in WS, eta_Au in ETAS, on three systems, self-energy on (types 1):
  field : the section-2.4 Au slab in 0.5 V/A (ms_alphasigma_overleaf/conductor_slab, deck adapted to `gaussian`)
  exp1  : the Au + water interface, per-molecule constraints (samqeq arm)
  exp2  : the 3x3 slab with +-1 point charges at z = 6 A
plus a Lanczos lambda_min run for field and exp1. Writes runs/wscan/... ; analysis: analyze_wscan.py.
"""
import os, re, subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import run_exp1, run_exp2

HERE = Path(__file__).resolve().parent
LMP = os.environ.get("LMP", str(Path.home() / "codes/lammps/lammps-release/build/lmp"))
SLAB = Path.home() / "manuscripts_todo/ms_alphasigma_overleaf/conductor_slab"
WS = ["0.3", "0.5", "0.75", "1.0", "1.5", "2.0", "2.378"]
ETAS = ["5.172", "0.3"]


def field_deck(w, lmin):
    s = (SLAB / "in.lammps").read_text()
    s = s.replace("coul/shield/intra 12.0 all pqeq", "coul/shield/intra 12.0 all gaussian")
    s = s.replace("fix_modify      chg shield pqeq", "fix_modify      chg shield gaussian")
    s = s.replace("fix_modify      chg gself on width 0.50", f"fix_modify      chg gself on width {w}")
    s = s.replace("read_data       data.au_ortho", f"read_data       {SLAB/'data.au_ortho'}")
    if lmin:
        s = s.replace("run             0", "fix_modify      chg ridge eig 1000.0 200\nrun             0")
    assert f"width {w}" in s and "gaussian" in s
    return s


def make():
    dirs = []
    for w in WS:
        for eta in ETAS:
            # field
            for lm in (False, True):
                d = HERE / f"runs/wscan/field_eta{eta}_w{w}" / ("lmin" if lm else "")
                d.mkdir(parents=True, exist_ok=True)
                (d / "in.lammps").write_text(field_deck(w, lm))
                (d / "au_pqeq.param").write_text(f"1  2.0  2.0  12.0  0.0  2  0.0\n1    5.3100   {eta}   1.6180  0.0  -7.0  3.0\n")
                dirs.append(d)
            # exp1 samqeq arm
            for lm in (False, True):
                d = HERE / f"runs/wscan/exp1_eta{eta}_w{w}" / ("lmin" if lm else "")
                d.mkdir(parents=True, exist_ok=True)
                (d / "param").write_text(run_exp1.PARAM.format(eta=eta))
                deck = run_exp1.DECK.format(arm="samqeq", eta=eta, nometal="",
                    cons="# per-molecule neutrality", gself=f"fix_modify      chg gself on width {w} types 1",
                    lmin="fix_modify      chg ridge eig 1000.0 200" if lm else "")
                (d / "in.lammps").write_text(deck.replace("../../../data.exp1", str(HERE / "data.exp1")))
                dirs.append(d)
            # exp2 at z = 6
            d = HERE / f"runs/wscan/exp2_eta{eta}_w{w}"
            d.mkdir(parents=True, exist_ok=True)
            run_exp2.write_data(d / "data", 6.0)
            (d / "param").write_text(run_exp2.PARAM.format(eta=eta))
            (d / "in.lammps").write_text(run_exp2.DECK.format(arm="gself", eta=eta, z="6.0", data="data",
                                         gself=f"fix_modify      chg gself on width {w}"))
            dirs.append(d)
    return dirs


def run(d):
    if (d / "pe.txt").exists(): return d, 0
    r = subprocess.run([LMP, "-in", "in.lammps", "-log", "log.lammps", "-screen", "none"], cwd=d,
                       capture_output=True, text=True)
    (d / "run.out").write_text(r.stdout + r.stderr)
    return d, r.returncode


if __name__ == "__main__":
    ds = make()
    with ThreadPoolExecutor(int(os.environ.get("NJ", "8"))) as ex:
        for d, rc in ex.map(run, ds):
            if rc: print("FAIL", rc, d.relative_to(HERE))
    print("done", len(ds))
