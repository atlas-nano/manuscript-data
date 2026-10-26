"""Reads the simulation tables that Figures 3 and S3 plot."""
from __future__ import annotations

import statistics as st

from paths import data

TABLE_A = data("simulation", "xpt_dgmix_anion_clamp_293K.csv")
TABLE_RAW = data("simulation", "xpt_dgmix_293K.csv")

BE_TABLE_S8 = {"FM":   [-26.1, -24.2, -22.5, -21.1],
               "Me₂O": [-37.0, -33.4, -30.5, -28.1],
               "DEE":  [-40.8, -36.4, -33.0, None],
               "ACN":  [-43.2, -39.5, -36.3, -33.6],
               "DME":  [-52.6, -44.5, None, None]}

NSPECIES_JSON = data("simulation", "snapshot_stats.json")
NSPECIES_CUT, NSPECIES_ROBUST_CUT = 3.0, 4.5


def nspecies_canonical(path=None):
    import json
    p = path or NSPECIES_JSON
    raw = json.load(open(p))
    key = {"litpfa": "LiTPFA", "lipf6": "LiPF6"}
    out = {}
    for tag, sysd in raw["systems"].items():
        ns = sysd["nspecies"]
        hc, rc = f"{NSPECIES_CUT:.1f}", f"{NSPECIES_ROBUST_CUT:.1f}"
        g = ns["across_replicas"][hc]
        reps = sorted(ns["per_replica"])
        out[key[tag]] = dict(
            mean=g["mean"], sd=g["sd"], lo=g["lo"], hi=g["hi"], n_replicas=g["n_replicas"],
            replicas={r: ns["per_replica"][r][hc]["mean"] for r in reps},
            sweep={c: ns["across_replicas"][c]["mean"] for c in sorted(ns["across_replicas"],
                                                                       key=float)},
            frac_ceiling={c: st.fmean(ns["per_replica"][r][c]["frac_ceiling"] for r in reps)
                          for c in sorted(ns["across_replicas"], key=float)},
            median_sweep={c: st.median([ns["per_replica"][r][c]["median"] for r in reps])
                          for c in sorted(ns["across_replicas"], key=float)},
            pn_replicas={c: {r: {int(k): v / ns["per_replica"][r]["n_kept"]
                                 for k, v in ns["per_replica"][r][c]["hist"].items()}
                             for r in reps}
                         for c in sorted(ns["across_replicas"], key=float)},
            closest=ns["closest"], closest_hi=ns["closest_hi"],
            closest_median_lo=ns["closest_median_lo"], closest_median_hi=ns["closest_median_hi"],
            ts_first=ns["ts_first"], ts_last=ns["ts_last"], frames=ns["frames_total"],
            traj=sysd["traj"])
    for salt, d in out.items():
        d["pn"] = {c: {n: st.fmean(d["pn_replicas"][c][r].get(n, 0.0)
                                   for r in d["pn_replicas"][c])
                       for n in sorted({k for r in d["pn_replicas"][c]
                                        for k in d["pn_replicas"][c][r]})}
                   for c in d["pn_replicas"]}
        for c, h in d["pn"].items():
            tot = sum(h.values())
    return out
N_IONS, N_PAIRS = 16, 8

DG_ASSOC_FM = {"LiPF6": (-24.1, 0.4), "LiTPFA": (+15.0, 0.7)}
DDG_ASSOC_FM = +39.1

LAB = {"LiTPFA": "LiTPFA", "LiPF6": "LiPF$_6$", "LiBF4": "LiBF$_4$",
       "LiFSI": "LiFSI", "LiTFSI": "LiTFSI"}
SALT_COL = {"LiTPFA": "#0072B2", "LiPF6": "#CC79A7", "LiBF4": "#E69F00",
            "LiFSI": "#56B4E9", "LiTFSI": "#009E73"}


def gated_diss_tds():
    import csv
    with open(data("tables", "crystal_dissociation.csv"), encoding="utf-8") as fh:
        diss = {r["salt"]: float(r["dG_diss"]) for r in csv.DictReader(fh)}
    with open(data("tables", "ideal_mixing.csv"), encoding="utf-8") as fh:
        tds = {r["solvent"].replace("Me2O", "Me₂O"): float(r["TdS_ideal"]) for r in csv.DictReader(fh)}
    return diss, tds


def _per_box(path):
    import csv
    out = {}
    with open(path, encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            row = {c: float(r[c]) for c in _COLS}
            row["sigma_ctl"] = float(r["sigma_ctl"]) if r["sigma_ctl"] else None
            out.setdefault(r["system"], []).append(row)
    return out


_COLS = ("c_ggla", "dH", "dS_Li", "dS_an", "dG_solv", "dG_mix", "band_lo", "band_hi")


def dgmix_systems(systems):
    A, RAW = _per_box(TABLE_A), _per_box(TABLE_RAW)
    out = {}
    for salt, solv in systems:
        key = f"{salt.lower()}_{solv.lower()}"
        a, r = A[key], RAW[key]
        sig = [x["sigma_ctl"] for x in a if x["sigma_ctl"] is not None]
        out[(salt, solv)] = dict(gmix=st.median([x["dG_mix"] for x in a]),
                                 gmix_raw=st.median([x["dG_mix"] for x in r]),
                                 gsolv=st.median([x["dG_solv"] for x in a]),
                                 ggla=[x["c_ggla"] for x in a],
                                 n=len(a),
                                 sigma_ctl=(sum(sig) / len(sig)) if sig else None)
    return out
