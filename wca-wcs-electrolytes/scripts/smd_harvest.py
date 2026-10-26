"""Reads chained Q-Chem SMD outputs (job 1 gas phase, later jobs the same geometry in each solvent) and writes the continuum solvation free energies of Li+ and the anions to data/tables/smd_solvation.csv and smd_solvation_table.csv."""
import argparse
import csv
import os
import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parents[1]
HA_KCAL = 627.509474
E_TOT = re.compile(r"(?:Total energy in the final basis set\s*=|^ Total energy =)\s*(-?\d+\.\d+)", re.M)
G_ENP = re.compile(r"(?:G-ENP\(liq\).*?|G_ENP\s*=\s*)(-?\d+\.\d+)\s*a\.u\.")
G_CDS = re.compile(r"G-CDS\(liq\)|G_CDS\s*=")
KCAL = re.compile(r"(-?\d+\.\d+)\s*kcal/mol")
SPECIES = [("li", "Li+"), ("bf4", "BF4-"), ("pf6", "PF6-"), ("fsi", "FSI-"), ("tfsi", "TFSI-"), ("tpfa", "TPFA-")]
SOLVENTS = [("acn", "ACN"), ("dme", "DME"), ("dee", "DEE"), ("me2o", "Me2O"), ("fm", "FM")]


def parse(path):
    txt = open(path, errors="replace").read()
    parts = re.split(r"Running Job (\d+) of (\d+)", txt)
    jobs = []
    for i in range(1, len(parts), 3):
        body = parts[i + 2]
        e = E_TOT.findall(body)
        enp = G_ENP.findall(body)
        cds = None
        m = G_CDS.search(body)
        if m:
            k = KCAL.search(body, m.end())
            if k:
                cds = float(k.group(1))
        jobs.append(dict(n=int(parts[i]), e_tot=float(e[-1]) if e else None,
                         g_enp=float(enp[-1]) if enp else None, g_cds=cds))
    return jobs


def solvation(path, solvents, tol):
    jobs = {j["n"]: j for j in parse(path)}
    e_gas = jobs[1]["e_tot"]
    if jobs[1]["g_enp"] is not None:
        raise SystemExit(f"{path}: job 1 is not a gas-phase reference")
    out = {}
    for idx, sv in enumerate(solvents):
        j = jobs.get(idx + 2)
        if sv == "-" or j is None or j["e_tot"] is None or j["g_enp"] is None or j["g_cds"] is None:
            continue
        dgA = (j["e_tot"] - e_gas) * HA_KCAL
        dgB = (j["g_enp"] - e_gas) * HA_KCAL + j["g_cds"]
        if abs(dgA - dgB) > tol * HA_KCAL:
            raise SystemExit(f"{path} {sv}: total-energy and decomposition routes disagree ({dgA:.4f} vs {dgB:.4f})")
        out[sv] = dict(dg=dgA, cds=j["g_cds"], enp=(j["g_enp"] - e_gas) * HA_KCAL)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run", action="append", required=True, metavar="ION=FILE:SOLV1,SOLV2,...",
                    help="one output file, the ion it holds and the solvent of jobs 2, 3, ... in order ('-' skips a job); "
                         "a later --run never overwrites a value an earlier one supplied")
    ap.add_argument("--tol", type=float, default=1e-5, help="agreement required between the two routes (Ha)")
    ap.add_argument("--out", default=str(ROOT / "data" / "tables"))
    a = ap.parse_args()
    rows = {}
    for spec in a.run:
        ion, rest = spec.split("=", 1)
        path, _, sl = rest.rpartition(":")
        for sv, v in solvation(path, sl.split(","), a.tol).items():
            rows.setdefault(ion, {}).setdefault(sv, v)
    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, "smd_solvation.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["ion", "solvent", "dG_solv_kcal", "electrostatic_kcal", "cavity_dispersion_kcal"])
        for k, lab in SPECIES:
            for sk, sl in SOLVENTS:
                d = rows.get(k, {}).get(sk)
                if d is not None:
                    w.writerow([lab, sl, repr(d["dg"]), repr(d["enp"]), repr(d["cds"])])
    with open(os.path.join(a.out, "smd_solvation_table.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Ion"] + [sl for _, sl in SOLVENTS])
        for k, lab in SPECIES:
            w.writerow([lab] + [f"{rows[k][sk]['dg']:.2f}" if sk in rows.get(k, {}) else "" for sk, _ in SOLVENTS])
    for k, lab in SPECIES:
        print(f"{lab:6s} " + " ".join(f"{sl}:{rows[k][sk]['dg']:8.2f}" for sk, sl in SOLVENTS if sk in rows.get(k, {})))
    print(f"-> {a.out}/smd_solvation.csv, smd_solvation_table.csv")


if __name__ == "__main__":
    main()
