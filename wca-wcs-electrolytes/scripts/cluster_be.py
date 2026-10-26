"""Computes the Li+-solvent cluster binding energies BE(n) per solvent molecule, with the ligand-ligand strain removed, from the Q-Chem cluster and ligand-only single points and writes data/tables/cluster_binding.csv and cluster_binding_table.csv."""
import argparse
import csv
import os
import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parents[1]
HA_KCAL = 627.5094740631
INFO = {"FM": ("CH3F", "Li-F, monodentate", 1.800), "Me2O": ("(CH3)2O", "Li-O, bisector", 1.750),
        "DEE": ("(C2H5)2O", "Li-O, bisector", 1.800), "ACN": ("CH3CN", "Li-N, monodentate", 1.950),
        "DME": ("CH3OCH2CH2OCH3", "Li-O, chelate (2 O)", 2.168)}


def energies(path):
    return [float(m) for m in re.findall(r"Total energy in the final basis set =\s*(-?\d+\.\d+)", open(path).read())]


def cluster_block(path):
    e = energies(path)
    if len(e) < 3:
        raise SystemExit(f"{path}: expected at least monomer, Li+ and n = 1")
    return e[0], e[1], {n: e[1 + n] for n in range(1, len(e) - 1)}


def ligand_only(path):
    txt = open(path).read()
    tags = re.findall(r"ligand-only, (\w+)_n(\d)", txt)
    vals = energies(path)
    return {(t[0], int(t[1])): v for t, v in zip(tags, vals)}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--solvent", action="append", required=True, metavar="NAME=CLUSTER:LIGAND_ONLY:TAG:NMAX",
                    help="cluster output (jobs: monomer, Li+, Li(L)n for n = 1, 2, ...), ligand-only output whose "
                         "$comment lines read 'ligand-only, <TAG>_n<n>', and the largest n reported")
    ap.add_argument("--out", default=str(ROOT / "data" / "tables"))
    a = ap.parse_args()
    rows, order, li_refs = {}, [], {}
    for spec in a.solvent:
        name, rest = spec.split("=", 1)
        clus_f, lig_f, tag, nmax = rest.rsplit(":", 3)
        nmax = int(nmax)
        order.append(name)
        lig = ligand_only(lig_f)
        mono, eli, clus = cluster_block(clus_f)
        li_refs[name] = eli
        prev = None
        for n in (1, 2, 3, 4):
            if n == 1:
                e_lig = mono
            elif n not in clus or n > nmax or (tag, n) not in lig:
                prev = None
                continue
            else:
                e_lig = lig[(tag, n)]
            tot = (clus[n] - eli - e_lig) * HA_KCAL
            rows[(name, n)] = dict(raw=(clus[n] - eli - n * mono) * HA_KCAL, tot=tot, strain=(e_lig - n * mono) * HA_KCAL,
                                   seq=tot - prev if prev is not None else tot, per=tot / n)
            prev = tot
    if len(set(f"{v:.10f}" for v in li_refs.values())) != 1:
        raise SystemExit("the Li+ reference energy differs between the cluster series")
    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, "cluster_binding.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["solvent", "n", "BE_raw", "BE_frozen", "strain", "BE_seq", "per_molecule"])
        for name in order:
            for n in (1, 2, 3, 4):
                d = rows.get((name, n))
                if d:
                    w.writerow([name, n] + [f"{d[k]:.2f}" for k in ("raw", "tot", "strain", "seq", "per")])
    with open(os.path.join(a.out, "cluster_binding_table.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Solvent", "Formula", "Coordination", "r (A)", "n = 1", "n = 2", "n = 3", "n = 4", "Strain, n = 4"])
        for name in sorted(order, key=lambda s: -rows[(s, 1)]["per"]):
            f, c, r = INFO.get(name, ("", "", float("nan")))
            cells = [f"{rows[(name, n)]['per']:.1f}" if (name, n) in rows else "" for n in (1, 2, 3, 4)]
            w.writerow([name, f, c, f"{r:.2f}"] + cells + [f"{rows[(name, 4)]['strain']:.1f}" if (name, 4) in rows else ""])
    for name in order:
        print(f"{name:5s} " + "  ".join(f"n={n}: {rows[(name, n)]['per']:7.2f}" for n in (1, 2, 3, 4) if (name, n) in rows))
    print(f"-> {a.out}/cluster_binding.csv, cluster_binding_table.csv")


if __name__ == "__main__":
    main()
