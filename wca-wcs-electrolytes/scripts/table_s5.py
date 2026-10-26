"""Computes the ideal mixing entropy of a 0.25 M salt in each solvent at 293.15 K (SI equation S3) from the solvent molar masses and experimental liquid densities and writes data/tables/table_s5.csv."""
import argparse
import csv
import math
import os
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
R = 8.314462618
KCAL = 4184.0
T = 293.15
C_SALT = 0.25
SOLVENTS = [("FM", 34.033, 0.5996), ("ACN", 41.053, 0.786), ("Me2O", 46.069, 0.670), ("DME", 90.121, 0.868),
            ("DEE", 74.123, 0.713)]


def s_mix(fracs, ntot):
    return -R * ntot * sum(x * math.log(x) for x in fracs if x > 0)


def ds_ideal(c_solv, c_salt=C_SALT):
    n0, n1 = c_salt + c_solv, 2 * c_salt + c_solv
    return (s_mix([c_salt / n1, c_salt / n1, c_solv / n1], n1) - s_mix([c_salt / n0, c_solv / n0], n0)) / c_salt / KCAL


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=str(ROOT / "data" / "tables" / "table_s5.csv"))
    a = ap.parse_args()
    rows = []
    for name, mw, rho in SOLVENTS:
        c = rho * 1000.0 / mw
        s = ds_ideal(c)
        rows.append([name, f"{mw:.3f}", f"{rho:.3f}", f"{c:.2f}", f"{s:.4f}", f"{T * s:.2f}"])
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Solvent", "M (g mol-1)", "rho (g cm-3)", "c(solv) (M)", "dS_ideal (kcal mol-1 K-1)",
                    "TdS_ideal at 293.15 K (kcal mol-1)"])
        w.writerows(rows)
    for r in rows:
        print(" ".join(r))
    print(f"-> {a.out}")


if __name__ == "__main__":
    main()
