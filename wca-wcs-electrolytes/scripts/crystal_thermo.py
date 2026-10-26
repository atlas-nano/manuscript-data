"""Assembles the crystal and gas-ion thermochemistry, the salt dissociation free energies and the ideal mixing entropies at 293.15 K from data/simulation/phonons and data/simulation/qchem_freq, writing data/tables/crystal_dissociation.csv, ideal_mixing.csv, crystal_thermo.csv and gas_ion_thermo.csv."""
import argparse
import csv
import gzip
import math
import os
import pathlib

import yaml

ROOT = pathlib.Path(__file__).resolve().parents[1]
R = 8.314462618
KCAL = 4184.0
T = 293.15
P = 1.0e5
KB = 1.380649e-23
H = 6.62607015e-34
AMU = 1.6605390666e-27
C = 2.99792458e10
C2 = H * C / KB
MASS = {"H": 1.008, "Li": 6.94, "B": 10.811, "C": 12.011, "N": 14.007, "O": 15.999,
        "F": 18.998403, "Al": 26.9815385, "P": 30.973762, "S": 32.06}
SIGMA = {"bf4": 12, "pf6": 24, "fsi": 2, "tfsi": 2, "tpfa": 1}
ANIONS = ("bf4", "pf6", "fsi", "tfsi", "tpfa")
SALT_OF = {"bf4": "libf4", "pf6": "lipf6", "fsi": "lifsi", "tfsi": "litfsi", "tpfa": "litpfa"}
SALT_NAME = {"libf4": "LiBF4", "lipf6": "LiPF6", "lifsi": "LiFSI", "litfsi": "LiTFSI", "litpfa": "LiTPFA"}
ANION_NAME = {"bf4": "BF4-", "pf6": "PF6-", "fsi": "FSI-", "tfsi": "TFSI-", "tpfa": "TPFA-"}
FU = {"libf4": 3, "lipf6": 3, "lifsi": 6, "litfsi": 4, "litpfa": 4}
ATOMS_PER_FU = {"libf4": 6, "lipf6": 8, "lifsi": 10, "litfsi": 16, "litpfa": 58}
SOLVENT_M = {"FM": (34.033, 0.5996), "ACN": (41.053, 0.786), "Me2O": (46.069, 0.670),
             "DME": (90.121, 0.868), "DEE": (74.123, 0.713)}
C_SALT = 0.25


def open_any(path):
    if os.path.exists(path):
        return open(path)
    return gzip.open(path + ".gz", "rt")


def ion_thermo(nu_list, geom, sigma):
    import numpy as np
    real = [x for x in nu_list if x > 0]
    zpe = 0.5 * sum(H * C * nu for nu in real) * 6.02214076e23 / KCAL
    evib = sum(R * C2 * nu / math.expm1(C2 * nu / T) for nu in real) / KCAL
    cvib = sum(R * (C2 * nu / T) ** 2 * math.exp(C2 * nu / T) / math.expm1(C2 * nu / T) ** 2 for nu in real) / KCAL
    svib = sum(R * ((C2 * nu / T) / math.expm1(C2 * nu / T) - math.log(-math.expm1(-C2 * nu / T))) for nu in real) / KCAL
    m = sum(MASS[e] for e, *_ in geom) * AMU
    q_t = (2 * math.pi * m * KB * T / H ** 2) ** 1.5 * (KB * T / P)
    s_tr = R * (math.log(q_t) + 2.5) / KCAL
    ms = [MASS[e] for e, *_ in geom]
    X = np.array([[x, y, z] for _, x, y, z in geom])
    Xc = X - np.average(X, axis=0, weights=ms)
    I = np.zeros((3, 3))
    for mi, r in zip(ms, Xc):
        I += mi * (np.dot(r, r) * np.eye(3) - np.outer(r, r))
    theta = [H ** 2 / (8 * math.pi ** 2 * max(v, 1e-12) * AMU * 1e-20 * KB) for v in np.linalg.eigvalsh(I)]
    q_r = (math.pi ** 0.5 / sigma) * (T ** 3 / (theta[0] * theta[1] * theta[2])) ** 0.5
    s_rot = R * (math.log(q_r) + 1.5) / KCAL
    return zpe, evib + 4.0 * R * T / KCAL, cvib + 4.0 * R / KCAL, s_tr + s_rot + svib, len(nu_list) - len(real)


def mesh_modes(path):
    modes, wq = [], None
    with open_any(path) as fh:
        for line in fh:
            t = line.strip()
            if t.startswith("weight:"):
                wq = int(t.split(":")[1])
            elif t.startswith("frequency:"):
                modes.append((float(t.split(":")[1]), wq))
    return modes


def harmonic(modes, nbands, temp):
    na = 6.02214076e23
    qw = sum(wq for _, wq in modes) / nbands
    zpe = eth = ent = cv = 0.0
    for nu, wq in modes:
        if nu <= 0:
            continue
        e = H * nu * 1e12
        zpe += wq * 0.5 * e
        x = e / (KB * temp)
        if x < 500:
            ex = math.exp(x)
            eth += wq * e / (ex - 1.0)
            ent += wq * (x / (ex - 1.0) - math.log(1.0 - math.exp(-x)))
            cv += wq * x * x * ex / (ex - 1.0) ** 2
    return zpe * na / 1000 / qw, eth * na / 1000 / qw, ent * KB * na / qw, cv * KB * na / qw


def soft_modes(modes, nbands):
    qw = sum(w for _, w in modes) / nbands
    dropped = sum(w for nu, w in modes if nu <= 0)

    def entropy(use_abs):
        ent = 0.0
        for nu, w in modes:
            v = abs(nu) if use_abs else nu
            if v <= 1e-6 or (nu <= 0 and not use_abs):
                continue
            x = H * v * 1e12 / (KB * T)
            if x > 500:
                continue
            ent += w * (x / (math.exp(x) - 1.0) - math.log(1.0 - math.exp(-x)))
        return ent * KB * 6.02214076e23 / qw

    s0, s1 = entropy(False), entropy(True)
    return 100.0 * dropped / sum(w for _, w in modes), 100 * (s1 - s0) / s0


def salt_thermo(d, salt):
    y = yaml.safe_load(open_any(os.path.join(d, "thermal_properties.yaml")))
    natom = int(y["natom"])
    nfu = round(natom / ATOMS_PER_FU[salt])
    modes = mesh_modes(os.path.join(d, "mesh.yaml"))
    grid = min(y["thermal_properties"], key=lambda r: abs(r["temperature"] - T))
    z_g, e_g, s_g, c_g = harmonic(modes, 3 * natom, grid["temperature"])
    f_g = z_g + e_g - grid["temperature"] * s_g / 1000
    for name, mine, ph in (("ZPE", z_g, float(y["zero_point_energy"])), ("F", f_g, grid["free_energy"]),
                           ("S", s_g, grid["entropy"]), ("Cv", c_g, grid["heat_capacity"])):
        if abs(mine - ph) > 1e-4 * max(abs(ph), 1.0):
            raise SystemExit(f"{salt}: mesh integration does not reproduce phonopy {name} at {grid['temperature']} K")
    zpe_k, eth_k, s_j, cv_j = harmonic(modes, 3 * natom, T)
    frac, dS_fold = soft_modes(modes, 3 * natom)
    return zpe_k / 4.184 / nfu, eth_k / 4.184 / nfu, cv_j / 4184.0 / nfu, s_j / 4184.0 / nfu, frac, dS_fold


def s_ideal_mix(c_solvent, c_salt=C_SALT):
    def s_mix(fracs, ntot):
        return -R * ntot * sum(x * math.log(x) for x in fracs if x > 0)
    n0 = c_salt + c_solvent
    n1 = 2 * c_salt + c_solvent
    return (s_mix([c_salt / n1, c_salt / n1, c_solvent / n1], n1) - s_mix([c_salt / n0, c_solvent / n0], n0)) / c_salt / KCAL


def read_qm(qdir, spectrum_set):
    freqs, geoms = {}, {}
    for row in csv.DictReader(open(os.path.join(qdir, "anion_frequencies.csv"))):
        freqs.setdefault((row["set"], row["anion"]), []).append(float(row["frequency_cm1"]))
    for row in csv.DictReader(open(os.path.join(qdir, "anion_geometries.csv"))):
        geoms.setdefault((row["set"], row["anion"]), []).append((row["element"], float(row["x"]), float(row["y"]), float(row["z"])))
    out = {}
    for a in ANIONS:
        key = (spectrum_set, a) if (spectrum_set, a) in freqs else ("assembly", a)
        out[a] = (freqs[key], geoms[key])
    return out


def write_csv(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--phonons", default=str(ROOT / "data" / "simulation" / "phonons"),
                    help="directory with <salt>/thermal_properties.yaml, mesh.yaml(.gz) and E_el.txt")
    ap.add_argument("--qchem", default=str(ROOT / "data" / "simulation" / "qchem_freq"))
    ap.add_argument("--anion-set", default="assembly", help="anion spectrum set in anion_frequencies.csv")
    ap.add_argument("--out", default=str(ROOT / "data" / "tables"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    q_t = (2 * math.pi * MASS["Li"] * AMU * KB * T / H ** 2) ** 1.5 * (KB * T / P)
    s_li = R * (math.log(q_t) + 2.5) / KCAL
    li_dh = 2.5 * R * T / KCAL
    gas = [["Li+", f"{0.0:.3f}", f"{li_dh:.3f}", f"{2.5 * R / KCAL:.5f}", f"{s_li:.4f}", 0]]
    ions = {}
    for an, (nu, g) in read_qm(a.qchem, a.anion_set).items():
        zpe, dh, cp, s, n_im = ion_thermo(nu, g, SIGMA[an])
        ions[an] = (zpe, dh, cp, s)
        gas.append([ANION_NAME[an], f"{zpe:.3f}", f"{dh:.3f}", f"{cp:.5f}", f"{s:.4f}", n_im])
    write_csv(os.path.join(a.out, "gas_ion_thermo.csv"),
              ["ion", "ZPE_kcal", "dH_0_to_T_kcal", "Cp_kcal_K", "S_kcal_K", "imaginary_modes_excluded"], gas)
    legs = {r["salt"]: r for r in csv.DictReader(open(os.path.join(a.qchem, "dissociation_legs.csv")))}
    rows, diss = [], []
    print(f"{'salt':8s} {'E_el/FU':>9s} {'ZPE':>8s} {'dH(0->T)':>9s} {'S':>7s} {'dH_diss':>8s} {'T*dS':>6s} {'dG_diss':>8s}")
    for an in ANIONS:
        salt = SALT_OF[an]
        d = os.path.join(a.phonons, salt)
        e_el = float(open(os.path.join(d, "E_el.txt")).read().strip()) / FU[salt]
        zpe_s, dh_s, cv_s, s_s, frac, dS_fold = salt_thermo(d, salt)
        L = legs[salt]
        d1, d2 = float(L["dE1_kcal"]), float(L["dE2_kcal"])
        dzpe = (ions[an][0] + 0.0) - zpe_s
        ddh = (ions[an][1] + li_dh) - dh_s
        dh_diss = d1 + d2 + dzpe + ddh
        ds = s_li + ions[an][3] - s_s
        dg = dh_diss - T * ds
        rows.append([SALT_NAME[salt], f"{e_el:.3f}", f"{zpe_s:.3f}", f"{dh_s:.3f}", f"{cv_s:.5f}", f"{s_s:.4f}",
                     f"{d1:.2f}", f"{d2:.2f}", f"{dzpe:.2f}", f"{ddh:.2f}", f"{dh_diss:.2f}",
                     f"{float(L['dE1_last_box_increment_kcal']):.2f}", f"{ds:.4f}", f"{T * ds:.2f}", f"{dg:.2f}", f"{frac:.2f}", f"{dS_fold:+.2f}"])
        diss.append([SALT_NAME[salt], f"{dg:.2f}"])
        print(f"{SALT_NAME[salt]:8s} {e_el:9.3f} {zpe_s:8.3f} {dh_s:9.3f} {s_s:7.4f} {dh_diss:8.2f} {T * ds:6.2f} {dg:8.2f}")
    write_csv(os.path.join(a.out, "crystal_thermo.csv"),
              ["salt", "E_el_per_FU_eV", "ZPE_kcal", "dH_0_to_T_kcal", "Cv_kcal_K", "S_kcal_K", "dE1_kcal", "dE2_kcal",
               "dZPE_kcal", "d_dH_kcal", "dH_diss_kcal", "dE1_box_increment_kcal", "dS_diss_kcal_K", "TdS_diss_kcal",
               "dG_diss_kcal", "imaginary_mode_fraction_pct", "dS_if_imaginary_folded_pct"], rows)
    write_csv(os.path.join(a.out, "crystal_dissociation.csv"), ["salt", "dG_diss"], diss)
    mix = []
    for name, (mw, rho) in SOLVENT_M.items():
        c = rho * 1000.0 / mw
        si = s_ideal_mix(c)
        mix.append([name, f"{mw:.3f}", f"{rho:.3f}", f"{c:.2f}", f"{si:.4f}", f"{T * si:.2f}"])
    write_csv(os.path.join(a.out, "ideal_mixing.csv"),
              ["solvent", "M_g_mol", "rho_g_cm3", "c_solvent_M", "dS_ideal_kcal_K", "TdS_ideal"], mix)
    print("ideal mixing (T*dS_ideal, kcal/mol): " + ", ".join(f"{r[0]} {r[5]}" for r in mix))
    print(f"-> {a.out}")


if __name__ == "__main__":
    main()
