"""Computes the per-box xPT salt solvation and mixing free energies behind Figures 3F and S3 from data/simulation/xpt and writes data/simulation/xpt_dgmix_293K.csv and xpt_dgmix_anion_clamp_293K.csv."""
import argparse
import csv
import gzip
import json
import math
import pathlib

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
XPT = ROOT / "data" / "simulation" / "xpt"
KB, PLANCK, EV_J, NA = 1.380649e-23, 6.62607015e-34, 1.602176634e-19, 6.02214076e23
EV2JMOL = 96485.33212
EV_KJ = 96.48533
EV_KCAL = 23.060548
KJ_KCAL = 1 / 4.184
R_GAS = 8.314462618
KCAL_J = 4184.0
GROUPS = ("ion", "anion", "bulk", "r226", "r260", "r300", "a260", "a300", "a350")
ISQ, IN, IDOF, IV, IEMD = range(5)
COLUMNS = ("system", "k", "c_ggla", "dH", "dS_Li", "dS_an", "dG_solv", "dG_mix", "band_lo", "band_hi", "branch", "sigma_ctl")


def open_text(path):
    p = pathlib.Path(path)
    if p.exists():
        return open(p, newline="")
    return gzip.open(str(p) + ".gz", "rt", newline="")


def read_rows(name):
    with open_text(XPT / name) as f:
        return list(csv.DictReader(f))


def load_inputs():
    with open_text(XPT / "constants.json") as f:
        const = json.load(f)
    with open_text(XPT / "shell_control.json") as f:
        control = json.load(f)
    boxes = read_rows("boxes.csv")
    windows = {}
    for r in read_rows("group_windows.csv"):
        windows.setdefault((r["label"], r["group"]), []).append(
            [float(r["S_q"]), float(r["N"]), float(r["DOF"]), float(r["V"]), float(r["E_md"])])
    windows = {k: np.array(v) for k, v in windows.items()}
    neat, neat_rotsym = {}, {}
    for r in read_rows("neat_windows.csv"):
        neat.setdefault(r["label"], []).append([float(r["V"]), float(r["Nmol"]), float(r["E_md"])])
        neat_rotsym[r["label"]] = r["rotsym"]
    neat = {k: np.array(v) for k, v in neat.items()}
    gg = {}
    for r in read_rows("group_group.csv"):
        gg.setdefault(r["label"], []).append([float(r["step"]), float(r["gg_ls"]), float(r["gg_as"]), float(r["gg_la"])])
    gg = {k: np.array(v) for k, v in gg.items()}
    vacuum = {}
    for r in read_rows("vacuum_windows.csv"):
        vacuum.setdefault(r["anion"], []).append((int(r["window"]), r))
    return const, control, boxes, windows, neat, neat_rotsym, gg, vacuum


def s_st(m_amu, v_a3, t):
    m = m_amu / NA / 1000.0
    lam = PLANCK / np.sqrt(2 * np.pi * m * KB * t)
    return KB * (np.log(v_a3 * 1e-30 / lam ** 3) + 2.5) / EV_J


def vacuum_reference(rows, anion, const):
    srot, svib, cages = [], [], []
    for _, r in sorted(rows, key=lambda x: x[0]):
        if int(float(r["Nmol"])) != 1:
            raise ValueError(f"{anion}: vacuum window with Nmol != 1")
        cages.append(abs(float(r["S_cage_r"])))
        srot.append(float(r["S_rot"]))
        svib.append(float(r["S_vib"]))
    if len(srot) < const["min_windows_vacuum"]:
        raise ValueError(f"{anion}: {len(srot)} vacuum windows")
    if cages and max(cages) > 0.0:
        raise ValueError(f"{anion}: vacuum windows carry a cage term")
    s_rot, s_vib = float(np.mean(srot)), float(np.mean(svib))
    target = const["rigid_rotor_S_rot_J_per_mol_K"][anion]
    if abs(s_rot * EV2JMOL - target) / target > const["rigid_rotor_tolerance"]:
        raise ValueError(f"{anion}: S_rot {s_rot * EV2JMOL:.1f} J/mol/K against a rigid rotor {target}")
    return s_rot, s_vib


def s_ideal_mix(c_solvent, c_salt):
    def s_mix(fracs, ntot):
        return -R_GAS * ntot * sum(x * math.log(x) for x in fracs if x > 0)
    n0 = c_salt + c_solvent
    s0 = s_mix([c_salt / n0, c_solvent / n0], n0)
    n1 = 2 * c_salt + c_solvent
    s1 = s_mix([c_salt / n1, c_salt / n1, c_solvent / n1], n1)
    return (s1 - s0) / c_salt / KCAL_J


def t_ds_ideal(solvent, const):
    key = const["solvent_key"][solvent]
    mw = const["solvent_molar_mass_g_per_mol"][key]
    rho = const["solvent_density_g_per_cm3"][key]
    return const["temperature_K"] * s_ideal_mix(rho * 1000.0 / mw, const["c_salt_mol_per_L"])


def solvation_enthalpy(box, windows, neat, neat_rotsym, gg, const):
    label = box["label"]
    s = windows[(label, "bulk")]
    n = neat[label]
    if len(s) < const["min_windows"] or len(n) < const["min_windows"]:
        raise ValueError(f"{label}: windows {len(s)} / {len(n)}")
    vs, vn = s[:, IV].mean(), n[:, 0].mean()
    ns, nn = s[:, IN].mean(), n[:, 1].mean()
    dev = abs(vn - vs) / vs
    if dev > const["volume_tolerance"]:
        raise ValueError(f"{label}: neat volume differs by {100 * dev:.2f} %")
    if abs(nn - ns) > 0.5:
        raise ValueError(f"{label}: neat has {nn:.0f} molecules, solution {ns:.0f}")
    if box["rotsym"] != neat_rotsym[label]:
        raise ValueError(f"{label}: xPT rotational symmetry differs between solution and neat runs")
    g = gg[label]
    ggls, ggas = g[:, 1].mean() * EV_KJ, g[:, 2].mean() * EV_KJ
    reorg = ns * (s[:, IEMD].mean() - n[:, 2].mean()) * EV_KJ
    inter = 0.5 * (ggls + ggas)
    return inter + reorg


def per_box(box, windows, gg, vac, const):
    label = box["label"]
    n_ion = int(box["n_ion"])
    t = const["temperature_K"]
    a_solv = int(box["atoms_per_solvent"])
    io, bk, an = (windows[(label, g)][1:] for g in ("ion", "bulk", "anion"))
    v = io[:, IV].mean()
    ds_ion = (io[:, ISQ].mean() - s_st(const["li_molar_mass_amu"], v, t)) * EV2JMOL
    sq_an = an[:, ISQ].mean()
    if np.isnan(sq_an):
        return None
    anion = const["anion_of_salt"][box["salt"]]
    s_rot, s_vib = vac[anion]
    s_gas = s_st(const["anion_molar_mass_amu"][anion], an[:, IV].mean(), t) + s_rot + s_vib
    ds_an_gas = (sq_an - s_gas) * EV2JMOL
    g = gg[label][1:]
    ggls, ggas, ggla = g[:, 1:4].mean(axis=0)
    s_bulk = bk[:, ISQ].mean()
    ds_li = {}
    for tag in const["li_shell_radii"]:
        sh = windows.get((label, tag))
        if sh is None:
            continue
        sh = sh[1:]
        if len(sh) < const["min_windows"]:
            continue
        nmol = sh[:, IN].mean()
        atoms = sh[:, IDOF].mean() / 3.0
        if not abs(atoms / nmol - a_solv) <= const["composition_tolerance"]:
            continue
        ds_li[tag] = ds_ion + (nmol / n_ion) * (sh[:, ISQ].mean() - s_bulk) * EV2JMOL
    shell_n, shell_term, ds_an = {}, {}, {}
    for tag in const["anion_shell_radii"]:
        sh = windows.get((label, tag))
        if sh is None:
            continue
        sh = sh[1:]
        if len(sh) < const["min_windows"]:
            continue
        nmol_a = sh[:, IN].mean()
        shell_n[tag] = float(nmol_a)
        shell_term[tag] = (nmol_a / n_ion) * (sh[:, ISQ].mean() - s_bulk) * EV2JMOL
        ds_an[tag] = ds_an_gas + shell_term[tag]
    if not ds_li or not ds_an:
        return None
    half = 0.5 * (g[:, 1] + g[:, 2]) * EV_KCAL / n_ion
    sem = float(np.std(half, ddof=1)) / math.sqrt(len(half))
    return dict(ggla=ggla * EV_KCAL, ds_li=ds_li, ds_an=ds_an, ds_an_gas=ds_an_gas, shell_term=shell_term,
                shell_n=shell_n, n_ion=n_ion, sem=sem)


def headline_pair(li_tags, an_tags, const):
    pairs = [(lr, ar) for lr in li_tags for ar in an_tags]
    preferred = [p for p in pairs if list(p) == const["headline_pair"]]
    return preferred[0] if preferred else pairs[len(pairs) // 2]


def interpolate(xs, ys, x):
    if x <= xs[0]:
        return ys[0]
    if x >= xs[-1]:
        return ys[-1]
    for i in range(1, len(xs)):
        if x <= xs[i]:
            frac = (x - xs[i - 1]) / (xs[i] - xs[i - 1])
            return ys[i - 1] + frac * (ys[i] - ys[i - 1])
    return ys[-1]


def control_at(curve, n):
    xs = [float(p["N"]) for p in curve]
    ys = [float(p["dS_per_mol_mean"]) for p in curve]
    sem = [float(p["sem_windows"]) if p.get("sem_windows") is not None else 0.0 for p in curve]
    sd = [p.get("sd_seeds") for p in curve]
    unc = [math.sqrt(s * s + float(d) * float(d)) for s, d in zip(sem, sd)] if all(d is not None for d in sd) else sem
    n = float(n)
    if n < xs[0] or n > xs[-1]:
        edge = 0 if n < xs[0] else -1
        return ys[edge], unc[edge]
    return interpolate(xs, ys, n), interpolate(xs, unc, n)


def apply_anion_control(terms, curve, const):
    t = const["temperature_K"]
    n_ion = terms["n_ion"]
    li_tags, an_tags = list(terms["ds_li"]), list(terms["shell_term"])
    lr, ar = headline_pair(li_tags, an_tags, const)
    unc_terms = {}
    for tag in an_tags:
        value, unc = control_at(curve, terms["shell_n"][tag])
        ctl_term = (terms["shell_n"][tag] / n_ion) * value
        unc_terms[tag] = (terms["shell_n"][tag] / n_ion) * (unc or 0.0)
        corrected = terms["shell_term"][tag] - ctl_term
        terms["shell_term"][tag] = corrected
        terms["ds_an"][tag] = terms["ds_an_gas"] + corrected
    uncs = [unc_terms[ar]]
    terms["sigma_ctl"] = t * math.sqrt(sum(u * u for u in uncs)) / KCAL_J
    return terms


def assemble_row(box, terms, dh, const):
    t = const["temperature_K"]
    dgd = const["dg_diss_kcal_per_mol"][box["salt"]]
    tsi = t_ds_ideal(box["solvent"], const)
    vals = []
    for lr, sli in terms["ds_li"].items():
        for ar, san in terms["ds_an"].items():
            ds = sli + san
            dg_solv = dh * KJ_KCAL - t * ds * KJ_KCAL / 1000.0
            vals.append((lr, ar, ds, dg_solv, dg_solv + dgd - tsi))
    preferred = [v for v in vals if v[0] == const["headline_pair"][0] and v[1] == const["headline_pair"][1]]
    lr, ar, ds, dg_solv, dg_mix = (preferred or [vals[len(vals) // 2]])[0]
    sem = terms["sem"]
    lo = min(min(v[4] for v in vals), dg_mix - 2 * sem)
    hi = max(max(v[4] for v in vals), dg_mix + 2 * sem)
    branch = "free" if terms["ggla"] > const["contact_threshold_kcal_per_mol"] else "paired"
    sigma = terms.get("sigma_ctl")
    return [box["system"], box["k"], f"{terms['ggla']:.1f}", f"{dh:.1f}", f"{terms['ds_li'][lr]:.1f}",
            f"{terms['ds_an'][ar]:.1f}", f"{dg_solv:.1f}", f"{dg_mix:.2f}", f"{lo:.2f}", f"{hi:.2f}", branch,
            "" if sigma is None else f"{sigma:.2f}"]


def compute(corrected):
    const, control, boxes, windows, neat, neat_rotsym, gg, vacuum = load_inputs()
    vac = {an: vacuum_reference(rows, an, const) for an, rows in vacuum.items()}
    out = []
    for box in boxes:
        terms = per_box(box, windows, gg, vac, const)
        if terms is None:
            continue
        try:
            dh = solvation_enthalpy(box, windows, neat, neat_rotsym, gg, const)
        except ValueError:
            continue
        if corrected:
            terms = apply_anion_control(terms, control[box["solvent"]], const)
        out.append(assemble_row(box, terms, dh, const))
    return out


def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(COLUMNS)
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", default=str(ROOT / "data" / "simulation"))
    a = ap.parse_args()
    out = pathlib.Path(a.out_dir)
    raw = compute(corrected=False)
    clamp = compute(corrected=True)
    write_csv(out / "xpt_dgmix_293K.csv", raw)
    write_csv(out / "xpt_dgmix_anion_clamp_293K.csv", clamp)
    print(f"{len(raw)} boxes -> {out / 'xpt_dgmix_293K.csv'}")
    print(f"{len(clamp)} boxes -> {out / 'xpt_dgmix_anion_clamp_293K.csv'}")


if __name__ == "__main__":
    main()
