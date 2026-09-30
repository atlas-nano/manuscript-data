#!/usr/bin/env python3
"""Route-2 (lr_ewald = 2, projected CG, reciprocal inside the matvec) gas-phase driver.

One isolated molecule, centred in a cubic periodic box of side L, pppm/samqeq with a pinned
g_ewald and an explicit mesh (the charges are zero at setup, so q2 = 0 and stock PPPM cannot size
the grid from the accuracy value: pppm.cpp:991 "Must use kspace_modify gewald for uncharged
system", and compute_df_kspace() = 0 breaks the grid loop at its coarsest step).

The external field is applied through fix efield's atom-style POTENTIAL variable, phi = -E*x, with an
atom-style field component: FixQEqBaseSam::init() (fix_qeq_base_sam.cpp:512-528) refuses a CONSTANT
field component along a periodic axis and admits only the atom-style-potential form, and
FixEfield::init() (fix_efield.cpp:235-236) sets varflag = ATOM only from an atom-style field
component, which is what the guard tests first. get_chi_field() (fix_qeq_base_sam.cpp:1428) then
takes chi_field[i] = qe2f*phi_i = -E*x_i (metal: qe2f = 1), identical to the constant-field branch
(:1416, chi_field = -(E.x)), and qeq_solve() adds it to the RHS as qb[i] -= chi_field[i]
(fix_qeq_sam_lr.cpp:2755 OpenMP / :2776 serial).

Usage: python3 run_r2.py <phase> [args]      (phases at the bottom)
Every run leaves in runs/<subdir>/{in,log,res}.<tag>; a run whose res file exists is not repeated.
"""
import concurrent.futures as cf
import json
import math
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
LMP = os.environ.get("LMP", "/home/tpascal/codes/lammps/lammps-dev2/build/lmp_parallel")
NP = int(os.environ.get("NP", "1"))            # ranks per run
WORKERS = int(os.environ.get("WORKERS", "2"))  # concurrent runs (cap: 8 cores in total)
QQRD2E = 14.399645                              # alpha[A^3] = QQRD2E * dmu[e A]/dE[V/A]

FIELDS7 = [-0.004, -0.002, -0.001, 0.0, 0.001, 0.002, 0.004]   # the published seven-point sweep

# ---------------------------------------------------------------------------------------- decks
ALK_DECK = """# route-2 alkane: cubic periodic box L={L} centred on the molecule, pppm/samqeq, field via potential
units           metal
atom_style      full
boundary        p p p
read_data       {datadir}/data.C{n}
variable        xc equal xcm(all,x)
variable        yc equal xcm(all,y)
variable        zc equal xcm(all,z)
variable        h  equal 0.5*{L}
change_box      all x final $(v_xc-v_h) $(v_xc+v_h) y final $(v_yc-v_h) $(v_yc+v_h) z final $(v_zc-v_h) $(v_zc+v_h) units box
variable        ext equal bound(all,xmax)-bound(all,xmin)
pair_style      hybrid/overlay lj/cut 12.0 coul/long 12.0 coul/shield/intra 12.0 all
pair_coeff      * * lj/cut 0.0 1.0
pair_coeff      * * coul/long
pair_coeff      * * coul/shield/intra 0.2499
kspace_style    pppm/samqeq 1.0e-6
kspace_modify   gewald {gewald} mesh {mesh} {mesh} {mesh} order {order}
neighbor        2.0 bin
fix             chg all qeq/sam 1 0.0 12.0 {tol} {param}
{solverline}
{selfline}
{modeline}
variable        E   equal {E}
variable        ex  atom v_E
variable        phi atom -1.0*v_E*x
fix             ef all efield v_ex 0.0 0.0 potential v_phi
compute         dip all dipole/samqeq
thermo_style    custom step pe c_dip[1]
run             0
print "RESULT_R2 n={n} L={L} E={E} mux=$(c_dip[1]:%.10f) muy=$(c_dip[2]:%.10f) pe=$(pe:%.8f) ecoul=$(ecoul:%.8f) elong=$(elong:%.8f) fchg=$(f_chg:%.8f) ext=$(v_ext:%.3f)"
"""

WAT_DECK = """# route-2 kernel comparison: TIP4P-FQ monomer, cubic periodic box L={L}, pppm/samqeq, kernel = {kernel}
units           metal
atom_style      full
boundary        p p p
variable        h equal 0.5*{L}
region          box block $(-v_h) $(v_h) $(-v_h) $(v_h) $(-v_h) $(v_h)
create_box      3 box
mass            1 15.9994
mass            2 1.008
mass            3 1.0e-6
molecule        h2o {datadir}/water4.mol
create_atoms    0 single 0.0 0.0 0.0 mol h2o 1
group           fq type 2 3
pair_style      hybrid/overlay lj/cut 12.0 coul/long 12.0 coul/shield/intra 12.0 all {kernelkw}
pair_coeff      * * lj/cut 0.0 1.0
pair_coeff      1 1 lj/cut 0.0124 3.159
pair_coeff      * * coul/long
pair_coeff      1 1 coul/shield/intra 3.0803
pair_coeff      2 2 coul/shield/intra 1.7007
pair_coeff      3 3 coul/shield/intra 3.0803
kspace_style    pppm/samqeq 1.0e-6
kspace_modify   gewald {gewald} mesh {mesh} {mesh} {mesh} order {order}
neighbor        2.0 bin
fix             chg fq qeq/sam 1 0.0 12.0 {tol} {param}
fix_modify      chg shield {kernelkw}
{solverline}
{selfline}
compute         dip all dipole/samqeq
run             0
variable        qM  equal q[4]
variable        qH1 equal q[2]
variable        qH2 equal q[3]
variable        mu  equal c_dip
print "RESULT_R2 kernel={kernel} L={L} qM=$(v_qM:%.6f) qH1=$(v_qH1:%.6f) qH2=$(v_qH2:%.6f) mu_D=$(v_mu:%.6f) pe=$(pe:%.6f) ecoul=$(ecoul:%.6f) elong=$(elong:%.6f) fchg=$(f_chg:%.6f) evdwl=$(evdwl:%.6f)"
write_dump      all custom charges_{kernel}_L{L}.dump id type q modify sort id format line "%d %d %.12g"
"""

SOLVER = "fix_modify      chg solver cg nofallback"
SELF = {"peratom": "fix_modify      chg recip_self peratom",
        "probes16": "fix_modify      chg recip_probes 16",
        "default": ""}


def run(subdir, tag, deck, np=None):
    """Run one deck; return the RESULT_R2 line (cached in res.<tag>)."""
    np = np or NP
    d = os.path.join(HERE, "runs", subdir)
    os.makedirs(d, exist_ok=True)
    res = os.path.join(d, f"res.{tag}")
    if os.path.exists(res):
        return open(res).read().strip()
    inp = os.path.join(d, f"in.{tag}")
    with open(inp, "w") as fh:
        fh.write(deck)
    r = subprocess.run(["mpirun", "-np", str(np), LMP, "-in", inp, "-log", f"log.{tag}"],
                       cwd=d, stdin=subprocess.DEVNULL, capture_output=True, text=True, timeout=7200)
    m = re.search(r"^RESULT_R2 .*$", r.stdout, re.M)
    if not m:
        tail = "\n  ".join((r.stdout + r.stderr).strip().splitlines()[-8:])
        raise RuntimeError(f"{subdir}/{tag}: no result\n  {tail}")
    # solver banner: route + iterations, for the record
    ban = re.search(r"^samqeq SOLVER-DIAG.*$", r.stdout, re.M)
    st = os.stat(LMP)   # the binary was rebuilt by another session mid-campaign (19:29): stamp every result
    line = m.group(0) + (f"  |  {ban.group(0)}" if ban else "") + f"  |  bin={int(st.st_mtime)}:{st.st_size}"
    with open(res, "w") as fh:
        fh.write(line + "\n")
    return line


def parse(line):
    out = {}
    for k, v in re.findall(r"(\w+)=(\S+)", line.split("|")[0]):
        try:
            out[k] = float(v)
        except ValueError:
            out[k] = v
    return out


def alk(n, L, E, mesh, subdir, tag, gewald=0.35, order=5, tol="1.0e-10", self="peratom",
        modeline="", param="alkane_r2.param", solver=True, np=None):
    deck = ALK_DECK.format(n=n, L=L, E=E, mesh=mesh, gewald=gewald, order=order, tol=tol,
                           param=os.path.join(HERE, "alkane", param),
                           datadir=os.path.join(HERE, "alkane"),
                           solverline=SOLVER if solver else "", selfline=SELF[self], modeline=modeline)
    return parse(run(subdir, tag, deck, np))


def wat(kernel, L, mesh, subdir, tag, gewald=0.35, order=5, tol="1.0e-10", self="peratom",
        param="water4_r2.param", np=None):
    kw = {"slater": "slater 2s 1 3", "cbrt": "cbrt"}[kernel]
    deck = WAT_DECK.format(kernel=kernel, kernelkw=kw, L=L, mesh=mesh, gewald=gewald, order=order,
                           tol=tol, param=os.path.join(HERE, "kernel", param),
                           datadir=os.path.join(HERE, "kernel"),
                           solverline=SOLVER, selfline=SELF[self])
    return parse(run(subdir, tag, deck, np))


def pmap(fn, jobs):
    """Run jobs (list of kwargs dicts) with at most WORKERS concurrent LAMMPS processes."""
    with cf.ThreadPoolExecutor(max_workers=WORKERS) as ex:
        futs = [ex.submit(fn, **j) for j in jobs]
        return [f.result() for f in futs]


def fit(xs, ys):
    k = len(xs)
    sx, sy = sum(xs), sum(ys)
    sxx = sum(x * x for x in xs)
    sxy = sum(x * y for x, y in zip(xs, ys))
    slope = (k * sxy - sx * sy) / (k * sxx - sx * sx)
    icpt = (sy - slope * sx) / k
    ybar = sy / k
    ss_res = sum((y - (slope * x + icpt)) ** 2 for x, y in zip(xs, ys))
    ss_tot = sum((y - ybar) ** 2 for y in ys)
    return slope, icpt, (1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan"))


def sweep(n, L, mesh, subdir, fields=FIELDS7, **kw):
    """Seven-point (or given) field sweep -> alpha, alpha/n, mu0, R^2."""
    tagbase = kw.pop("tagbase", f"C{n}_L{L}_m{mesh}")
    jobs = [dict(n=n, L=L, E=E, mesh=mesh, subdir=subdir, tag=f"{tagbase}_E{E}", **kw) for E in fields]
    res = pmap(alk, jobs)
    mus = [r["mux"] for r in res]
    slope, _, r2 = fit(fields, mus)
    mu0 = mus[fields.index(0.0)] if 0.0 in fields else float("nan")
    return {"n": n, "L": L, "mesh": mesh, "alpha": QQRD2E * slope, "alpha_per_n": QQRD2E * slope / n,
            "mu0": mu0, "r2": r2, "mus": dict(zip(fields, mus)), "pe": {f: r["pe"] for f, r in zip(fields, res)}}


def loglog(ns, alphas):
    return fit([math.log(n) for n in ns], [math.log(a) for a in alphas])[0]


# ---------------------------------------------------------------------------------------- phases
def phase_smoke():
    print(wat("slater", 40, 80, "smoke", "wat_slater_L40"))
    print(wat("cbrt", 40, 80, "smoke", "wat_cbrt_L40"))
    print(alk(8, 60, 0.004, 120, "smoke", "C8_L60_E0.004"))
    print(alk(8, 60, 0.004, 120, "smoke", "C8_L60_E0.004_bs", modeline="fix_modify      chg bondsoft 10.0 3.0"))


def _show(label, r):
    print(f"{label:<34} alpha={r['alpha']:.4f}  alpha/n={r['alpha_per_n']:.4f}  mu0={r['mu0']:.2e}  R2={r['r2']:.9f}")


def phase_protocol():
    """Box, mesh, g_ewald, self-term and header-inertness checks (alkane C24 / C8)."""
    print("== C24 box convergence, mesh spacing 0.5 A, alpha from ±0.004 (2-pt) with mu(0) ==")
    for L in (40, 60, 80, 100):
        _show(f"C24 L={L} mesh={2*L}", sweep(24, L, 2 * L, "protocol", fields=[-0.004, 0.0, 0.004],
                                            tagbase=f"C24_L{L}_m{2*L}"))
    print("== C24 mesh convergence at L=80 ==")
    for mesh, order in ((80, 5), (120, 5), (160, 5), (200, 5), (120, 7)):
        _show(f"C24 L=80 mesh={mesh} order={order}", sweep(24, 80, mesh, "protocol", fields=[-0.004, 0.004],
                                                          tagbase=f"C24_L80_m{mesh}_o{order}", order=order))
    print("== C8 g_ewald independence (rc = 12 A), L=60 mesh=120 ==")
    for g in (0.30, 0.35, 0.40):
        _show(f"C8 gewald={g}", sweep(8, 60, 120, "protocol", fields=[-0.004, 0.004],
                                      tagbase=f"C8_L60_m120_g{g}", gewald=g))
    print("== C8 self-term choice, E=+0.004 ==")
    for s in ("peratom", "probes16", "default"):
        r = alk(8, 60, 0.004, 120, "protocol", f"C8_L60_m120_self_{s}", self=s)
        print(f"  self={s:<9} mux={r['mux']:.10f} pe={r['pe']:.8f}")
    print("== C8 header fields 1-4 inert on route 2, E=+0.004 ==")
    for p in ("alkane_r2.param", "hdr_a.param", "hdr_b.param", "hdr_c.param"):
        r = alk(8, 60, 0.004, 120, "protocol", f"C8_L60_m120_{p}", param=p)
        print(f"  {p:<16} mux={r['mux']:.10f} pe={r['pe']:.8f}")


def phase_efield():
    """C8: mu(0), antisymmetry, step-independence of the finite-difference alpha."""
    r = sweep(8, 60, 120, "efield", tagbase="C8_L60_m120")
    _show("C8 7-point", r)
    m = r["mus"]
    for E in (0.001, 0.002, 0.004):
        a = QQRD2E * (m[E] - m[-E]) / (2 * E)
        print(f"  central diff ±{E}: alpha={a:.4f}   mu(+E)+mu(-E)={m[E]+m[-E]:.2e}")
    print(f"  mu(0)={m[0.0]:.3e} e A")
    # pe(E) is even in E to the order of the linear response: E(E) - E(0) = -alpha E^2/2 (in e A V/A = eV)
    for E in (0.002, 0.004):
        print(f"  pe(+E)-pe(0)={r['pe'][E]-r['pe'][0.0]:.6e}  pe(-E)-pe(0)={r['pe'][-E]-r['pe'][0.0]:.6e}  "
              f"-alpha E^2/(2 QQRD2E)={-r['alpha']*E*E/(2*QQRD2E):.6e}")


def phase_water():
    """Kernel comparison: box, mesh, self-term and g_ewald convergence for both kernels; writes water_route2.json.
    The box series is extrapolated in L^-3 (the Lorentz term is the leading one; the molecule is 1 A)."""
    out = {}
    for k in ("slater", "cbrt"):
        box = {}
        for L in (40, 60, 80, 100, 120):
            mesh = mesh_of(L, H_WAT)
            r = wat(k, L, mesh, "water", f"{k}_L{L}_m{mesh}")
            box[L] = r
            print(f"{k:<6} L={L:<3} mesh={mesh:<3} qM={r['qM']:.6f} qH={r['qH1']:.6f} mu={r['mu_D']:.6f} pe={r['pe']:.6f} "
                  f"ecoul={r['ecoul']:.6f} elong={r['elong']:.6f} fchg={r['fchg']:.6f}")
        # two-box L^-3 extrapolation from the two largest boxes, checked against the third
        def x3(key, La, Lb):
            f = Lb ** -3 / (La ** -3 - Lb ** -3)
            return box[Lb][key] + (box[Lb][key] - box[La][key]) * f
        ext = {key: x3(key, 100, 120) for key in ("qM", "qH1", "mu_D", "pe", "ecoul", "elong", "fchg")}
        chk = {key: x3(key, 80, 100) for key in ("qM", "mu_D", "pe")}
        print(f"{k:<6} L->inf (100,120): qM={ext['qM']:.6f} qH={ext['qH1']:.6f} mu={ext['mu_D']:.6f} pe={ext['pe']:.6f}   "
              f"[from (80,100): qM={chk['qM']:.6f} mu={chk['mu_D']:.6f} pe={chk['pe']:.6f}]")
        out[k] = {"box": {L: {kk: vv for kk, vv in r.items() if kk != "kernel"} for L, r in box.items()},
                  "inf": ext, "inf_check": chk}
        for mesh in (40, 80, 120, 160):
            r = wat(k, 40, mesh, "water", f"{k}_L40_m{mesh}")
            print(f"{k:<6} L=40  mesh={mesh:<3} qM={r['qM']:.6f} mu={r['mu_D']:.6f} pe={r['pe']:.6f}")
        for s in ("probes16", "default"):
            try:
                r = wat(k, 60, 120, "water", f"{k}_L60_m120_self_{s}", self=s)
                print(f"{k:<6} L=60  self={s:<8} qM={r['qM']:.6f} mu={r['mu_D']:.6f} pe={r['pe']:.6f}")
                out[k][f"self_{s}"] = {kk: vv for kk, vv in r.items() if kk != "kernel"}
            except RuntimeError as e:
                msg = [l for l in str(e).splitlines() if "ERROR" in l]
                print(f"{k:<6} L=60  self={s:<8} REFUSED: {msg[0].strip() if msg else str(e)[:200]}")
                out[k][f"self_{s}"] = "refused: " + (msg[0].strip() if msg else "")
        for g in (0.30, 0.40):
            r = wat(k, 60, 120, "water", f"{k}_L60_m120_g{g}", gewald=g)
            print(f"{k:<6} L=60  gewald={g} qM={r['qM']:.6f} mu={r['mu_D']:.6f} pe={r['pe']:.6f}")
            out[k][f"gewald_{g}"] = {kk: vv for kk, vv in r.items() if kk != "kernel"}
    with open(os.path.join(HERE, "water_route2.json"), "w") as fh:
        json.dump(out, fh, indent=1)


def phase_boxcheck(order="5"):
    """Third-box consistency of the Lorentz + L^-5 protocol (C24 plain at L=160; C8 at 60/120/200) and the
    box-independence of the bondsoft saddle (C24 bondsoft at L=200 vs 120)."""
    order = int(order)
    out = {}
    ac = {}
    for L in (120, 160, 200):
        mesh = mesh_of(L, H_ALK)
        r = sweep(24, L, mesh, "series_plain", tagbase=f"C24_L{L}_m{mesh}_o{order}", order=order)
        ac[L] = lorentz(r["alpha"], L)
        print(f"C24 plain L={L} box alpha={r['alpha']:.4f}  Lorentz-corrected={ac[L]:.4f}")
    for La, Lb in ((120, 160), (160, 200), (120, 200)):
        print(f"  extrap5({La},{Lb}) = {extrap5(ac[La], ac[Lb], La, Lb):.4f}")
    out["C24_plain_corrected"] = ac
    a8 = {}
    for L in (60, 120, 200):
        mesh = mesh_of(L, H_ALK)
        r = sweep(8, L, mesh, "series_plain" if L > 60 else "protocol", tagbase=f"C8_L{L}_m{mesh}_o{order}", order=order)
        a8[L] = (r["alpha"], lorentz(r["alpha"], L))
        print(f"C8 plain L={L} box alpha={r['alpha']:.5f}  Lorentz-corrected={a8[L][1]:.5f}")
    out["C8_plain"] = a8
    bs = {}
    for L in (120, 200):
        mesh = mesh_of(L, H_ALK)
        r = sweep(24, L, mesh, "series_bondsoft", tagbase=f"C24_L{L}_m{mesh}_o{order}", order=order, modeline=BS)
        bs[L] = r["alpha"]
        print(f"C24 bondsoft L={L} alpha={r['alpha']:.6f}")
    out["C24_bondsoft"] = bs
    with open(os.path.join(HERE, "boxcheck_route2.json"), "w") as fh:
        json.dump(out, fh, indent=1)


NS = [4, 6, 8, 10, 12, 14, 17, 20, 24]
BS = "fix_modify      chg bondsoft 10.0 3.0"


def lorentz(a_box, L):
    """Remove the tinfoil image term at dipole order. Under conducting (tinfoil) boundary conditions a cubic
    lattice of identical dipoles feels the Lorentz cavity field 4 pi mu/(3 V), so the box polarizability is
    alpha_box = alpha/(1 - 4 pi alpha/(3 L^3)), i.e. 1/alpha = 1/alpha_box + 4 pi/(3 L^3). Exact for the
    dipole term; what remains is the higher-multipole coupling of an extended molecule, O(L^-5)."""
    return 1.0 / (1.0 / a_box + 4.0 * math.pi / (3.0 * L ** 3))


def extrap5(a1, a2, L1, L2):
    """Two-box Richardson extrapolation of the Lorentz-corrected values assuming a(L) = a_inf + c L^-5."""
    f = L2 ** -5 / (L1 ** -5 - L2 ** -5)
    return a2 + (a2 - a1) * f


def series_fit(rows):
    """rows: list of (n, alpha). Returns full-series exponent, long-chain-half exponent, 4.5-fold ratio."""
    rows = sorted(rows)
    ns = [r[0] for r in rows]
    al = [r[1] for r in rows]
    p_all = loglog(ns, al)
    h = len(rows) // 2
    p_hi = loglog(ns[h:], al[h:])
    apn = [a / n for n, a in rows]
    return p_all, p_hi, ns[h], max(apn) / min(apn), apn[-1] / apn[0]


H_ALK = 2.0 / 3.0     # mesh spacing for the alkanes (C8 at L=60: h=1.0 -> 5e-6 rel. off, h=0.67 -> 1.6e-7)
H_WAT = 0.5           # mesh spacing for the water monomer (h=1.0 -> 5e-5 rel. off in q_M, h=0.5 -> converged)


def mesh_of(L, h):
    return int(round(L / h))


def phase_series(L1="120", L2="200", order="5"):
    """Plain (atomic-charge) series at two boxes, Lorentz-corrected and L^-5-extrapolated; bondsoft series
    at L1 (no reciprocal term in that saddle). Writes alpha_series_route2.csv and series_route2.json."""
    L1, L2, order = int(L1), int(L2), int(order)
    import csv
    res = {}
    for n in NS:
        for L in (L1, L2):
            mesh = mesh_of(L, H_ALK)
            res[("plain", n, L)] = sweep(n, L, mesh, "series_plain", tagbase=f"C{n}_L{L}_m{mesh}_o{order}", order=order)
        mesh = mesh_of(L1, H_ALK)
        res[("bondsoft", n, L1)] = sweep(n, L1, mesh, "series_bondsoft", tagbase=f"C{n}_L{L1}_m{mesh}_o{order}",
                                         order=order, modeline=BS)
    rows = []
    for n in NS:
        r1, r2 = res[("plain", n, L1)], res[("plain", n, L2)]
        a1c, a2c = lorentz(r1["alpha"], L1), lorentz(r2["alpha"], L2)
        a_inf = extrap5(a1c, a2c, L1, L2)
        rows.append({"mode": "plain", "n": n, "npts": 7, "alpha": a_inf, "alpha_per_n": a_inf / n,
                     "mu0": r2["mu0"], "r2": r2["r2"], "alpha_L1": r1["alpha"], "alpha_L2": r2["alpha"],
                     "alpha_L1c": a1c, "alpha_L2c": a2c, "mu0_L1": r1["mu0"], "r2_L1": r1["r2"], "L1": L1, "L2": L2})
    for n in NS:
        b = res[("bondsoft", n, L1)]
        rows.append({"mode": "bondsoft", "n": n, "npts": 7, "alpha": b["alpha"], "alpha_per_n": b["alpha_per_n"],
                     "mu0": b["mu0"], "r2": b["r2"], "alpha_L1": b["alpha"], "alpha_L2": "", "alpha_L1c": "",
                     "alpha_L2c": "", "mu0_L1": b["mu0"], "r2_L1": b["r2"], "L1": L1, "L2": ""})
    with open(os.path.join(HERE, "alpha_series_route2.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    summary = {}
    for mode in ("plain", "bondsoft"):
        sub = [r for r in rows if r["mode"] == mode]
        for key in (("alpha_L1", "alpha_L2", "alpha_L1c", "alpha_L2c", "alpha") if mode == "plain" else ("alpha",)):
            p_all, p_hi, nh, ratio, ratio_ends = series_fit([(r["n"], r[key]) for r in sub])
            a8 = next(r[key] for r in sub if r["n"] == 8)
            a17 = next(r[key] for r in sub if r["n"] == 17)
            p_loc = math.log(a17 / a8) / math.log(17 / 8)
            flat = 100 * (max(r[key] / r["n"] for r in sub) / min(r[key] / r["n"] for r in sub) - 1)
            summary[f"{mode}:{key}"] = dict(p_all=p_all, p_hi=p_hi, n_hi=nh, apn_ratio=ratio, apn_ratio_ends=ratio_ends,
                                            a8=a8, a17=a17, p_local=p_loc, flat_pct=flat,
                                            min_r2=min(r["r2"] for r in sub), max_mu0=max(abs(r["mu0"]) for r in sub))
            print(f"{mode:>8} {key:<8}: n^{p_all:.3f} (full)  n^{p_hi:.3f} (n>={nh})  alpha/n ratio max/min {ratio:.2f} "
                  f"(C24/C4 {ratio_ends:.2f})  a8={a8:.2f} a17={a17:.2f} p(8->17)={p_loc:.3f}  flat={flat:.1f}%  "
                  f"minR2={summary[f'{mode}:{key}']['min_r2']:.8f} max|mu0|={summary[f'{mode}:{key}']['max_mu0']:.1e}")
    with open(os.path.join(HERE, "series_route2.json"), "w") as fh:
        json.dump({"summary": summary, "rows": rows}, fh, indent=1)
    for r in rows:
        print(f"{r['mode']:>8} n={r['n']:>2} alpha_L1={r['alpha_L1']!s:>12} alpha_L2={r['alpha_L2']!s:>12} "
              f"alpha_inf={r['alpha']:.4f} alpha/n={r['alpha_per_n']:.4f} mu0={r['mu0']:.1e} R2={r['r2']:.9f}")


def phase_kappa(L1="120", L2="200", order="5"):
    """SI scan rows on route 2: the plain (±0.002 central difference) reference p and alpha(C8), the bold
    bondsoft row (3.0, 10), and the kappa -> infinity limit (1e4, 1e5, 1e6 at bcut 3.0)."""
    L1, L2, order = int(L1), int(L2), int(order)
    E = 0.002
    out = {}

    def cd(n, L, modeline, tagbase, sub):
        mesh = mesh_of(L, H_ALK)
        rp = alk(n, L, +E, mesh, sub, f"{tagbase}_L{L}_m{mesh}_E{E}", order=order, modeline=modeline)
        rm = alk(n, L, -E, mesh, sub, f"{tagbase}_L{L}_m{mesh}_E{-E}", order=order, modeline=modeline)
        return rp["mux"], rm["mux"]

    # plain reference: box value -> Lorentz-corrected -> L^-5-extrapolated
    a = {}
    for n in (8, 17):
        p1, m1 = cd(n, L1, "", f"plain_C{n}", "scan")
        p2, m2 = cd(n, L2, "", f"plain_C{n}", "scan")
        a[(n, L1)] = QQRD2E * (p1 - m1) / (2 * E)
        a[(n, L2)] = QQRD2E * (p2 - m2) / (2 * E)
        a[n] = extrap5(lorentz(a[(n, L1)], L1), lorentz(a[(n, L2)], L2), L1, L2)
    p_loc = math.log(a[17] / a[8]) / math.log(17 / 8)
    out["plain"] = dict(a8=a[8], a17=a[17], p=p_loc, a8_L1=a[(8, L1)], a8_L2=a[(8, L2)], a17_L1=a[(17, L1)], a17_L2=a[(17, L2)],
                        a8_L1c=lorentz(a[(8, L1)], L1), a8_L2c=lorentz(a[(8, L2)], L2),
                        a17_L1c=lorentz(a[(17, L1)], L1), a17_L2c=lorentz(a[(17, L2)], L2))
    print(f"plain ±0.002 reference (corrected+extrapolated): alpha(C8)={a[8]:.4f} alpha(C17)={a[17]:.4f} p={p_loc:.3f}   "
          f"[box L1: {a[(8, L1)]:.4f}/{a[(17, L1)]:.4f}  box L2: {a[(8, L2)]:.4f}/{a[(17, L2)]:.4f}  "
          f"corrected L1: {out['plain']['a8_L1c']:.4f}/{out['plain']['a17_L1c']:.4f}  L2: {out['plain']['a8_L2c']:.4f}/{out['plain']['a17_L2c']:.4f}]")
    # bondsoft rows (no reciprocal in the bondsoft saddle: one box suffices)
    for kappa in (10.0, 1e4, 1e5, 1e6):
        ml = f"fix_modify      chg bondsoft {kappa:g} 3.0"
        a8 = QQRD2E * (lambda p, m: p - m)(*cd(8, L1, ml, f"bs_k{kappa:g}_C8", "scan")) / (2 * E)
        a17 = QQRD2E * (lambda p, m: p - m)(*cd(17, L1, ml, f"bs_k{kappa:g}_C17", "scan")) / (2 * E)
        pk = math.log(a17 / a8) / math.log(17 / 8)
        out[f"bondsoft_k{kappa:g}"] = dict(a8=a8, a17=a17, p=pk)
        print(f"bondsoft bcut=3.0 kappa={kappa:<8g} alpha(C8)={a8:.4f} alpha(C17)={a17:.4f} p={pk:.3f}")
    with open(os.path.join(HERE, "scan_route2.json"), "w") as fh:
        json.dump(out, fh, indent=1)


if __name__ == "__main__":
    globals()["phase_" + sys.argv[1]](*sys.argv[2:])
