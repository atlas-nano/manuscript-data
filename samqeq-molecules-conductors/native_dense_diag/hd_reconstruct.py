#!/usr/bin/env python3
"""Reconstruct the constrained samQEq charge operator of ONE frame as fix qeq/sam applies it on the
lr_ewald=2 route, validate the reconstruction against the charges LAMMPS solved on that frame, and
diagonalize the constrained operator densely.

OPERATOR (fix qeq/sam, lr_ewald=2, lr_alpha=0 header -> lr_alpha := kspace g_ewald at run time)
  variables      atoms of the fix group; solve variable qs = q - q0[type]          fix_qeq_sam_lr.cpp qeq_solve
  real space     H_ij = qqrd2e [ J(r) - erf(a r)/r ],  r <= swb, group-group        fix_qeq_sam.cpp compute_H
                 J = 1/cbrt(r^3 + 1/g^3), g = sqrt(gamma_i gamma_j)      (shield cbrt, default)
                 J = erf(alpha_ij r)/r, alpha_i = lam/(2 Rc_i^2),
                     alpha_ij = sqrt(alpha_i alpha_j/(alpha_i+alpha_j))  (shield pqeq)   fix_qeq_sam.cpp shielded_coulomb
  reciprocal     out_i += qqrd2e [ sum_j K_ij x_j - rs x_i ]  (PPPM k!=0 potential incl. j=i; rs =
                 recip_self, calibrated or pinned)                                    fix_qeq_sam_lr.cpp add_reciprocal
                 -> here K_ij := exact Ewald k-space sum psi_k(r_ij) (converged)
  diagonal       eta[type] (+ (c4/6) q_lag^2 on the quartic group, Picard-lagged)     fix_qeq_sam_levels.cpp apply_quartic_eta
  constraint     P = per-molecule mean subtraction over group atoms                   fix_qeq_sam_lr.cpp project_neutral
  RHS            b = -(chi + field_fixed + q0field);  P b                            fix_qeq_sam_lr.cpp qeq_solve
                 field_fixed: shielded real-space + reciprocal field of the NON-group fixed charges,
                 own bonded shell = spring (its J skipped, its reciprocal cancelled)  fix_qeq_sam_lr.cpp add_fixed_charge_field
                 q0field = P[(H_real + qqrd2e (K - rs I)) q0]                         fix_qeq_sam_lr.cpp coulomb_field
The solve is  P (D + H_real + qqrd2e (K - rs I)) P qs = P b;  q = qs + q0.
"""
import argparse, json, sys, time
import numpy as np
from math import erf, pi, sqrt

QQRD2E = 14.399645          # LAMMPS metal units, eV A / e^2


# ----------------------------------------------------------------------------- inputs
def read_dump(path):
    with open(path) as fh:
        lines = fh.readlines()
    i = 0
    box = []
    while not lines[i].startswith("ITEM: BOX BOUNDS"):
        i += 1
    for k in range(3):
        lo, hi = (float(v) for v in lines[i + 1 + k].split()[:2])
        box.append((lo, hi))
    i += 4
    hdr = lines[i].split()[2:]
    col = {k: j for j, k in enumerate(hdr)}
    rows = [l.split() for l in lines[i + 1:] if l.strip()]
    ids = np.array([int(r[col["id"]]) for r in rows])
    mol = np.array([int(r[col["mol"]]) for r in rows])
    typ = np.array([int(r[col["type"]]) for r in rows])
    q = np.array([float(r[col["q"]]) for r in rows])
    x = np.array([[float(r[col["x"]]), float(r[col["y"]]), float(r[col["z"]])] for r in rows])
    o = np.argsort(ids)
    return ids[o], mol[o], typ[o], q[o], x[o], np.array([hi - lo for lo, hi in box])


def read_data(path):
    """Atoms (full style) -> {id: (mol, type, q)}; Bonds -> list of (a1, a2)."""
    atoms, bonds = {}, []
    sec = None
    with open(path) as fh:
        for l in fh:
            s = l.split()
            if not s:
                continue
            if s[0] in ("Atoms", "Bonds", "Velocities", "Angles", "Masses", "Pair", "Bond", "Angle"):
                sec = s[0]
                continue
            if sec == "Atoms" and s[0].lstrip("-").isdigit():
                atoms[int(s[0])] = (int(s[1]), int(s[2]), float(s[3]))
            elif sec == "Bonds" and s[0].isdigit():
                bonds.append((int(s[2]), int(s[3])))
    return atoms, bonds


def read_param(path):
    with open(path) as fh:
        lines = [l.split() for l in fh if l.strip() and not l.startswith("#")]
    hdr = [float(v) for v in lines[0]]
    header = dict(gamma_align=hdr[0], kappa_bond=hdr[1], r_ov=hdr[2],
                  r_loc=hdr[3] if len(hdr) > 3 else 2.5,
                  lr_alpha=hdr[4] if len(hdr) > 4 else 0.0,
                  lr_ewald=int(hdr[5] + 0.5) if len(hdr) > 5 else 0,
                  lr_ridge=hdr[6] if len(hdr) > 6 else 0.0)
    per = {}
    for f in lines[1:]:
        per[int(f[0])] = dict(chi=float(f[1]), eta=float(f[2]), col4=float(f[3]), q0=float(f[4]))
    return header, per


# ----------------------------------------------------------------------------- kernels
def kernel_pair_params(kind, per, lam):
    """Return f(ti, tj) -> (J(r) callable on arrays, J(0))."""
    def make(ti, tj):
        if kind == "cbrt":
            g = sqrt(per[ti]["col4"] * per[tj]["col4"])
            return (lambda r: 1.0 / np.cbrt(r ** 3 + 1.0 / g ** 3)), g
        elif kind == "pqeq":
            ai = lam * 0.5 / per[ti]["col4"] ** 2
            aj = lam * 0.5 / per[tj]["col4"] ** 2
            aij = sqrt(ai * aj / (ai + aj))
            from scipy.special import erf as verf
            return (lambda r: verf(aij * r) / r), 2.0 * aij / sqrt(pi)
        raise ValueError(kind)
    return make


def min_image(x, box):
    d = x[:, None, :] - x[None, :, :]
    d -= box * np.round(d / box)
    return np.sqrt((d ** 2).sum(-1))


def ewald_k_matrix(x, box, a, kexp):
    """psi_k(r_ij) = (4 pi / V) sum_{k != 0} exp(-k^2/4a^2)/k^2 cos(k.r_ij), all i, j (incl. i=j),
    k on the reciprocal lattice of the orthogonal box, spherical cutoff k^2/(4a^2) <= kexp."""
    V = float(np.prod(box))
    nmax = np.ceil(np.sqrt(kexp) * 2 * a * box / (2 * pi)).astype(int)
    rng = [np.arange(-n, n + 1) for n in nmax]
    N = np.array(np.meshgrid(*rng, indexing="ij")).reshape(3, -1).T
    # half space: first nonzero component positive
    keep = (N[:, 0] > 0) | ((N[:, 0] == 0) & (N[:, 1] > 0)) | ((N[:, 0] == 0) & (N[:, 1] == 0) & (N[:, 2] > 0))
    N = N[keep]
    K = 2 * pi * N / box
    k2 = (K ** 2).sum(1)
    sel = k2 / (4 * a * a) <= kexp
    K, k2 = K[sel], k2[sel]
    w = (4 * pi / V) * np.exp(-k2 / (4 * a * a)) / k2
    ph = x @ K.T                                   # (N, Nk)
    C = np.cos(ph) * np.sqrt(w)
    S = np.sin(ph) * np.sqrt(w)
    psi = 2.0 * (C @ C.T + S @ S.T)
    return psi, len(k2)


def neutral_basis(mols):
    """Orthonormal block basis B (n x (n - nfrag)) of the per-fragment-neutral subspace."""
    n = len(mols)
    um = np.unique(mols)
    B = np.zeros((n, n - len(um)))
    c = 0
    for m in um:
        idx = np.where(mols == m)[0]
        k = len(idx)
        A = np.eye(k) - np.ones((k, k)) / k
        w, V = np.linalg.eigh(A)
        Vk = V[:, w > 0.5]
        B[np.ix_(idx, np.arange(c, c + k - 1))] = Vk
        c += k - 1
    return B


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", required=True)
    ap.add_argument("--param", required=True)
    ap.add_argument("--data", help="LAMMPS data file: starting charges (quartic Picard) and bonds (own shell)")
    ap.add_argument("--kernel", choices=["cbrt", "pqeq"], required=True)
    ap.add_argument("--lam", type=float, default=0.462770)
    ap.add_argument("--gewald", type=float, required=True, help="kspace g_ewald of the run (1/A)")
    ap.add_argument("--rs", type=float, required=True, help="recip_self actually used (raw, 1/A)")
    ap.add_argument("--swb", type=float, required=True)
    ap.add_argument("--solve-types", type=int, nargs="+", required=True)
    ap.add_argument("--c4", type=float, default=0.0)
    ap.add_argument("--c4-types", type=int, nargs="*", default=[])
    ap.add_argument("--mix", type=float, default=0.5)
    ap.add_argument("--niter", type=int, default=8)
    ap.add_argument("--scf-tol", type=float, default=1.0e-5)
    ap.add_argument("--kexp", type=float, default=40.0)
    ap.add_argument("--out", required=True, help="JSON summary path")
    ap.add_argument("--label", default="")
    args = ap.parse_args()

    t0 = time.time()
    ids, mol, typ, qL, x, box = read_dump(args.dump)
    header, per = read_param(args.param)
    if header["lr_ewald"] != 2:
        sys.exit(f"this script reconstructs the lr_ewald=2 route; header says {header['lr_ewald']}")
    a, rs, swb = args.gewald, args.rs, args.swb
    V = float(np.prod(box))
    natoms = len(ids)

    var = np.isin(typ, args.solve_types)
    src = (~var) & (qL != 0.0)
    iv, isrc = np.where(var)[0], np.where(src)[0]
    nv = len(iv)
    mols_v = mol[iv]
    nfrag = len(np.unique(mols_v))
    eta = np.array([per[t]["eta"] for t in typ[iv]])
    chi = np.array([per[t]["chi"] for t in typ[iv]])
    q0 = np.array([per[t]["q0"] for t in typ[iv]])

    # own bonded shell of each variable (Drude spring, not a Coulomb pair)
    own = {}
    q_start = None
    if args.data:
        atoms, bonds = read_data(args.data)
        idpos = {int(i): k for k, i in enumerate(ids)}
        for a1, a2 in bonds:
            k1, k2 = idpos[a1], idpos[a2]
            if var[k1] and src[k2]:
                own[k1] = k2
            elif var[k2] and src[k1]:
                own[k2] = k1
        q_start = np.array([atoms[int(i)][2] for i in ids])
    if len(isrc) and not own:
        print("NOTE: fixed sources present but no core-shell bonds found (no --data?) -> no own-shell rule")

    # geometry
    r = min_image(x, box)
    make = kernel_pair_params(args.kernel, per, args.lam)
    types = sorted(set(typ.tolist()))
    Jfun, J0 = {}, {}
    for ti in types:
        for tj in types:
            Jfun[(ti, tj)], J0[(ti, tj)] = make(ti, tj)

    # real-space matrix over ALL atoms (kernel - erf(ar)/r), r <= swb, r > 0
    Hreal = np.zeros((natoms, natoms))
    Jker = np.zeros((natoms, natoms))          # bare kernel J(r) within swb (kernel-only diagnostic)
    with np.errstate(divide="ignore", invalid="ignore"):
        from scipy.special import erf as verf
        for ti in types:
            for tj in types:
                m = (typ[:, None] == ti) & (typ[None, :] == tj) & (r <= swb) & (r > 0)
                if not m.any():
                    continue
                rr = r[m]
                Jv = Jfun[(ti, tj)](rr)
                Hreal[m] = QQRD2E * (Jv - verf(a * rr) / rr)
                Jker[m] = QQRD2E * Jv
    # reciprocal (exact Ewald k-space, includes i=j)
    psi, nk = ewald_k_matrix(x, box, a, args.kexp)
    phi0 = psi[0, 0]
    # convergence check of the k-sum: redo with a smaller cutoff on the diagonal only
    psi_chk, nk_chk = ewald_k_matrix(x[:1], box, a, args.kexp * 0.6)
    two_a_sqrtpi = 2 * a / sqrt(pi)
    xi_img = phi0 - two_a_sqrtpi - pi / (a * a * V)      # Ewald image self-term (should be ~ -2.837/L for cubic)

    # fixed-charge field on the variables
    field = np.zeros(nv)
    if len(isrc):
        qs_src = qL[isrc]
        Hr_vs = Hreal[np.ix_(iv, isrc)].copy()
        psi_vs = psi[np.ix_(iv, isrc)]
        for kv, ks in own.items():
            jv = np.where(iv == kv)[0][0]
            js = np.where(isrc == ks)[0][0]
            rr = r[kv, ks]
            Hr_vs[jv, js] = -QQRD2E * erf(a * rr) / rr if rr <= swb else 0.0   # own shell: cancel reciprocal only
        field = Hr_vs @ qs_src + QQRD2E * (psi_vs @ qs_src)

    # variable-variable operator pieces
    Hr_vv = Hreal[np.ix_(iv, iv)]
    K_vv = QQRD2E * psi[np.ix_(iv, iv)]
    Joff = Hr_vv + K_vv - QQRD2E * rs * np.eye(nv)        # J_offdiag as coulomb_field applies it (no eta)
    B = neutral_basis(mols_v)
    P = B @ B.T

    q0field = P @ (Joff @ q0) if np.any(q0 != 0) else np.zeros(nv)
    b = -(chi + field + q0field)
    Bb = B.T @ b

    def solve(eta_diag):
        H = Joff + np.diag(eta_diag)
        A = B.T @ H @ B
        y = np.linalg.solve(A, Bb)
        return B @ y + q0

    qgrp = np.isin(typ[iv], args.c4_types) if args.c4 else np.zeros(nv, bool)

    def eta_eff(qlag):
        e = eta.copy()
        if args.c4:
            e[qgrp] += args.c4 / 6.0 * qlag[qgrp] ** 2
        return e

    res = dict(label=args.label, dump=args.dump, param=args.param, kernel=args.kernel, gewald=a, rs=rs,
               swb=swb, box=box.tolist(), V=V, natoms=natoms, nvar=nv, nsrc=int(len(isrc)), nfrag=nfrag,
               nk_halfspace=nk, phi0=phi0, two_a_sqrtpi=two_a_sqrtpi, xi_image=xi_img,
               xi_cubic_estimate=-2.837297 / box[0], phi0_check_lower_kexp=psi_chk[0, 0],
               uniform_shift_eV=QQRD2E * (two_a_sqrtpi - rs), header=header,
               per_type={t: per[t] for t in types}, own_shells=len(own))

    # ---- validation against the LAMMPS charges
    qLv = qL[iv]
    molsum = np.array([qLv[mols_v == m].sum() for m in np.unique(mols_v)])
    res["lammps_molsum_max_dev"] = float(np.abs(molsum - q0[mols_v == np.unique(mols_v)[0]].sum()).max())
    if args.c4:
        if q_start is None:
            sys.exit("--c4 needs --data for the Picard starting charges")
        qlag = q_start[iv].copy()
        hist = []
        for it in range(args.niter):
            qn = solve(eta_eff(qlag))
            dq = float(np.abs(qn - qlag).max())
            qc = args.mix * qn + (1 - args.mix) * qlag
            hist.append(dq)
            if dq < args.scf_tol:
                break
            qlag = qc
        qrec = qc
        res["picard_iters"] = len(hist)
        res["picard_dq_hist"] = hist
        # exact fixed point (undamped, tight)
        qf = qrec.copy()
        for it in range(200):
            qn = solve(eta_eff(qf))
            d = float(np.abs(qn - qf).max())
            qf = qn
            if d < 1e-13:
                break
        res["fixedpoint_iters"] = it + 1
        res["fixedpoint_vs_lammps_max"] = float(np.abs(qf - qLv).max())
        res["fixedpoint_vs_committed_max"] = float(np.abs(qf - qrec).max())
        eta_code = eta_eff(qLv)
    else:
        qrec = solve(eta)
        eta_code = eta
    d = qrec - qLv
    res["val_max_abs_dq"] = float(np.abs(d).max())
    res["val_rms_dq"] = float(np.sqrt((d ** 2).mean()))
    res["val_max_abs_dq_by_type"] = {int(t): float(np.abs(d[typ[iv] == t]).max()) for t in args.solve_types}
    res["q_lammps_mean_by_type"] = {int(t): float(qLv[typ[iv] == t].mean()) for t in args.solve_types}
    res["q_lammps_min_max_by_type"] = {int(t): [float(qLv[typ[iv] == t].min()), float(qLv[typ[iv] == t].max())]
                                       for t in args.solve_types}

    # ---- spectra on the constrained subspace
    def spectrum(H, name):
        A = B.T @ H @ B
        A = 0.5 * (A + A.T)
        w = np.linalg.eigvalsh(A)
        # null count of P H P in the full space
        wfull = np.linalg.eigvalsh(0.5 * ((P @ H @ P) + (P @ H @ P).T))
        nnull = int((np.abs(wfull) < 1e-8 * np.abs(wfull).max()).sum())
        out = dict(lam_min=float(w[0]), lam_2=float(w[1]), lam_max=float(w[-1]),
                   kappa=float(w[-1] / w[0]) if w[0] > 0 else None,
                   n_neg=int((w < 0).sum()), dim=int(len(w)), n_null_full=nnull)
        print(f"  {name:38s} lam_min {w[0]:10.5f}  lam_max {w[-1]:10.5f}  kappa {out['kappa']}  "
              f"neg {out['n_neg']}  null {nnull}")
        return out

    Ipsi0 = QQRD2E * (phi0 - rs)
    print("spectra (eV):")
    spec = {}
    spec["code_operator"] = spectrum(Joff + np.diag(eta_code), "code operator (as solved)")
    if args.c4:
        spec["code_operator_bare_eta"] = spectrum(Joff + np.diag(eta), "code operator, bare eta (no quartic)")
    spec["rs_equals_2a_sqrtpi"] = spectrum(Hr_vv + K_vv - QQRD2E * two_a_sqrtpi * np.eye(nv) + np.diag(eta),
                                            "rs := 2a/sqrt(pi) (no pin offset)")
    spec["kernel_only_no_ewald"] = spectrum(Jker[np.ix_(iv, iv)] + np.diag(eta), "kernel only (no reciprocal), bare eta")
    res["spectra"] = spec

    # ---- contact values, closest approaches
    st = args.solve_types
    res["J0_eV"] = {f"{ti}-{tj}": QQRD2E * J0[(ti, tj)] for ti in st for tj in st if ti <= tj}
    res["eta_min"] = float(eta.min())
    rvv = r[np.ix_(iv, iv)]
    inter = mols_v[:, None] != mols_v[None, :]
    intra = (~inter) & (rvv > 0)
    close = {}
    for ti in st:
        for tj in st:
            if ti > tj:
                continue
            m = (typ[iv][:, None] == ti) & (typ[iv][None, :] == tj)
            rin = rvv[m & inter].min() if (m & inter).any() else None
            ria = rvv[m & intra].min() if (m & intra).any() else None
            close[f"{ti}-{tj}"] = dict(r_inter_min=float(rin) if rin else None,
                                      J_at_r_inter_min=float(QQRD2E * Jfun[(ti, tj)](np.array([rin]))[0]) if rin else None,
                                      r_intra_min=float(ria) if ria else None,
                                      J_at_r_intra_min=float(QQRD2E * Jfun[(ti, tj)](np.array([ria]))[0]) if ria else None)
    res["closest"] = close
    res["seconds"] = time.time() - t0

    print(f"\nframe {args.dump}: {natoms} atoms, {nv} variables in {nfrag} fragments, {len(isrc)} fixed sources, "
          f"{len(own)} own shells; box {box[0]:.6f}")
    print(f"g_ewald {a}  rs {rs}  2a/sqrt(pi) {two_a_sqrtpi:.6f}  phi0 {phi0:.6f} (kexp check {psi_chk[0,0]:.6f})  "
          f"xi_image {xi_img:.6f} (cubic est {-2.837297/box[0]:.6f})  uniform shift {QQRD2E*(two_a_sqrtpi-rs):+.5f} eV")
    print(f"VALIDATION vs LAMMPS: max|dq| {res['val_max_abs_dq']:.3e} e   rms {res['val_rms_dq']:.3e} e   "
          f"by type {res['val_max_abs_dq_by_type']}   (LAMMPS per-mol sum dev {res['lammps_molsum_max_dev']:.2e})")
    if args.c4:
        print(f"  Picard replicated: {res['picard_iters']} iters, dq hist {['%.2e' % h for h in hist]}; "
              f"fixed point vs LAMMPS {res['fixedpoint_vs_lammps_max']:.2e}, vs committed {res['fixedpoint_vs_committed_max']:.2e}")
    print(f"J(0) eV: {res['J0_eV']}   eta_min {res['eta_min']}")
    print(f"closest: {json.dumps(close)}")
    with open(args.out, "w") as fh:
        json.dump(res, fh, indent=1, default=str)
    print(f"wrote {args.out}  ({res['seconds']:.1f} s)")


if __name__ == "__main__":
    main()
