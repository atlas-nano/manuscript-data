"""Standalone kernel + cage-entropy helpers for the 3PT manuscript figures.

Vendored from the py-xPT thermodynamics module so the figure scripts in this
repository run without installing py-xPT (dependencies: numpy, scipy only).

Provides:
  invert_kernel_matrix  -- Form-B discrete Volterra inversion of the VACF -> K(t)
  cage_memory_entropy   -- the parameter-free 3PT cage entropy correction
  packing_from_f_dgen, hs_excess_entropy_dgen, hs_contact_dgen
                        -- d-dimensional hard-sphere references
  safe_matrix_inverse   -- condition-number-guarded matrix inverse

The d=4 hard-sphere branches are omitted; the manuscript figures use d=1,2,3.
"""
from __future__ import annotations
import logging
import math
import numpy as np

log = logging.getLogger(__name__)

# physical constants (SI)
NA     = 6.0221367e23       # Avogadro (mol^-1)
KB     = 1.380658e-23       # Boltzmann (J/K)
H      = 6.62606896e-34     # Planck (J*s)
R      = KB * NA            # gas constant (J/mol/K)
VLIGHT = 2.99792458e8       # speed of light (m/s)
PI     = math.pi


def safe_matrix_inverse(A: np.ndarray, *,
                        cond_warn: float = 1e10,
                        rcond_pinv: float = 1e-10,
                        label: str | None = None,
                        ) -> tuple[np.ndarray, float, bool]:
    """Condition-number-guarded matrix inverse with pseudo-inverse fallback.

    Returns ``np.linalg.pinv(A, rcond=rcond_pinv)`` (with a WARNING) when the
    condition number exceeds ``cond_warn`` or the direct inverse raises;
    otherwise ``np.linalg.inv(A)``.

    Returns
    -------
    A_inv : (m, m) np.ndarray -- the inverse or pseudo-inverse.
    cond  : float             -- the condition number (inf if exactly singular).
    fallback : bool           -- True if the pinv path was taken.
    """
    try:
        cond = float(np.linalg.cond(A))
    except (np.linalg.LinAlgError, ValueError):
        cond = float('inf')
    if not np.isfinite(cond) or cond > cond_warn:
        if label is not None:
            log.warning("safe_matrix_inverse(%s): condition number %.3g exceeds "
                        "%.0e -- using pseudo-inverse with rcond=%.0e.",
                        label, cond, cond_warn, rcond_pinv)
        return np.linalg.pinv(A, rcond=rcond_pinv), cond, True
    try:
        return np.linalg.inv(A), cond, False
    except np.linalg.LinAlgError:
        if label is not None:
            log.warning("safe_matrix_inverse(%s): direct inverse failed; using "
                        "pseudo-inverse with rcond=%.0e.", label, rcond_pinv)
        return np.linalg.pinv(A, rcond=rcond_pinv), cond, True


def hs_excess_entropy_dgen(eta: float, d: int) -> float:
    """Rigorous hard-sphere excess entropy S^ex/Nk_B in d dimensions
    (closed forms d=1,2,3)."""
    if eta <= 0.0:
        return 0.0
    if d == 3:
        return eta * (3.0 * eta - 4.0) / (1.0 - eta)**2                          # Carnahan-Starling
    if d == 2:
        return -(9.0 / 8.0) * eta / (1.0 - eta) + (7.0 / 8.0) * math.log(1.0 - eta)  # Henderson disk
    if d == 1:
        return math.log(1.0 - eta)                                               # Tonks rod
    raise NotImplementedError(f"dimension {d} not bundled (d=1,2,3 only)")


def hs_contact_dgen(gamma: float, d: int) -> float:
    """Hard-sphere contact value g_d(eta+) in d dimensions (d=1,2,3)."""
    if d == 3:
        return (1.0 - gamma / 2.0) / (1.0 - gamma)**3        # Carnahan-Starling
    if d == 2:
        return (1.0 - 7.0 * gamma / 16.0) / (1.0 - gamma)**2  # Henderson disk
    if d == 1:
        return 1.0 / (1.0 - gamma)                            # Tonks rod
    raise NotImplementedError(f"dimension {d} not bundled (d=1,2,3 only)")


def packing_from_f_dgen(f: float, d: int) -> float:
    """Packing fraction gamma from fluidicity f by inverting the Enskog contact
    relation 1/f = g_d(gamma).  At d=3 reproduces the standard 2PT packing y."""
    from scipy.optimize import brentq
    f = float(max(1e-9, min(f, 1.0 - 1e-9)))
    gmax = {1: 0.999, 2: 0.9, 3: 0.74}[d]
    try:
        return float(brentq(lambda g: hs_contact_dgen(g, d) - 1.0 / f,
                            1e-12, gmax, xtol=1e-13, rtol=1e-12, maxiter=200))
    except ValueError:
        return gmax


def invert_kernel_matrix(dt: float, C_matrix: np.ndarray,
                          *,
                          noise_floor: float = 1e-3,
                          divergence_factor: float = 1e6,
                          swing_threshold: float = 10.0,
                          swing_window: int = 5,
                          nf_run: int = 1,
                          label: str | None = None,
                          ) -> np.ndarray:
    """Discrete Volterra inversion of the matrix memory kernel via the
    second-derivative (Form B) GLE.

    Normalising the VACF by C(0) so Chat(0)=I, the kernel obeys
        Chat''(t) = -K(t) - int_0^t K(t-s) Chat'(s) ds,
    with the even-extension boundary value K(0) = 2(Chat(0)-Chat(dt))/dt^2 and
    the central-difference recursion
        K[n] = -Chat''[n] - sum_{j=1..n} K[n-j] (Chat[j+1]-Chat[j-1])/2.

    Parameters
    ----------
    C_matrix : (N, m, m) ndarray -- velocity autocorrelation matrix.
    dt : float -- time step [ps].
    noise_floor, divergence_factor, swing_threshold, swing_window :
        per-lag truncation guards (stop at the noise floor / divergence /
        instability swing-streak).
    nf_run : int -- consecutive sub-noise-floor lags required before truncating
        (>=3 makes the guard envelope-aware for oscillatory/librational VACFs;
        1 = first-lag behaviour).

    Returns
    -------
    M : (N, m, m) ndarray -- memory kernel [1/ps^2]; zero past truncation.
    """
    n_lags = C_matrix.shape[0]
    m      = C_matrix.shape[1]
    M = np.zeros((n_lags, m, m))

    C0_inv, _C0_cond, _C0_pinv = safe_matrix_inverse(
        C_matrix[0],
        label=(f"C(0) for invert_kernel_matrix({label})" if label else None))

    C_hat = np.einsum('nij,jk->nik', C_matrix, C0_inv)

    # boundary K(0) from even-extension second difference (Form B)
    M[0] = 2.0 * (C_hat[0] - C_hat[1]) / (dt**2)

    M0_max = float(np.max(np.abs(M[0])))
    if M0_max <= 0.0:
        M0_max = 1.0

    _swing_streak = 0
    _nf_streak    = 0
    _M_prev_norm  = float(np.linalg.norm(M[0]))
    _truncation = None

    for n in range(1, n_lags - 1):
        # (1) Frobenius noise floor on Chat[n]; with nf_run>1 require a streak.
        c_norm = float(np.linalg.norm(C_hat[n]))
        if c_norm / m < noise_floor:
            _nf_streak += 1
            if _nf_streak >= nf_run:
                n0 = n - _nf_streak + 1
                M[n0:n] = 0.0
                _truncation = (n0, "noise_floor", c_norm, 0.0)
                break
        else:
            _nf_streak = 0

        Cdd_n = (C_hat[n-1] - 2.0 * C_hat[n] + C_hat[n+1]) / (dt ** 2)

        conv = np.zeros((m, m))
        for j in range(1, n + 1):
            Cdot_j = (C_hat[j+1] - C_hat[j-1]) / 2.0
            conv += M[n-j] @ Cdot_j

        # (2) conv-finite check
        if not np.isfinite(conv).all():
            _truncation = (n, "conv_nonfinite", c_norm, 0.0)
            break

        with np.errstate(over='ignore', invalid='ignore'):
            M_n = -Cdd_n - conv

        if not np.isfinite(M_n).all():
            _truncation = (n, "k_nonfinite", c_norm, 0.0)
            break

        # (3) per-element divergence clamp
        k_max = float(np.max(np.abs(M_n)))
        if k_max > divergence_factor * M0_max:
            _truncation = (n, "divergence_clamp", c_norm, k_max)
            break

        # (4) adaptive swing-streak on Frobenius norm
        M_n_norm = float(np.linalg.norm(M_n))
        if n > 1 and _M_prev_norm > 1e-30:
            ratio = M_n_norm / _M_prev_norm
            if ratio > swing_threshold or ratio < 1.0 / swing_threshold:
                _swing_streak += 1
                if _swing_streak >= swing_window:
                    _truncation = (n, "swing_instability", c_norm, k_max)
                    break
            else:
                _swing_streak = 0

        M[n] = M_n
        _M_prev_norm = M_n_norm

    if _truncation is not None and label is not None:
        idx, reason, c_val, k_val = _truncation
        log.info("invert_kernel_matrix(%s) truncated at lag %d (t=%.4g ps): %s",
                 label, idx, idx * dt, reason)
    return M


def cage_memory_entropy(dt: float, C_scalar: np.ndarray,
                        nu_cm: np.ndarray, dos_total: np.ndarray,
                        dos_gas: np.ndarray, T_K: float,
                        mass_amu: float, vol_per_atom_A3: float,
                        *,
                        prefactor: float | None = 1.0 / 3.0,
                        dimension: int = 3,
                        gate_f0: float = 0.01,
                        clip_eps: float = 1e-3,
                        nuc_scale: float = 1.0,
                        mainlobe_alpha: float = 0.02,
                        tail_tol: float = 1.0,
                        Wg_override: float | None = None,
                        label: str | None = None,
                        cage_out: list | None = None) -> float | None:
    """Parameter-free cage-memory entropy correction dS [k_B per atom] to the
    rigorous-HS 2PT entropy:

        dS = p * g(f) * int cage(nu) * (1 - w(nu)) * (W_g - W_s(nu)) d_nu ,

    where cage(nu) = const*(F_K - F_M) is the memory excess of the friction
    kernel (F_K = Re[1/(i*omega + Ktilde)] minus its Markovian equivalent
    F_M = Re[1/(i*omega + gamma)], gamma = Ktilde(0)), clipped to [0, solid2]
    with solid2 = max(dos_total - dos_gas, 0); w(nu) = nu^2/(nu^2 + nu_c^2)
    high-passes the harmonic tail (nu_c = Omega_0/2*pi*c, the Einstein
    frequency); W_g is the rigorous-HS gas weight per DoF (Sackur-Tetrode +
    Carnahan-Starling excess, no ln z) and W_s = 1 - ln(hc*nu/kT) is harmonic.
    The prefactor defaults to p = 1/d; g(f) = f^2/(f^2 + f0^2) gates the cage off
    in the crystalline (f->0) limit.

    For a non-translational (e.g. rotational) channel, pass ``Wg_override`` to
    supply the per-DoF gas weight directly (the free-rotor weight) instead of
    the hard-sphere Sackur-Tetrode form.

    Returns dS, or None if the kernel inversion is degenerate or f is
    unphysical.  If ``cage_out`` is given, the clipped cage DoS array is
    appended to it on success.
    """
    from scipy.signal import czt

    nu = np.asarray(nu_cm, dtype=float)
    tot = np.asarray(dos_total, dtype=float)
    gas = np.asarray(dos_gas, dtype=float)
    C = np.asarray(C_scalar, dtype=float)
    if nu.size < 2 or C.size < 2 or C[0] == 0.0:
        return None
    dnu = nu[1] - nu[0]

    # smooth (C-infinity) clip via softplus/LogSumExp; width = clip_eps*S(0)
    # (clip_eps -> 0 recovers the hard clip max[min(.,solid2),0]).
    _beta = (1.0 / (clip_eps * tot[0])) if (clip_eps > 0.0 and tot[0] > 0.0) else None
    def _smax0(x):            # smooth max(x, 0)
        return x if _beta is None else np.logaddexp(0.0, _beta * x) / _beta
    def _smin(a, b):          # smooth min(a, b)
        return np.minimum(a, b) if _beta is None else -np.logaddexp(-_beta * a, -_beta * b) / _beta

    solid2 = _smax0(tot - gas)

    hc_k = 100.0 * H * VLIGHT / KB                # hc/k_B [cm*K] ~ 1.43877
    u = np.where(nu > 0, hc_k * nu / T_K, 1e-9)
    Ws = np.where(nu > 0, 1.0 - np.log(u), 0.0)

    c_cm_ps = VLIGHT * 1e-10                       # speed of light [cm/ps]
    wa = 2.0 * PI * nu * c_cm_ps                   # angular frequency [1/ps]

    # scalar (trace/3) Volterra kernel Ktilde(omega) via the time-domain round-trip
    cn = C / C[0]
    Cm = np.zeros((cn.size, 3, 3))
    for i in range(3):
        Cm[:, i, i] = cn
    K = np.einsum('tii->t', invert_kernel_matrix(dt, Cm)) / 3.0
    nz = np.nonzero(K)[0]
    if not nz.size:
        if label:
            log.warning("cage_memory_entropy(%s): degenerate kernel", label)
        return None
    nK_auto = int(nz[-1]) + 1

    # Auto support is the right limit for a cleanly-decaying kernel, but a
    # smooth spurious tail beyond the main lobe inflates gamma=Ktilde(0) and
    # hence F_M, over-counting the cage.  Detect it by comparing the friction at
    # the auto cutoff with that at the main-lobe cutoff (first |K|<alpha*|K(0)|)
    # and fall back to the main-lobe cutoff when the tail shifts gamma by more
    # than tail_tol.
    def _gamma_dc(nk):                       # Ktilde(0) = int K dt
        if nk < 2:
            return float(K[0]) * dt
        kv = K[:nk].copy(); kv[0] *= 0.5; kv[-1] *= 0.5
        return float(kv.sum()) * dt
    K0 = abs(K[0])
    _below = np.where(np.abs(K[3:]) < mainlobe_alpha * K0)[0] if K0 > 0.0 else np.array([], dtype=int)
    nK_main = (int(_below[0]) + 3) if _below.size else nK_auto
    nK = nK_auto
    if nK_main < nK_auto and tail_tol is not None:
        g_auto, g_main = _gamma_dc(nK_auto), _gamma_dc(nK_main)
        if abs(g_main) > 0.0 and abs(g_auto - g_main) > tail_tol * abs(g_main):
            log.warning("cage_memory_entropy%s: auto-cutoff friction inflated "
                        "%.2gx by a smooth post-main-lobe kernel tail; falling "
                        "back to the main-lobe cutoff.",
                        f"({label})" if label else "", g_auto / g_main)
            nK = nK_main

    Kv = K[:nK].copy()
    Kv[0] *= 0.5
    Kv[-1] *= 0.5
    Ktil = czt(Kv, m=wa.size,
               w=np.exp(-1j * 2.0 * PI * dnu * c_cm_ps * dt), a=1.0) * dt

    gamma = float(Ktil[0].real)
    const = tot[0] / np.trapezoid(cn, dx=dt)
    F_K = np.real(1.0 / (1j * wa + Ktil))
    Om0sq = float(np.trapezoid(wa**2 * tot, dx=1.0) / np.trapezoid(tot, dx=1.0))
    F_ref = gamma / (gamma**2 + wa**2)
    cage = _smax0(_smin(const * (F_K - F_ref), solid2))

    Om0 = np.sqrt(Om0sq)
    nuc = nuc_scale * Om0 / (2.0 * PI * c_cm_ps)   # nuc_scale=1: Einstein cutoff
    w = nu**2 / (nu**2 + nuc**2)

    f = float(np.trapezoid(gas, dx=dnu) / dimension)   # gas DoF fraction (int gas = d*f)
    if not (0.0 < f < 1.0):
        if label:
            log.warning("cage_memory_entropy(%s): unphysical fluidicity f=%.3f",
                        label, f)
        return None
    if Wg_override is not None:
        # non-translational channel: per-DoF gas weight supplied directly
        Wg = float(Wg_override)
    else:
        y = packing_from_f_dgen(f, dimension)      # HS packing (d-dim contact relation)
        m_kg = mass_amu * 1e-3 / NA
        lam_invd = (2.0 * PI * m_kg * KB * T_K / H**2)**(dimension / 2.0)   # lambda^-d [1/m^d]
        Vpa = vol_per_atom_A3 * (1e-10)**dimension     # [m^d]
        Wg = ((dimension / 2.0 + 1.0 + np.log(lam_invd * Vpa / f)) / dimension
              + hs_excess_entropy_dgen(y, dimension) / dimension)

    p = f if prefactor is None else prefactor
    if cage_out is not None:
        cage_out.append(cage)          # per-atom cage DoS for the .pwr writer

    # fluidicity gate (enforces the f->0 Debye-Einstein crystalline limit)
    gate = (f * f) / (f * f + gate_f0 * gate_f0)
    dS = float(p * gate * np.trapezoid(cage * (1.0 - w) * (Wg - Ws), dx=dnu))
    return dS
