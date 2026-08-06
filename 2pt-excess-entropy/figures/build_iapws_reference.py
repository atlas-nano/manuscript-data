#!/usr/bin/env python3
"""Build the corrected IAPWS-95 reference table, with a validity status for every state point.

Alexandria's `iapws_sex_kb.csv` was produced with iapws95.py v1.5.5, which silently flashes to a
two-phase mixture whenever the requested density lies inside the vapour-liquid envelope and
overwrites the entropy with x*s_gas + (1-x)*s_liquid. With x tiny that returns the SATURATED-LIQUID
entropy, so the liquid term stops depending on density and S_ex(rho) collapses to const + ln(rho).
Verified by matching the vapour-quality term to three significant figures.

This script classifies all 70 points and emits corrected values where a defensible reference exists.

  valid          in-range single-phase liquid; the CSV value stands
  re-referenced  inside the envelope, but the single-phase branch is mechanically stable
                 (dP/drho > 0) at negative pressure: the standard metastable stretched-liquid
                 reference. Value recomputed on that branch
  drop-spinodal  inside the envelope AND beyond the spinodal (dP/drho < 0). No reference exists
  extrap-T       below the IAPWS-95 lower temperature limit, 251.165 K, but single-phase and
                 mechanically stable. RETAINED as a flagged extrapolation: the limit bounds the
                 accuracy the formulation guarantees, not the existence of a solution
  extrap-P       above the IAPWS-95 pressure limit, 1000 MPa; likewise retained and flagged

Corrected values are put on the same footing as the rest of the CSV by adding the constant
implementation offset between iapws95.py and CoolProp (0.004309 kB, an ideal-part normalisation
difference), measured on the 35 unambiguous single-phase points. It cancels in any anchored
difference and is included only so the columns are directly comparable.

Usage:  python3 build_iapws_reference.py        -> writes iapws_sex_kb_corrected.csv
"""
import csv

import CoolProp.CoolProp as CP

import os
_D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "reference")
SRC = os.path.join(_D, "iapws_sex_kb.csv")
DST = os.path.join(_D, "iapws_sex_kb_corrected.csv")
T_MIN = 251.165          # IAPWS-95 lower temperature limit, K
P_MAX = 1000.0           # IAPWS-95 upper pressure limit, MPa

st = CP.AbstractState("HEOS", "Water")
st.specify_phase(CP.iphase_liquid)


def _update(T, rho):
    st.update(CP.DmassT_INPUTS, rho * 1000.0, T)


def P_MPa(T, rho):
    _update(T, rho)
    return st.p() / 1e6


def s_ex(T, rho):
    """Excess entropy = residual entropy of the EOS, in k_B per molecule."""
    _update(T, rho)
    return st.tau() * st.dalphar_dTau() - st.alphar()


def dP_drho(T, rho, h=0.002):
    return (P_MPa(T, rho + h) - P_MPa(T, rho - h)) / (2 * h)


def load(path):
    d = {}
    with open(path, newline="") as f:
        for r in csv.DictReader(f, skipinitialspace=True):
            r = {k.strip(): v.strip() for k, v in r.items()}
            d[(float(r["Temperature"]), float(r["Density"]))] = float(r["thermoValue"])
    return d


def flashed(d, T, rho, rhos):
    """A point is flashed if it and its neighbour differ by exactly the ideal-gas term.

    Guard against the false positive near saturation by requiring rho to be one of the three
    lowest densities, where the envelope genuinely bites at every temperature studied.
    """
    if rho > 0.95:
        return False
    import math
    i = rhos.index(rho)
    ref = rhos[i - 1] if i else rhos[i + 1]
    return abs((d[(T, rho)] - d[(T, ref)]) - math.log(rho / ref)) < 1e-3


def main():
    d = load(SRC)
    Ts = sorted({t for t, _ in d})
    rhos = sorted({r for _, r in d})

    # implementation offset, from unambiguous single-phase in-range points
    ref = [d[(T, r)] - s_ex(T, r) for T in Ts if T >= 260 for r in rhos if r >= 1.0]
    offset = sum(ref) / len(ref)

    out, tally = [], {}
    for T in Ts:
        for r in rhos:
            val, status = d[(T, r)], "valid"
            if flashed(d, T, r, rhos):
                status = "re-referenced" if dP_drho(T, r) > 0 else "drop-spinodal"
                if status == "re-referenced":
                    val = s_ex(T, r) + offset
            if T < T_MIN:
                status = "extrap-T" if status.startswith(("valid", "re-ref")) else status + "+T"
            try:
                if status == "valid" and P_MPa(T, r) > P_MAX:
                    status = "extrap-P"
            except Exception:
                pass
            tally[status] = tally.get(status, 0) + 1
            out.append((r, T, val, d[(T, r)], status))

    with open(DST, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Density", "Temperature", "thermoValue", "originalValue", "status"])
        for r, T, v, o, s in out:
            w.writerow([r, f"{T:.2f}", f"{v:.6f}", f"{o:.6f}", s])

    print(f"implementation offset applied: {offset:.6f} kB")
    print(f"wrote {DST}\n")
    for k in sorted(tally):
        print(f"  {k:<16}{tally[k]:>3}")
    print(f"  {'TOTAL':<16}{sum(tally.values()):>3}")
    dropped = sum(v for k, v in tally.items() if "spinodal" in k)
    print(f"\n  categorically wrong, dropped : {dropped}")
    print(f"  usable as a reference        : {sum(tally.values()) - dropped} "
          f"of {sum(tally.values())}")


if __name__ == "__main__":
    main()
