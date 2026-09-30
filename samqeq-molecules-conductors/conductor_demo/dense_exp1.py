#!/usr/bin/env python3
"""Dense spectrum of the experiment-I operators (referee item 2).

opdump.txt (written by a temporary SAMQEQ_OPDUMP diagnostic of the release code, not shipped) holds the columns of
P*H, H being exactly the operator qeq_matvec applies (on-site diagonal, real-space kernel, reciprocal term with the
grid self-term removed). Here: A = Q^T (P H) Q on an orthonormal basis Q of the constrained subspace, dense eigh,
and a check that A reproduces the charges the code solved (b = -chi, no fixed charges, q0 = 0).
Also the analytic Gram bound: lambda_min >= min_i (eta_i + s_i - J_ii(0)), J_ii(0) = K_e (2 alpha_i / pi)^1/2.
"""
import csv, math
from pathlib import Path
import gzip
import numpy as np

HERE = Path(__file__).resolve().parent
KE, LAM = 14.3996454784, 0.462770
RC = {1: 1.618, 3: 0.404355, 4: 0.315328}
CHI = {1: 5.31, 3: 4.50, 4: 7.47}
ETA_W = {3: 15.31, 4: 16.11}
W_SELF = 0.5


def jself(rc):
    a = LAM / (2 * rc * rc)
    return KE * math.sqrt(2 * a / math.pi)


def data_atoms():
    t, m, go = {}, {}, False
    for l in (HERE / "data.exp1").read_text().splitlines():
        s = l.strip()
        if s.startswith("Atoms"): go = True; continue
        if go and s.startswith("Bonds"): break
        if go and s:
            a = s.split(); t[int(a[0])] = int(a[2]); m[int(a[0])] = int(a[1])
    return t, m


def basis(frags):
    """orthonormal basis of {v : sum over each fragment = 0}; frags = list of index lists"""
    n = sum(len(f) for f in frags)
    cols = []
    for f in frags:
        for k in range(1, len(f)):          # Helmert vectors
            v = np.zeros(n); v[f[:k]] = 1.0; v[f[k]] = -k
            cols.append(v / np.linalg.norm(v))
    return np.array(cols).T


def load_charges(d, tags):
    q, go = {}, False
    for l in (d / "charges.dump").read_text().splitlines():
        if l.startswith("ITEM: ATOMS"): go = True; continue
        if go and l.strip():
            a = l.split(); q[int(a[0])] = float(a[3])
    return np.array([q[t] for t in tags])


if __name__ == "__main__":
    typ, mol = data_atoms()
    rows = []
    for arm in ["qeq", "qeq_gself", "mol", "samqeq"]:
        for eta in ["5.172", "2.5", "1.0", "0.3", "0.05"]:
            d = HERE / f"runs/exp1/{arm}_eta{eta}"
            f = d / "opdump.txt"
            lines = (f.read_text() if f.exists() else gzip.open(f.with_suffix(".txt.gz"), "rt").read()).split("\n")
            n = int(lines[0]); tags = [int(x) for x in lines[1].split()]
            PH = np.array([[float(x) for x in lines[2 + j].split()] for j in range(n)]).T   # column j = P H e_j
            if arm.startswith("qeq"):
                frags = [list(range(n))]
            else:
                by = {}
                for k, t in enumerate(tags): by.setdefault(mol[t], []).append(k)
                frags = list(by.values())
            Q = basis(frags)
            A = Q.T @ PH @ Q
            asym = np.abs(A - A.T).max() / np.abs(A).max()
            w = np.linalg.eigvalsh(0.5 * (A + A.T))
            # reproduce the solved charges: minimize over the subspace -> A y = Q^T (-chi)
            chi = np.array([CHI[typ[t]] for t in tags])
            y = np.linalg.solve(0.5 * (A + A.T), Q.T @ (-chi))
            qd = Q @ y
            qc = load_charges(d, tags)
            dq = np.abs(qd - qc).max()
            # Gram bound
            gs = "gself" in arm or arm == "samqeq"
            br = []
            for tt in (1, 3, 4):
                e = float(eta) if tt == 1 else ETA_W[tt]
                s = (KE / (math.sqrt(math.pi) * W_SELF)) if (tt == 1 and gs) else 0.0
                br.append(e + s - jself(RC[tt]))
            rows.append(dict(arm=arm, eta=eta, dim=A.shape[0], lmin=w[0], lmax=w[-1], nneg=int((w < 0).sum()),
                             cond=(w[-1] / w[0] if w[0] > 0 else float("nan")), asym=asym, dq=dq,
                             bound_Au=br[0], bound_H=br[1], bound_M=br[2]))
            r = rows[-1]
            print(f"{arm:<10}{eta:>6}  dim {r['dim']}  lmin {r['lmin']:+8.4f}  lmax {r['lmax']:7.2f}  neg {r['nneg']:3d}"
                  f"  asym {asym:.1e}  |q_dense-q_code| {dq:.1e}   brackets Au/H/M {br[0]:+.3f}/{br[1]:+.3f}/{br[2]:+.3f}")
    with open(HERE / "dense_exp1.csv", "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0])); wr.writeheader(); wr.writerows(rows)
