"""Computes Li-F radial distribution functions and running coordination numbers from data/simulation/conc_ladder."""
import gzip

import numpy as np

from paths import data

RUNGS = [("litpfa_fm_n1", 0.029), ("litpfa_fm_n2", 0.059), ("litpfa_fm_n4", 0.114), ("litpfa_fm_n8", 0.232)]
CN_SIGMA = 0.2
NBINS, RMAX = 400, 12.0
EDGES = np.linspace(0.0, RMAX, NBINS + 1)
CENTERS = 0.5 * (EDGES[:-1] + EDGES[1:])
SHELL_VOL = (4.0 / 3.0) * np.pi * (EDGES[1:] ** 3 - EDGES[:-1] ** 3)
TYPE_LI, TYPE_F_ANION = 11, 9


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def rdf_blocks(path):
    step, rows = None, []
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            t = line.split()
            if len(t) == 2:
                if step is not None and rows:
                    yield step, np.array(rows)
                step, rows = int(t[0]), []
            else:
                rows.append([float(x) for x in t[1:4]])
    if step is not None and rows:
        yield step, np.array(rows)


def fm_cn(tag, skip=1):
    blocks = list(rdf_blocks(data("simulation", "conc_ladder", tag, "rdf.txt")))[skip:]
    a = np.mean([b[1] for b in blocks], axis=0)
    r, g, n = a[:, 0], a[:, 1], a[:, 2]
    w = (r >= 1.5) & (r <= 4.0)
    ipk = np.where(w)[0][np.argmax(g[w])]
    tail = np.where((np.arange(len(r)) > ipk) & (g < g[ipk] / 2.0))[0]
    ihalf = tail[0] if len(tail) else len(r) - 1
    return dict(r=r, g=g, n=n, nblocks=len(blocks), r_half=r[ihalf], cn=n[ihalf])


def frames(path):
    with _open(path) as f:
        while True:
            line = f.readline()
            if not line:
                return
            if not line.startswith("ITEM: TIMESTEP"):
                continue
            ts = int(f.readline())
            f.readline()
            n = int(f.readline())
            f.readline()
            L = []
            for _ in range(3):
                lo, hi = (float(x) for x in f.readline().split()[:2])
                L.append(hi - lo)
            hdr = f.readline().split()[2:]
            col = {k: i for i, k in enumerate(hdr)}
            typ = np.empty(n, int)
            r = np.empty((n, 3))
            for k in range(n):
                p = f.readline().split()
                typ[k] = int(p[col["type"]])
                r[k] = [float(p[col[c]]) for c in ("xu", "yu", "zu")]
            yield ts, (np.array(L), typ, r)


def load_traj(tag):
    keep, order = {}, []
    for ts, fr in frames(data("simulation", "conc_ladder", tag, "ions.lammpstrj")):
        if ts not in keep:
            order.append(ts)
        keep[ts] = fr
    out = []
    for ts in sorted(order):
        L, typ, r = keep[ts]
        out.append((ts, L, r[typ == TYPE_LI], r[typ == TYPE_F_ANION]))
    return out


def min_image_dist(a, b, L):
    d = a[:, None, :] - b[None, :, :]
    d -= L * np.round(d / L)
    return np.linalg.norm(d, axis=-1)


def hist_over_frames(frame_list):
    hist = np.zeros(NBINS)
    n_li_tot, n_an_tot, vols, closest = 0, 0, [], np.inf
    for ts, L, li, an in frame_list:
        d = min_image_dist(li, an, L)
        hist += np.histogram(d.ravel(), bins=EDGES)[0]
        n_li_tot += len(li)
        n_an_tot += len(an)
        vols.append(float(np.prod(L)))
        closest = min(closest, float(d.min()) if d.size else np.inf)
    nframes = len(frame_list)
    return hist, nframes, n_li_tot / nframes, n_an_tot / nframes, float(np.mean(vols)), closest


def gr_nr_from_hist(hist, nframes, n_li, n_an, V):
    per_li = hist / (nframes * n_li)
    g = per_li / ((n_an / V) * SHELL_VOL)
    return g, np.cumsum(per_li)
