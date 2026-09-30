#!/usr/bin/env python3
"""Build the geometries of the conductor demonstration (see DESIGN.md).

exp1: the 240-atom Au(111) slab of the paper's section 2.4 with a TIP4P-FQ water film on its top face.
Types: 1 Au, 2 O, 3 H, 4 M. Water molecules keep mol IDs 1..N; the slab is molecule N+1.
"""
import math, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
AU = Path.home() / "manuscripts_todo/ms_alphasigma_overleaf/conductor_slab/data.au_ortho"
LIQ = Path.home() / "codes/samQEq/examples/tip4pfq_liquid/tip4pfq_liquid.data"

LX, LY, LZ = 17.65065805, 20.38122435, 60.0
Z_TOP = 19.60780130            # top Au plane
GAP = 3.0                      # closest O above the top Au plane
FILM = 7.0                     # thickness of the O window cut from the liquid
OO_MIN = 2.5                   # drop a molecule whose O sits closer than this to a kept O (min image, lateral)


def read_sections(path):
    lines = Path(path).read_text().splitlines()
    sec, out = None, {}
    for l in lines:
        s = l.strip()
        w = s.split()[0] if s else ""
        if w in ("Masses", "Atoms", "Bonds", "Angles", "Velocities"):
            sec = w; out[sec] = []; continue
        if sec and s and not s.startswith("#"):
            out[sec].append(s.split())
    return out


def box_len(path):
    for l in Path(path).read_text().splitlines():
        if "xlo xhi" in l:
            a, b = map(float, l.split()[:2]); return b - a
    raise SystemExit("no box")


def water_film():
    L = box_len(LIQ)
    d = read_sections(LIQ)
    mols = {}
    for a in d["Atoms"]:
        i, m, t = int(a[0]), int(a[1]), int(a[2])
        x, y, z = map(float, a[4:7])
        mols.setdefault(m, []).append((i, t, x, y, z))
    # whole molecules, positions relative to their O (unwrap H/M by minimum image)
    rigid = []
    for m, at in mols.items():
        o = [a for a in at if a[1] == 1][0]
        rel = []
        for (i, t, x, y, z) in sorted(at, key=lambda a: a[0]):
            dx = [(p - q) - L * round((p - q) / L) for p, q in ((x, o[2]), (y, o[3]), (z, o[4]))]
            rel.append((t, dx))
        rigid.append(((o[2] % L, o[3] % L, o[4] % L), rel))
    # tile laterally 2x2, take O inside the Au cell and a FILM-thick window in z
    z0 = 4.0
    cand = []
    for (ox, oy, oz), rel in rigid:
        for ix in range(2):
            for iy in range(2):
                X, Y = ox + ix * L, oy + iy * L
                if X < LX and Y < LY and z0 <= oz < z0 + FILM:
                    cand.append(((X, Y, oz), rel))
    kept = []
    for (p, rel) in sorted(cand, key=lambda c: c[0][2]):
        ok = True
        for (q, _) in kept:
            dx = p[0] - q[0]; dx -= LX * round(dx / LX)
            dy = p[1] - q[1]; dy -= LY * round(dy / LY)
            if dx * dx + dy * dy + (p[2] - q[2]) ** 2 < OO_MIN ** 2:
                ok = False; break
        if ok:
            kept.append((p, rel))
    # shift so the lowest O is GAP above the top Au plane
    zmin = min(p[2] for p, _ in kept)
    shift = Z_TOP + GAP - zmin
    return [((p[0], p[1], p[2] + shift), rel) for p, rel in kept]


def write_exp1(out):
    au = read_sections(AU)["Atoms"]
    film = water_film()
    atoms, bonds, angles = [], [], []
    aid = 0
    for n, (p, rel) in enumerate(film, start=1):
        ids = {}
        for (t, d) in rel:
            aid += 1
            typ = {1: 2, 2: 3, 3: 4}[t]
            atoms.append((aid, n, typ, 0.0, p[0] + d[0], p[1] + d[1], p[2] + d[2]))
            ids.setdefault(t, []).append(aid)
        o, (h1, h2) = ids[1][0], ids[2]
        bonds += [(1, o, h1), (1, o, h2)]
        angles += [(1, h1, o, h2)]
    nw = len(film)
    for a in au:
        aid += 1
        atoms.append((aid, nw + 1, 1, 0.0, float(a[4]), float(a[5]), float(a[6])))
    # sanity: closest water atom to any Au atom
    wat = [a for a in atoms if a[2] != 1]
    dmin = min(math.dist(w[4:7], (float(a[4]), float(a[5]), float(a[6]))) for w in wat for a in au
               if abs(w[6] - float(a[6])) < 4.0)
    with open(out, "w") as f:
        f.write(f"# Au(111) 5L slab + TIP4P-FQ water film ({nw} molecules); built by build.py\n\n")
        f.write(f"{len(atoms)} atoms\n{len(bonds)} bonds\n{len(angles)} angles\n\n")
        f.write("4 atom types\n1 bond types\n1 angle types\n\n")
        f.write(f"0.0 {LX:.8f} xlo xhi\n0.0 {LY:.8f} ylo yhi\n0.0 {LZ:.8f} zlo zhi\n\n")
        f.write("Masses\n\n1 196.966569\n2 15.9994\n3 1.008\n4 1.0e-6\n\nAtoms # full\n\n")
        for a in atoms:
            f.write("%d %d %d %.6f %.8f %.8f %.8f\n" % a)
        f.write("\nBonds\n\n")
        for k, b in enumerate(bonds, 1):
            f.write("%d %d %d %d\n" % (k, *b))
        f.write("\nAngles\n\n")
        for k, b in enumerate(angles, 1):
            f.write("%d %d %d %d %d\n" % (k, *b))
    zw = [a[6] for a in wat]
    print(f"exp1: {nw} waters, {len(atoms)} atoms, water z {min(zw):.2f}-{max(zw):.2f}, "
          f"closest water-Au {dmin:.2f} A")
    return nw


if __name__ == "__main__":
    write_exp1(HERE / "data.exp1")
