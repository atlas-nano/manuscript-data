"""Computes the van der Waals volumes, equal-volume radii and size ratios k of tables S6 and S7 from data/simulation/vdw_structures/*.xyz."""
import argparse, collections, csv, math, pathlib, subprocess, tempfile
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
GEOM = ROOT / "data" / "simulation" / "vdw_structures"
RADII = {"H": 1.20, "Li": 1.82, "B": 1.92, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47, "Al": 1.84, "P": 1.80, "S": 1.80}
ANIONS = ["BF4", "PF6", "FSI", "TFSI", "TPFA"]
SHELLS = ["Li_ACN", "Li_DME", "Li_Me2O", "Li_DEE", "Li_FM"]
SPACINGS = [0.12, 0.08, 0.05, 0.035]


def read_xyz(path):
    lines = path.read_text().splitlines()
    n = int(lines[0])
    rows = [x.split() for x in lines[2:] if x.strip()]
    if len(rows) != n:
        raise ValueError(f"{path.name}: {len(rows)} atom records, header says {n}")
    elems = [r[0] for r in rows]
    unknown = set(elems) - RADII.keys()
    if unknown:
        raise ValueError(f"{path.name}: no radius for {sorted(unknown)}")
    coords = np.array([[float(v) for v in r[1:4]] for r in rows])
    if not np.isfinite(coords).all() or len({(e, *p) for e, p in zip(elems, map(tuple, coords))}) != n:
        raise ValueError(f"{path.name}: non-finite or duplicate atom records")
    return elems, coords


def voxel_volume(coords, radii, h):
    low = np.min(coords - radii[:, None], axis=0) - h * 2
    high = np.max(coords + radii[:, None], axis=0) + h * 2
    dims = np.ceil((high - low) / h).astype(int)
    grid = np.zeros(tuple(dims), dtype=bool)
    for pos, r in zip(coords, radii):
        a = np.maximum(0, np.floor((pos - r - low) / h).astype(int))
        b = np.minimum(dims, np.ceil((pos + r - low) / h).astype(int) + 1)
        axes = [low[j] + (np.arange(a[j], b[j]) + .5) * h - pos[j] for j in range(3)]
        mask = axes[0][:, None, None] ** 2 + axes[1][None, :, None] ** 2 + axes[2][None, None, :] ** 2 < r * r
        grid[a[0]:b[0], a[1]:b[1], a[2]:b[2]] |= mask
    return float(grid.sum() * h ** 3)


def grid100_volume(coords, radii):
    low = coords.min(0) - 5
    high = coords.max(0) + 5
    axes = [np.linspace(low[j], high[j], 100) for j in range(3)]
    occ = np.zeros((100, 100, 100), bool)
    for p, r in zip(coords, radii):
        occ |= ((axes[0][:, None, None] - p[0]) ** 2 + (axes[1][None, :, None] - p[1]) ** 2
                + (axes[2][None, None, :] - p[2]) ** 2) < r * r
    return float(occ.sum() * np.prod((high - low) / 99))


def fragments(obabel, elems, coords):
    with tempfile.TemporaryDirectory() as d:
        xyz, sdf = pathlib.Path(d) / "s.xyz", pathlib.Path(d) / "s.sdf"
        xyz.write_text(f"{len(elems)}\n\n" + "".join(f"{e} {p[0]:.9f} {p[1]:.9f} {p[2]:.9f}\n" for e, p in zip(elems, coords)))
        subprocess.run([obabel, str(xyz), "-O", str(sdf)], check=True, capture_output=True)
        lines = sdf.read_text().splitlines()
    na, nb = int(lines[3][:3]), int(lines[3][3:6])
    adj = [set() for _ in elems]
    for b in lines[4 + na:4 + na + nb]:
        i, j = int(b[:3]) - 1, int(b[3:6]) - 1
        if elems[i] != "Li" and elems[j] != "Li":
            adj[i].add(j)
            adj[j].add(i)
    seen, out = set(), []
    for i in range(na):
        if i in seen:
            continue
        todo, group = [i], []
        seen.add(i)
        while todo:
            j = todo.pop()
            group.append(j)
            for k in adj[j] - seen:
                seen.add(k)
                todo.append(k)
        out.append(dict(sorted(collections.Counter(elems[j] for j in group).items())))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=pathlib.Path, default=ROOT / "data" / "tables")
    ap.add_argument("--openbabel", help="Open Babel executable; when given, the bonded fragments of each structure are printed")
    a = ap.parse_args()
    a.out_dir.mkdir(parents=True, exist_ok=True)
    results, conv = {}, []
    for sid in ANIONS + SHELLS:
        elems, coords = read_xyz(GEOM / f"{sid}.xyz")
        radii = np.array([RADII[e] for e in elems])
        vv = []
        for h in SPACINGS:
            v = voxel_volume(coords, radii, h)
            vv.append(v)
            conv.append([sid, h, v, (3 * v / (4 * math.pi)) ** (1 / 3)])
        results[sid] = dict(species=sid, volume_A3=vv[-1], effective_radius_A=(3 * vv[-1] / (4 * math.pi)) ** (1 / 3),
                            volume_grid100_A3=grid100_volume(coords, radii),
                            change_from_h005_percent=abs(vv[-1] / vv[-2] - 1) * 100)
        line = f"{sid:8s} V = {vv[-1]:9.3f} A^3  r = {results[sid]['effective_radius_A']:.4f} A"
        if a.openbabel:
            line += f"  fragments {fragments(a.openbabel, elems, coords)}"
        print(line, flush=True)
    with open(a.out_dir / "vdw_volumes.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(results["BF4"]), lineterminator="\n")
        w.writeheader()
        w.writerows(results.values())
    with open(a.out_dir / "vdw_grid_convergence.csv", "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["species", "spacing_A", "volume_A3", "effective_radius_A"])
        w.writerows(conv)
    with open(a.out_dir / "vdw_k.csv", "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["anion", "solvent", "r_minus_A", "r_plus_A", "k", "ideal_k_over_k_plus_1"])
        for an in ANIONS:
            for sh in SHELLS:
                ra, rs = results[an]["effective_radius_A"], results[sh]["effective_radius_A"]
                k = ra / rs
                w.writerow([an, sh[3:], ra, rs, k, k / (1 + k)])


if __name__ == "__main__":
    main()
