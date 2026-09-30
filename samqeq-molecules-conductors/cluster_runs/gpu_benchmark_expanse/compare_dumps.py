#!/usr/bin/env python3
"""Compare two LAMMPS custom dumps (columns: id, then numeric columns) frame by frame.

usage: compare_dumps.py A.dump B.dump [--faces L]
Prints, per shared timestep, max |A-B| for every non-id column. With --faces L (cubic box edge,
dump must carry x y z as its last three columns) also bins the force error (columns fx fy fz) by
distance to the nearest box face, which is the signature of a missing ghost-grid exchange.
"""
import sys


def read(fn):
    frames, cols = {}, None
    lines = open(fn).read().split('\n')
    i = 0
    while i < len(lines):
        if lines[i].startswith('ITEM: TIMESTEP'):
            ts = int(lines[i + 1])
            n = int(lines[i + 3])
            j = i + 1
            while not lines[j].startswith('ITEM: ATOMS'):
                j += 1
            cols = lines[j].split()[2:]
            frames[ts] = {int(l.split()[0]): [float(v) for v in l.split()[1:]]
                          for l in lines[j + 1:j + 1 + n]}
            i = j + 1 + n
        else:
            i += 1
    return frames, cols


def main():
    a, cols = read(sys.argv[1])
    b, _ = read(sys.argv[2])
    names = cols[1:]
    for ts in sorted(set(a) & set(b)):
        A, B = a[ts], b[ts]
        out = []
        for c, name in enumerate(names):
            if name in ('x', 'y', 'z'):
                continue
            out.append(f"{name} {max(abs(A[i][c] - B[i][c]) for i in A):.3e}")
        if 'fx' in names:
            ks = [names.index(k) for k in ('fx', 'fy', 'fz')]
            net = [sum(B[i][k] for i in B) - sum(A[i][k] for i in A) for k in ks]
            out.append("dnetF " + " ".join(f"{v:.3e}" for v in net))
        print(f"step {ts} natoms {len(A)} max|A-B|: " + "  ".join(out))
        if '--faces' in sys.argv and 'fx' in names and 'x' in names:
            L = float(sys.argv[sys.argv.index('--faces') + 1])
            ks = [names.index(k) for k in ('fx', 'fy', 'fz')]
            xs = [names.index(k) for k in ('x', 'y', 'z')]
            near, far = 0.0, 0.0
            for i in A:
                df = max(abs(A[i][k] - B[i][k]) for k in ks)
                d = min(min(A[i][x] % L, L - A[i][x] % L) for x in xs)
                if d < 2.0:
                    near = max(near, df)
                else:
                    far = max(far, df)
            print(f"   force error, atoms < 2 A from a box face: {near:.3e}; farther: {far:.3e}")


if __name__ == '__main__':
    main()
