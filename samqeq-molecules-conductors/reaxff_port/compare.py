#!/usr/bin/env python3
"""Per-atom charge agreement between a reference dump and one or more samQEq dumps."""
import sys, math

def read(p):
    q, on = {}, False
    for l in open(p):
        if l.startswith('ITEM: ATOMS'):
            on = True
            continue
        if on:
            a = l.split()
            q[int(a[0])] = (int(a[1]), float(a[2]))
    return q

ref = read(sys.argv[1])
for p in sys.argv[2:]:
    c = read(p)
    assert set(c) == set(ref)
    d = [c[k][1] - ref[k][1] for k in ref]
    rel = max(abs(c[k][1] - ref[k][1]) / abs(ref[k][1]) for k in ref)
    print(f"{p}: N={len(d)}  max|dq|={max(map(abs, d)):.3e} e  rms={math.sqrt(sum(x*x for x in d)/len(d)):.3e} e  "
          f"max rel={rel:.3e}  sum q ref={sum(v[1] for v in ref.values()):.2e} sam={sum(v[1] for v in c.values()):.2e}  "
          f"ref q range [{min(v[1] for v in ref.values()):.4f}, {max(v[1] for v in ref.values()):.4f}]")
