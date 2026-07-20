#!/usr/bin/env python3
"""Regenerate Table S1 (complete 71-point LJ grid, 5 functionals) for both the
.tex longtable and the .md table, from lj_fullgrid_cage.csv. Columns:
rho*, T*, region, ref(MBWR/vdH), rig-HS, Lin-2003, Desjarlais, R2PT, 3PT, Delta_3PT.
R2PT is blank for near-solid rho*>1.02 (shown as a dash)."""
import csv
from pathlib import Path

CSV = Path(__file__).resolve().parent.parent / "data" / "lj_fullgrid_cage.csv"
OUT = Path(__file__).resolve().parent

def region(rho):
    if rho <= 0.20: return "gas"
    if rho <= 0.55: return "SC-fluid"
    if rho <= 1.02: return "liquid"
    return "solid"

rows = []
for r in csv.DictReader(open(CSV)):
    rho = float(r["rho_star"]); T = float(r["T_star"])
    rows.append((rho, T, r))
rows.sort(key=lambda x: (x[0], x[1]))

def dash(v): return v if v not in ("", None, "nan", "NaN", "NAN") else None

tex, md = [], []
for rho, T, r in rows:
    reg = region(rho)
    coex = (abs(rho-0.40) < 1e-6 and abs(T-0.90) < 1e-6)   # LV-coexistence outlier
    r2 = dash(r["R2PT"])
    d3 = r["d_3PT"]
    d3_tex = d3.replace("-", "$-$").replace("+", "$+$")
    d3_md  = d3.replace("-", "−")
    Ttex = f"{T:.3f}" + ("$^\\dagger$" if coex else "")
    Tmd  = f"{T:.3f}" + ("$^\\dagger$" if coex else "")
    r2_tex = (r2 if r2 else "---")
    r2_md  = (r2 if r2 else "—")
    tex.append(f"{rho:.3f} & {Ttex} & {reg} & {r['S_ref']} & {r['rigHS']} & "
               f"{r['Lin2003']} & {r['Desj']} & {r2_tex} & {r['3PT']} & {d3_tex} \\\\")
    md.append(f"| {rho:.3f} | {Tmd} | {reg} | {r['S_ref']} | {r['rigHS']} | "
              f"{r['Lin2003']} | {r['Desj']} | {r2_md} | **{r['3PT']}** | {d3_md} |")

n_liq = sum(1 for rho, T, r in rows if 0.70 <= rho <= 1.02)
print(f"# {len(rows)} rows, {n_liq} liquid")
(OUT / "table_s1.tex").write_text("\n".join(tex) + "\n")
(OUT / "table_s1.md").write_text("\n".join(md) + "\n")
print("wrote table_s1.tex and table_s1.md")
