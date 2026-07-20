#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Na+-only solvation analysis using your `solvation_analysis` package.

Outputs (exactly 3 Na+-wide CSVs):
  1) Na+_rdf.csv       -> r_A plus g(r) columns for each solvent site
  2) Na+_density.csv   -> r_A plus density ρ(r) columns for each solvent site
  3) Na+_residence.csv -> per-solvent summary: solvent, CN, solvation_radius_A, residence_time_cutoff, residence_time_fit

Additional per-site artifacts:
  - Na+_rdf_<site>.png
  - Na+_density_<site>.png
  - Na+_residence_<site>.png        (C(t) line plot)
  - Na+_residenceCt_<site>.csv      (frame, C_t)

Usage
-----
python solvation_cli_na_only.py \
  --data-file /global/cfs/cdirs/m4248/xiaoxusr/2025/na_conductor/tetramer_glymes/classical_fix_momentum/glyme_tetramer_8_1/results_restart/lammps.data \
  --traj /global/cfs/cdirs/m4248/xiaoxusr/2025/na_conductor/tetramer_glymes/classical_fix_momentum/glyme_tetramer_8_1/results_restart/lammps.298K.prod.lammpsdump \
  --step 100 \
  --outdir ./solvation_out_2
"""

from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import MDAnalysis as mda

# Your package
import solvation_analysis
from solvation_analysis.solute import Solute
from solvation_analysis.residence import Residence

# ---------------- CLI ----------------
def parse_args():
    p = argparse.ArgumentParser(description="Na+-only solvation analysis (RDF, density, residence) with `solvation_analysis`.")
    p.add_argument("--data-file", default="./lammps.data", help="Topology file (LAMMPS data)")
    p.add_argument("--traj", nargs="+",
                   default=["/global/cfs/cdirs/m4248/xiaoxusr/na_conductor/tetramer_glymes/classical_fix_momentum/glyme_tetramer_8_1/results_restart/lammps.298K.prod.lammpsdump",
                            "./lammps.298K.prod.lammpsdump"],
                   help="Trajectory file(s)")
    p.add_argument("--step", type=int, default=100, help="Stride for RDF/residence analysis in frames")
    p.add_argument("--outdir", default="./solvation_sa_out", help="Directory to save outputs")
    return p.parse_args()

# ------------- Helpers --------------
def ensure_outdir(path: str|Path) -> Path:
    p = Path(path); p.mkdir(parents=True, exist_ok=True); return p

def detect_resid_groups(u: mda.Universe):
    """
    Notebook heuristic:
      glyme => residue with 16 atoms
      Na+   => residue with 1 atom
      else  => tetramer
    Returns dict of lists of RESID (1-based) for 'glyme', 'Na+', 'tetramer'.
    """
    groups = {"glyme": [], "Na+": [], "tetramer": []}
    for ix in u.residues.ix:  # 0-based indices
        resid = int(ix) + 1
        n_atoms = u.select_atoms(f"resid {resid}").n_atoms
        if n_atoms == 16:
            groups["glyme"].append(resid)
        elif n_atoms == 1:
            groups["Na+"].append(resid)
        else:
            groups["tetramer"].append(resid)
    return groups

def build_atomgroups(u: mda.Universe, resid_groups: dict):
    glyme_resid = " ".join(map(str, resid_groups["glyme"])) if resid_groups["glyme"] else ""
    tetramer_resid = " ".join(map(str, resid_groups["tetramer"])) if resid_groups["tetramer"] else ""

    selections = {
        # tetramer sites
        "tetramer_N":    f"(type 8) and (resid {tetramer_resid})" if tetramer_resid else "none",
        "tetramer_OS1":  f"(type 10) and (charge -0.529 to -0.528) and (resid {tetramer_resid})" if tetramer_resid else "none",
        "tetramer_OS2":  f"(type 10) and (charge -0.520 to -0.519) and (resid {tetramer_resid})" if tetramer_resid else "none",
        "tetramer_O":    f"(type 2) and (charge -0.05 to -0.04) and (resid {tetramer_resid})" if tetramer_resid else "none",
        # "tetramer_S6":   f"type 9 and (resid {tetramer_resid})" if tetramer_resid else "none",
        # "tetramer_SY":   f"type 7 and (resid {tetramer_resid})" if tetramer_resid else "none",
        # "tetramer_SS":   f"type 4 and (resid {tetramer_resid})" if tetramer_resid else "none",
        "tetramer_F":    "type 5",
        # glyme O
        "glyme_O":       f"type 2 and (resid {glyme_resid})" if glyme_resid else "none",
        # solute
        "Na+":           "type 11",
    }

    ag = {}
    for key, sel in selections.items():
        if sel == "none":
            continue
        grp = u.select_atoms(sel)
        if grp.n_atoms > 0:
            ag[key] = grp
    return ag

def plot_multi(x, curves_dict, xlabel, ylabel, title, outpath):
    plt.figure()
    for label, y in curves_dict.items():
        plt.plot(x, y, lw=1.8, label=label)
    plt.xlabel(xlabel); plt.ylabel(ylabel); plt.title(title)
    if len(curves_dict) > 1:
        plt.legend()
    plt.tight_layout()
    plt.savefig(outpath, dpi=200); plt.close()

# ------------- Main --------------
def main():
    args = parse_args()
    outdir = ensure_outdir(args.outdir)

    # Load system
    u = mda.Universe(args.data_file, *args.traj)

    # Build residue groups & selections
    grp = detect_resid_groups(u)
    ag = build_atomgroups(u, grp)

    if "Na+" not in ag:
        raise ValueError("No Na+ atoms found. Check topology/trajectory and type mapping.")

    # --------- Prepare Na+ analysis ---------
    # Solvents to include (present in ag)
    na_solvent_keys = ["tetramer_N","tetramer_OS1","tetramer_OS2","tetramer_O","tetramer_F","glyme_O"]
    na_solvents = {k: ag[k] for k in na_solvent_keys if k in ag}

    # Radii dicts (only keys present are used by solvation_analysis)
    na_radii_density = {"tetramer_N":5.5,"tetramer_O":7.5,"tetramer_F":5.5}
    na_radii_rdf     = {"tetramer_N":5.5,"tetramer_O":7.5,"tetramer_F":5.5}

    # --------- Density mode ---------
    sol_den = Solute.from_atoms(ag["Na+"], na_solvents, radii=na_radii_density,
                                solute_name="Na+", rdf_init_kwargs={"norm":"density"})
    sol_den.run(step=args.step)
    data_density = sol_den.rdf_data["Na+"] if "Na+" in sol_den.rdf_data else sol_den.rdf_data["solute_0"]

    r_d = None
    density_cols = {}
    for key, (r_arr, y_arr) in data_density.items():
        if r_d is None:
            r_d = np.asarray(r_arr)
        density_cols[key] = np.asarray(y_arr)

    df_density = pd.DataFrame({"r_A": r_d, **density_cols})
    df_density.to_csv(Path(outdir) / "Na+_density.csv", index=False)
    # Per-site density plots
    for _site, _y in density_cols.items():
        plt.figure()
        plt.plot(r_d, _y, lw=1.8)
        plt.xlabel("r (Å)"); plt.ylabel(r"ρ(r) (Å$^{-3}$)")
        plt.title(f"Na+ density: {_site}")
        plt.tight_layout()
        plt.savefig(Path(outdir) / f"Na+_density_{_site}.png", dpi=200)
        plt.close()

    # --------- RDF mode ---------
    sol_rdf = Solute.from_atoms(ag["Na+"], na_solvents, radii=na_radii_rdf,
                                solute_name="Na+", rdf_init_kwargs={"norm":"rdf"})
    sol_rdf.run(step=args.step)
    data_rdf = sol_rdf.rdf_data["Na+"] if "Na+" in sol_rdf.rdf_data else sol_rdf.rdf_data["solute_0"]

    r_g = None
    rdf_cols = {}
    for key, (r_arr, y_arr) in data_rdf.items():
        if r_g is None:
            r_g = np.asarray(r_arr)
        rdf_cols[key] = np.asarray(y_arr)

    df_rdf = pd.DataFrame({"r_A": r_g, **rdf_cols})
    df_rdf.to_csv(Path(outdir) / "Na+_rdf.csv", index=False)
    # Per-site RDF plots
    for _site, _y in rdf_cols.items():
        plt.figure()
        plt.plot(r_g, _y, lw=1.8)
        plt.xlabel("r (Å)"); plt.ylabel("g(r)")
        plt.title(f"Na+ RDF: {_site}")
        plt.tight_layout()
        plt.savefig(Path(outdir) / f"Na+_rdf_{_site}.png", dpi=200)
        plt.close()

    # --------- Residence + CN + Solvation radii summary ---------
    # Use the density-based Solute object for Residence
    res = Residence.from_solute(sol_den)
    
    # Collect CN, radii, and residence times into one table
    cn_map     = getattr(sol_den.coordination, "coordination_numbers", {})
    radii_map  = getattr(sol_den, "radii", {})
    rt_cut_map = getattr(res, "residence_times_cutoff", {})
    rt_fit_map = getattr(res, "residence_times_fit", {})

    rows = []
    # Only include solvents that appear in density data (keeps alignment with CSV columns)
    for solvent in data_density.keys():
        rows.append({
            "solvent": solvent,
            "CN": float(cn_map.get(solvent, np.nan)),
            "solvation_radius_A": float(radii_map.get(solvent, np.nan)),
            "residence_time_cutoff": float(rt_cut_map.get(solvent, np.nan)),
            "residence_time_fit": float(rt_fit_map.get(solvent, np.nan)),
        })

    df_res = pd.DataFrame(rows, columns=["solvent","CN","solvation_radius_A","residence_time_cutoff","residence_time_fit"])
    df_res.to_csv(Path(outdir) / "Na+_solvation.csv", index=False)

    # Overall residence-time bar plot (fit preferred, else cutoff)
    y = df_res["residence_time_fit"].to_numpy()
    if np.all(np.isnan(y)):
        y = df_res["residence_time_cutoff"].to_numpy()
        ylabel = "Residence time (cutoff units)"
    else:
        ylabel = "Residence time (fit units)"
    plt.figure(figsize=(8, 4.5))
    plt.bar(df_res["solvent"], y)
    plt.ylabel(ylabel); plt.title("Na+ residence times by solvent site")
    plt.xticks(rotation=30, ha="right")
    plt.tight_layout()
    plt.savefig(Path(outdir) / "Na+_residence_times.png", dpi=200)
    plt.close()

    # --- Per-site residence C(t) line plots + CSV export ---
    # Each site gets: Na+_residence_<site>.png (C(t) line) and Na+_residenceCt_<site>.csv
    autocovariances = res.auto_covariances
    for site in df_res["solvent"].tolist():
        ct = autocovariances[site]
        t = np.arange(len(ct))  # frame index; adjust externally if you know dt
        # Save CSV
        pd.DataFrame({"frame": t, "C_t": np.asarray(ct, dtype=float)}).to_csv(Path(outdir) / f"Na+_residenceCt_{site}.csv", index=False)
        # Save plot
        plt.figure()
        plt.plot(t, ct, lw=1.8)
        plt.xlabel("Frame index"); plt.ylabel("C(t)")
        plt.title(f"Na+ residence C(t): {site}")
        plt.ylim([0,1])
        plt.tight_layout()
        plt.savefig(Path(outdir) / f"Na+_residence_{site}.png", dpi=200)
        plt.close()

    print(f"Done. Wrote Na+_rdf.csv, Na+_density.csv, Na+_residence.csv and per-site plots/CSVs to {args.outdir}")

if __name__ == "__main__":
    main()
