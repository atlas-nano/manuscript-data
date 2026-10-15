"""Run every figure script in figures/ and write the PDFs to figures/out/."""
import os
import subprocess
import sys

HERE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
SCRIPTS = ["fig_lj_bias.py", "fig_metals_accuracy.py", "fig_cage_dos_partition.py",
           "fig_cage_efficiency.py", "fig_path_independence.py", "fig_SI_cage_sensitivity.py"]

for s in SCRIPTS:
    print(s)
    subprocess.run([sys.executable, s], cwd=HERE, check=True)
