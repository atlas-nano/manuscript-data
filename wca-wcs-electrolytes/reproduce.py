"""Regenerates every figure into figures/out/ by running the figure scripts in order."""
import argparse
import os
import pathlib
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parent
ORDER = ["fig1", "fig2", "fig3_panels_ae", "fig3f", "fig3", "fig4", "figS3", "figS5", "figS6", "figS8", "figS9",
         "figS11", "figS14", "figS15", "figS17", "figS19"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=str(ROOT / "figures" / "out"))
    ap.add_argument("--only", nargs="*")
    a = ap.parse_args()
    env = dict(os.environ, FIGURE_OUTDIR=str(pathlib.Path(a.outdir).resolve()), PYTHONDONTWRITEBYTECODE="1")
    for name in a.only or ORDER:
        r = subprocess.run([sys.executable, f"{name}.py"], cwd=ROOT / "figures", env=env)
        if r.returncode:
            sys.exit(f"{name}.py failed")


if __name__ == "__main__":
    main()
