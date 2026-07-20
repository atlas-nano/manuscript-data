# Manuscript data

Reproducibility data and analysis scripts for manuscripts from the ATLAS
Materials Physics Laboratory, organized one subdirectory per paper.

- [`3pt-construction/`](3pt-construction/) — data, figure/table scripts, and
  control-file templates for the Three-Phase Explicit Anharmonic Thermodynamics
  (3PT) entropy paper, *"An anharmonic liquid-entropy functional from the
  Mori–Zwanzig memory kernel."* See its README for the figure/table → script map
  and `reproduce.py`.

Molecular-dynamics trajectories are excluded for size; each subdirectory ships
the post-processed analysis data needed to regenerate its figures, plus LAMMPS
and analysis control templates. Entropy analysis uses the companion
[py-xPT](https://github.com/atlas-nano/codes/tree/main/py-xPT) code.
