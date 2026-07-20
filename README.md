# Manuscript data

Reproducibility data and analysis scripts for manuscripts from the ATLAS
Materials Physics Laboratory (UC San Diego), organized one subdirectory per
paper.

| Subdirectory | Paper / topic |
|---|---|
| [`3pt-construction/`](3pt-construction/) | *An anharmonic liquid-entropy functional from the Mori–Zwanzig memory kernel* — data, figure/table scripts, and control-file templates for the Three-Phase Explicit Anharmonic Thermodynamics (3PT) method. |
| [`solvation-spectra/`](solvation-spectra/) | *Practical Corrections for Finite-Concentration Molecular Dynamics Simulations* — inputs and analysis scripts (Zenodo [10.5281/zenodo.15377288](https://doi.org/10.5281/zenodo.15377288)). |
| [`na-orion-conductor/`](na-orion-conductor/) | *Tailoring Energy Landscapes for Vehicular Transport in Single-Ion Conducting Organo Ionic Solids* — molecular-dynamics data. |
| [`copper-electrodeposition-data/`](copper-electrodeposition-data/) | *Direct Spectroscopic Evidence of an Atom-Transfer Mechanism during Copper Electrodeposition on Gold* — computational data (XAS spectra on Zenodo [10.5281/zenodo.17657023](https://doi.org/10.5281/zenodo.17657023)). |
| [`sxshg-air-water/`](sxshg-air-water/) | Supporting information for the second-harmonic-generation (SXSHG) study of the air/water interface. |
| [`interfacial-carbonates/`](interfacial-carbonates/) | Supporting information for the interfacial carbonates manuscript (JACS). |
| [`lithiated-graphite-xas/`](lithiated-graphite-xas/) | VASP inputs and structures for molecular-dynamics / XAS of lithiated graphite (LiC_x). |
| [`llto-xray-shg/`](llto-xray-shg/) | Real-time VG-RTTDFT SHG simulation inputs for LLTO (X-ray SHG). |

Large molecular-dynamics trajectories are excluded for size where noted; each
subdirectory ships the post-processed analysis data and control templates needed
to regenerate its figures. Where a subdirectory also has a Zenodo archive, the
DOI above is the citable snapshot.

Companion research codes (py-xPT, lj-4d-md, DMAx) live in the
[`atlas-nano/codes`](https://github.com/atlas-nano/codes) repository.
