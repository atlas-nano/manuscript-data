# Manuscript data

Reproducibility data and analysis scripts for manuscripts from the ATLAS
Materials Physics Laboratory (UC San Diego), organized one subdirectory per
paper.

| Subdirectory | Paper / topic |
|---|---|
| [`3pt-construction/`](3pt-construction/) | *An anharmonic liquid-entropy functional from the Mori–Zwanzig memory kernel* — data, figure/table scripts, and control-file templates for the Three-Phase Explicit Anharmonic Thermodynamics (3PT) method (Zenodo [10.5281/zenodo.21447740](https://doi.org/10.5281/zenodo.21447740)). |
| [`solvation-spectra/`](solvation-spectra/) | *Practical Corrections for Finite-Concentration Molecular Dynamics Simulations* — inputs and analysis scripts (Zenodo [10.5281/zenodo.15377288](https://doi.org/10.5281/zenodo.15377288)). |
| [`na-orion-conductor/`](na-orion-conductor/) | *Tailoring Energy Landscapes for Vehicular Transport in Single-Ion Conducting Organo Ionic Solids* (Zellmann-Parrotta *et al.*, JACS, in press) — molecular-dynamics data (Zenodo [10.5281/zenodo.21448069](https://doi.org/10.5281/zenodo.21448069)). |
| [`copper-electrodeposition-data/`](copper-electrodeposition-data/) | *Direct Spectroscopic Evidence of an Atom-Transfer Mechanism during Copper Electrodeposition on Gold* — computational data (XAS spectra on Zenodo [10.5281/zenodo.17657023](https://doi.org/10.5281/zenodo.17657023)). |
| [`sxshg-air-water/`](sxshg-air-water/) | *Surface structure of water from soft X-ray second harmonic generation* (Hoffman *et al.*, Nat. Commun. **16**, 2025, [10.1038/s41467-025-65514-4](https://doi.org/10.1038/s41467-025-65514-4)) — supporting information (Zenodo [10.5281/zenodo.17230045](https://doi.org/10.5281/zenodo.17230045)). |
| [`interfacial-carbonates/`](interfacial-carbonates/) | *Agglomeration Drives the Reversed Fractionation of Aqueous Carbonate and Bicarbonate at the Air–Water Interface* (Devlin *et al.*, JACS **145**, 2023, [10.1021/jacs.3c05093](https://doi.org/10.1021/jacs.3c05093)) — supporting information (Zenodo [10.5281/zenodo.8327142](https://doi.org/10.5281/zenodo.8327142)). |
| [`lithiated-graphite-xas/`](lithiated-graphite-xas/) | *Electronic signatures of Lorentzian dynamics and charge fluctuations in lithiated graphite structures* (Jamnuch *et al.*, Nat. Commun. **14**, 2023, [10.1038/s41467-023-37857-3](https://doi.org/10.1038/s41467-023-37857-3)) — VASP inputs / XAS data (Zenodo [10.5281/zenodo.7739094](https://doi.org/10.5281/zenodo.7739094)). |
| [`llto-xray-shg/`](llto-xray-shg/) | *Probing lithium mobility at a solid electrolyte surface* (Woodahl *et al.*, Nat. Mater. **22**, 2023, [10.1038/s41563-023-01535-y](https://doi.org/10.1038/s41563-023-01535-y)) — real-time VG-RTTDFT SHG inputs for LLTO (Dryad [10.6078/D1N41X](https://doi.org/10.6078/D1N41X)). |

Large molecular-dynamics trajectories are excluded for size where noted; each
subdirectory ships the post-processed analysis data and control templates needed
to regenerate its figures. Where a subdirectory also has a data archive (Zenodo
or Dryad), that DOI is the citable snapshot.

Companion research codes (py-xPT, lj-4d-md, DMAx) live in the
[`atlas-nano/codes`](https://github.com/atlas-nano/codes) repository.
