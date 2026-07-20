# copper-electrodeposition-data
This repository contains the computational data for the paper: "Direct Spectroscopic Evidence of an Atom-Transfer Mechanism during Copper Electrodeposition on Gold". The information is organized as follows:
1. The input and output files for the XAS calculated spectra are publicly available at Zenodo: [https://doi.org/10.5281/zenodo.17657023](https://doi.org/10.5281/zenodo.17657023)
2. The [structures](./structures/) folder contains the structures for the five kinds of systems studied:
* bulk-solids
* bulk-liquid-water
* molecules-at-au-interface
* solvated-ions
* layered-heterostructures
3. The [semiempiricalMD](./semiempiricalMD/) folder contains the trajectories of the dynamics using the PM6-FM method. A sample input [Au_111_Cu2O_atom_watbox.cp2k.in](./semiempiricalMD/Au_111_Cu2O_atom_watbox.cp2k.in) is also provided.
4. The [AIMD](./AIMD/) folder has inputs an outputs of ab initio molecular dynamics.
5. The [psp](./psp/) folder contains the pseudopotentials used, including the ones constructed with the excited core hole (XCH).
6. The [figures_data](./figures_data/) folder contains the data and python scripts used for plotting, organized in two folders. The table below ralates the figures of the manuscript with the scripts used for creating them.  

📁 `bulk/`

| Figure     | Related Scripts                                                                 |
|------------|----------------------------------------------------------------------------------|
| figS7-bulk | cu-l-edge/norm_area_david_paper.py<br>o-k-edge/norm_area_david_refBulkCuO_paper.py |
| figS14-bandDos | bandDos/bulkCu2O/plot_band_dos_combined.py<br>bandDos/bulkCuO/plot_band_dos_combined.py                                      |

---

📁 `interface/`

| Figure                               | Related Scripts                                                                                                       |
|--------------------------------------|------------------------------------------------------------------------------------------------------------------------|
| fig3-cu-vs-dist                      | cu-l-edge/au_cuso4_vs_distance_paper.py                                                                               |
| fig4-main-comp                       | cu-l-edge/norm_area_david_paper.py<br>o-k-edge/norm_area_david_refBulkCuO_refWater_paper.py                                     |
| figS8-compare0V                      | o-k-edge/compare0V_paper.py                                                                                                |
| figS9-GCS                            | o-k-edge/compare0V_paper.py                                                                                                |
| figS10-LCA                           | o-k-edge/LCA.py                                                                                                       |
| figS11-o-k-edge-interfacial-systems  | o-k-edge/norm_area_david_refBulkCuO_atom_paper.py<br>o-k-edge/norm_area_david_refBulkCuO_layer_solv_paper.py         |
| figS12-cu-l-edge-interfacial-systems | cu-l-edge/norm_area_david_layer_solv_paper.py<br>cu-l-edge/norm_area_david_atom_paper.py                               |
| figS13-compareIntermediary           | o-k-edge/compareIntermediarySignal.py<br>cu-l-edge/compareIntermediarySignal.py                                       |
