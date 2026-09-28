# Donor-allocation effects vary with reference-cell budget in MuSiC PBMC simulations

Code and results accompanying the manuscript by Yunhao Jiang.

## Files

- `analysis/`: data preparation, fitting, summary and plotting scripts used in the study.
- `results/`: selected full-precision summaries, figure source data and displayed manuscript tables.
- `environment/`: recorded dependencies and MuSiC source version.
- `reproduction_integrated/`: parameterized code for the five initial PBMC simulation stages.

Processed inputs, simulated targets and truths, reference selections and saved predictions are in [Zenodo v1.2.0](https://doi.org/10.5281/zenodo.23012998). Download `data.zip` and extract it into this repository. Public source accessions and software citations are in [SOURCES.md](SOURCES.md); original GEO reads and installed software are not redistributed.

## Analyses

| Analysis | Entry points | Saved results |
|---|---|---|
| Initial PBMC simulations | `reproduction_integrated/code/run.py`, `summarize.py` | `results/music_*_review/` in the data archive |
| Additional reference draws | `analysis/run_music_mc_extension_v1.R`, `review_music_mc_extension_v1.py` | `results/music_mc_extension_review/` |
| Component interventions | `analysis/run_music_components.R`, `review_music_components.py` | `results/music_components/` |
| Exploratory pancreas extension | `analysis/run_music_pancreas.R`, `review_music_pancreas.py` | `results/music_pancreas/` |
| Measured bulk and input diagnostics | `analysis/run_real_bulk*.R`, `run_salmon_pilot_music.R` | `results/real_bulk*` |

See [run_commands.md](run_commands.md) for arguments and plotting commands. The repository retains necessary source dependencies; older exploratory work and historical releases remain available in v1.1.0.

The pancreas primary comparison was not reproduced, and measured-bulk assessment did not validate transfer of the PBMC trend. Their saved results are included. Component interventions provide computational attribution, not a biological causal mechanism. Monte Carlo blocks quantify conditional sampling variation, not population uncertainty.

## Versions and licence

v1.2.0 adds the pancreas and component analyses and updates displayed tables and figure sources without rerunning the scientific analyses. The existing [v1.1.0 release](https://github.com/Re1nhard-1/Donor-allocation-effects-vary-with-reference-cell-budget-in-MuSiC-PBMC-simulations/releases/tag/v1.1.0) and [data DOI](https://doi.org/10.5281/zenodo.22902701) remain available for earlier manuscript versions. The preprint is [Research Square v1](https://doi.org/10.21203/rs.3.rs-11156252/v1).

Project code is GPL-3.0-or-later; project results and documentation are CC BY 4.0. Third-party material retains its original terms; see LICENSE.md and SOURCES.md.
