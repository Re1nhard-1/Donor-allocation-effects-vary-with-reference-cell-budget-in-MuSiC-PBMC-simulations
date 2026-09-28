# Run commands

Run from the repository root after extracting data.zip. Python packages are recorded in environment/requirements-deconvolution.lock.txt and requirements-source-tools.txt. R scripts use environment/music/runtime.R and a local .tools/R-library with the recorded MuSiC version. The integrated runner accepts explicit Rscript, library and upstream-source paths. Source scripts retain their original study paths; run them in a separate working copy to avoid overwriting saved results.

## Five initial PBMC stages

```
python reproduction_integrated/code/run.py --help
python reproduction_integrated/code/summarize.py --endpoints work/predictions/endpoints --gradient work/predictions/gradient --thinning work/predictions/thinning --external work/predictions/external --alternative work/predictions/alternative --output recomputed
```

## Component interventions

The public PBMC data layout stores donor matrices and cells.tsv in work/inputs/external. Copy counts_*.mtx.gz and cells.tsv from that directory into data/processed/music_external/20260917T122709190802Z before invoking the original source scripts. The frozen cases, reference indices and targets are in results/music_mechanism_selection/20260917T205706431200Z. Replace TRIPLE below with each recorded triple_key in references.tsv.

```
Rscript analysis/run_music_components.R results/music_mechanism_selection/20260917T205706431200Z/selection.json results/music_components TRIPLE full
python analysis/review_music_components.py
python analysis/plot_music_components.py
```

The saved predictions.csv supports the Python summaries directly. If refitting all cases in a clean output directory, analysis/export_music_components.R converts the newly generated RDS objects into predictions.csv. Intermediate RDS objects are not bundled. The frozen intervention protocol is analysis/MUSIC_COMPONENT_INTERVENTION_PROTOCOL.md.

## Pancreas

Prepared matrices, gene/cell annotations, target truths and reference selections are in data/processed/music_pancreas. For a full refit, run the following for each of human1, human2, human3 and human4:

```
Rscript analysis/run_music_pancreas.R human1 full
python analysis/review_music_pancreas.py --mode full
```

analysis/prepare_music_pancreas.py documents input construction. To rebuild rather than use the supplied matrices, use a fresh working copy without data/processed/music_pancreas, obtain the four original GEO count files at the URLs and raw_path locations in results/deconv_input/20260916T210631005068Z/input_manifest.json, and retain the supplied processed Baron inputs. The original parser is analysis/prepare_baron.py; it writes a new timestamped run.

## Figures and further analyses

```
python analysis/plot_manuscript_figures.py
python analysis/plot_profiles_from_source.py
python analysis/plot_real_bulk.py
python analysis/summarize_salmon_pilot.py --plot-only
python analysis/plot_music_practical_effects.py
```

The pancreas review produces Fig. S7; the practical-effects plot produces Fig. S8; the component plot produces Fig. S9. The main figure script writes into the figures directory beside the repository. Additional-draw and measured-bulk scripts retain their recorded parameters and input paths in their accompanying analysis protocols. Raw read quantification requires the separate Salmon environment and original public sequencing data.

These are source entry points and existing saved outputs. Release assembly checks file integrity and required paths; it does not claim a fresh end-to-end scientific rerun or independent reproduction.
