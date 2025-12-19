Master Thesis Repository

From songbirds to insect chirps: acoustic indices capture diel and seasonal dynamics across Dutch habitats

This repository contains all data files and R scripts used for the analyses presented in the Master Thesis “From songbirds to insect chirps: acoustic indices capture diel and seasonal dynamics across Dutch habitats”. The materials enable transparency, reproducibility, and reuse of the analytical workflow.

Repository Structure and Contents

R Scripts:

Annotated_analysis.R
R script used for the analysis of soundscape composition based on annotated audio data.

Indices analysis.R
Full R script used to analyse the outputs of the acoustic indices, including data processing, statistical analyses, and visualisation.

Indices were computed by Dr. C. Bernard. Python code can be found via the following link:
https://github.com/BenMcEwen1/TABMON-Classification-Pipeline/blob/main/acoustic_indices/compute_acoustic_indices.py

Data Files:

annotated_audio_files_extra.csv
Raw annotated soundscape data used in the soundscape composition analyses.

glmmTMB_summary.xlsx
Summary output of the generalized linear mixed-effects models (GLMMs) fitted using the glmmTMB package.

acoustic_indices_bugg_1000000*.csv
Raw acoustic index output files.

Each file name includes:

The Sensor ID

The time period over which acoustic indices were calculated

Additional Outputs:

acoustic_indices_bugg/
Directory containing repeated acoustic index outputs. These files are included to compensate for missing original outputs.

Notes:

All analyses were conducted in R; package versions and dependencies should be checked within the scripts if full reproducibility is required.

Citation:

If you use any data or code from this repository, please cite the corresponding Master Thesis.
