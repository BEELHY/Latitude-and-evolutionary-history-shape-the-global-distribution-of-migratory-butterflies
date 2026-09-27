# Analysis code

Code for the analyses and figures of the manuscript on the macroecology of migratory butterflies.
All scripts read secondary data from `updatedata/` (available from the authors / data repository)
or files written by earlier scripts into `output/`. Run everything from the project root.

## Run order

| Script | Purpose | Main output |
|---|---|---|
| 01_range_tropics.R | Range size vs tropicality; phylogenetic models | 
| 02_wingspan_traits.R | Join wingspan and range |
| 03_wingspan_latitude.R | Wingspan vs latitude; interaction models | 
| 04_phylo_bpmm.R | Climate-baseline Bayesian phylogenetic mixed models (H²) | 
| 05_h2_exact_match.R | H² sensitivity, exact phylogeny matches only |
| 06_attrition_tables.R | Sample attrition (Tables S1–S2) | 
| 07_hemisphere_sensitivity.R | Northern vs southern hemisphere sensitivity | 
| 08_hpc_run.sh | Spatial richness GAM on HPC (calls 08_hpc_install.R, 08_richness_gam_hpc.R) | 
| 09_richness_gam_figures.R | GAM diagnostic figures | 
| 10_nri_nti_data.R | NRI/NTI per 1° cell and model data | 
| 11_richness_gam_linear.R | Linear-covariate GAMs with spatial smooth | 
| 12_forest_data.R | Forest plot data, built-up reference | 
| 13_region_hotspots.R | Hotspot metrics by region | 
| 14_realm_hotspots.R | Hotspot metrics by zoogeographic realm | 
| 15_figure2_panels.R, 16_figure2_compose.py | Figure 2 | 
| 17_figure3_panels.R, 18_figure3_compose.py | Figure 3 | 

- 01–04 run in one R session, in order (02–04 reuse objects created by 01).
- 08 needs a large-memory node (see SBATCH header); copy `output/BAM/output/` back before running 09.
- All other scripts run standalone once their inputs exist.
- Bayesian models use `seed = 1`; other random steps use `set.seed(1)`.

## Software

R 4.4+ with: ape, bayesplot, brms, broom, broom.mixed, car, datawizard, dplyr, effsize, GGally,
ggeffects, ggimage, ggplot2, ggtree, ggtreeExtra, gratia, knitr, lme4, lmerTest, mgcv,
MuMIn, patchwork, performance, phytools, picante, png, purrr, ragg, readr, spdep, stringr, tibble,
tidybayes, tidyr, tidyverse, viridis.

Python 3 with: Pillow, numpy.
