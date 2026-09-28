# Analysis code

Code for the analyses and figures of the manuscript on the macroecology of migratory butterflies.
## Run order

| Script | Purpose | 
|---|---|
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
| 13_realm_area_s2.py | Equal-area high-richness (≥10 species) area by zoogeographic realm; Fig. S2 | 
| 14_figure1_panels.py | Figure 1 panels (Robinson maps, latitude profile, colour bars) | 
| 15_figure2_panels.R, 16_figure2_compose.py | Figure 2 | 
| 17_figure3_panels.R, 18_figure3_compose.py | Figure 3 | 
| 19_variance_partition.R | Integrated model and variance partitioning | 

- 01–04 run in one R session, in order (02–04 reuse objects created by 01).
- 08 needs a large-memory node (see SBATCH header); copy `output/BAM/output/` back before running 09.
- 13 and 14 read the raw seasonal suitability maps (Chowdhury et al. 2021), which are not redistributed; their realm-level output is provided in `updatedata/`.
- All other scripts run standalone once their inputs exist.
- Bayesian models use `seed = 1`; other random steps use `set.seed(1)`.

## Software

R 4.4+ with: ape, bayesplot, brms, broom, broom.mixed, car, datawizard, dplyr, effsize, GGally,
ggeffects, ggimage, ggplot2, ggtree, ggtreeExtra, gratia, knitr, lme4, lmerTest, mgcv,
MuMIn, patchwork, performance, phytools, picante, png, purrr, ragg, readr, spdep, stringr, tibble,
tidybayes, tidyr, tidyverse, viridis.

Python 3 with: Pillow, numpy, pandas, matplotlib, seaborn, geopandas, rasterio, cartopy, tqdm.
