# Analysis code

Code for the analyses and figures of the manuscript on the macroecology of migratory butterflies.

## Run order

| Script | Purpose | Main output |
|---|---|---|
| 01_range_tropics.R | Range size vs tropicality; phylogenetic models | (session objects) |
| 02_wingspan_traits.R | Join wingspan and range | output/trait_range.csv |
| 03_wingspan_latitude.R | Wingspan vs latitude; interaction models | output/df_lat_sp.csv |
| 04_phylo_bpmm.R | Climate-baseline Bayesian phylogenetic mixed models (H²) | (session objects) |
| 05_h2_exact_match.R | H² sensitivity, exact phylogeny matches only | output/h2_sensitivity/ |
| 06_attrition_tables.R | Sample attrition (Tables S1–S2) | output/TableS1_attrition.csv, output/TableS2_nested.csv |
| 07_hemisphere_sensitivity.R | Northern vs southern hemisphere sensitivity | console |
| 08_hpc_run.sh | Spatial richness GAM on HPC (calls 08_hpc_install.R, 08_richness_gam_hpc.R) | output/BAM/output/ |
| 09_richness_gam_figures.R | GAM diagnostic figures | output/BAM/figures/ |
| 10_nri_nti_data.R | NRI/NTI and covariates per equal-area cell (Behrmann, ~110 km) | output/phylo_export/ |
| 11_richness_gam_linear.R | Richness GAMs on equal-area cells (main model, 247 species, + NRI/NTI) | output/phylo_export/ |
| 12_forest_data.R | Fig. 2d data, Table S3 GAM columns, Fig. S5 (built-up reference) | output/phylo_export/, output/SI/ |
| 13_realm_area_s2.py | Equal-area richness and high-richness (≥10 species) area by zoogeographic realm; Fig. S2 | output/realm_area/ |
| 14_figure1_panels.py | Figure 1 (Robinson maps, latitude profile, seasonal net change) | output/Manuscript/figure1_panels/Figure1.png, .pdf |
| 15_figure2_panels.R | Figure 2 (and Fig. S3) | output/Manuscript/reproducibility_code/Figure2_final.png, .pdf |
| 17_figure3_panels.R, 18_figure3_compose.R | Figure 3 | output/Manuscript/reproducibility_code/ |
| 19_variance_partition.R | Integrated phylogenetic model and variance partitioning (Table S6, Fig. S9) | output/phylo_export/ |

- 01–04 run in one R session, in order (02–04 reuse objects created by 01).
- 08 needs a large-memory node (see SBATCH header); copy `output/BAM/output/` back before running 09.
- 13 and 14 are the only scripts that read the raw seasonal suitability maps (Chowdhury et al. 2025, *Conservation Biology*; method of Chowdhury et al. 2021, *Ecology Letters*), which are not redistributed here; they also need the CMEC zoogeographic regions (Holt et al. 2013) and a world boundaries shapefile. Their realm-level output is provided as `updatedata/realm_area_hotspot_summary.csv`.
- All other scripts run standalone once their inputs exist.
- Bayesian models use `seed = 1`; other random steps use `set.seed(1)`.

## Software

R 4.4+ with: ape, bayesplot, brms, broom, broom.mixed, car, datawizard, dplyr, effsize, GGally,
ggeffects, ggnewscale, ggplot2, ggtext, ggtree, ggtreeExtra, gratia, knitr, lme4, lmerTest, mgcv,
MuMIn, patchwork, performance, phytools, picante, png, purrr, ragg, readr, rsvg, scales, spdep, stringr, svglite, tibble,
tidybayes, tidyr, tidyverse, viridis.

Python 3 with: numpy, pandas, matplotlib, seaborn, geopandas, rasterio, cartopy, tqdm.
