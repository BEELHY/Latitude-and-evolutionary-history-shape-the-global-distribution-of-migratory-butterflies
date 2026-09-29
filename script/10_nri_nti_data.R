# NRI/NTI and covariates on an equal-area grid (Behrmann, ~110 km cells).

suppressMessages({
  library(ape)
  library(picante)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(terra)
})

set.seed(1)
dir.create("output/phylo_export", showWarnings = FALSE, recursive = TRUE)

# Tree, drop duplicate tip
tree <- read.tree("updatedata/phylogeny_matched.tre")
dup_idx <- which(duplicated(tree$tip.label))
if (length(dup_idx) > 0) {
  cat(sprintf("Dropping %d duplicated genus-proxy tip(s): %s\n",
              length(dup_idx), paste(tree$tip.label[dup_idx], collapse = ", ")))
  tree <- drop.tip(tree, dup_idx)
}
stopifnot(!any(duplicated(tree$tip.label)))
tree_sp <- tree$tip.label
cat(sprintf("Tree species after dedupe: %d\n", length(tree_sp)))

phylo_dist_full <- cophenetic(tree)

# Presence in equal-area cells (species present in any season)
comm <- readRDS("updatedata/community_ea.rds")
comm_tree    <- comm$comm_tree
cell_xy      <- comm$cell_xy
richness_247 <- comm$richness_247
cat(sprintf("Cells with >=1 tree-covered species present: %d / %d\n",
            sum(richness_247 > 0), nrow(cell_xy)))

# NRI/NTI, cells with >=3 species
MIN_RICHNESS_FOR_NRI <- 3
keep <- richness_247 >= MIN_RICHNESS_FOR_NRI
comm_sub <- comm_tree[keep, , drop = FALSE]
cell_xy_sub <- cell_xy[keep, , drop = FALSE]
cat(sprintf("Cells with >=%d tree-covered species (used for NRI/NTI): %d\n",
            MIN_RICHNESS_FOR_NRI, sum(keep)))

phylo_dist <- phylo_dist_full[colnames(comm_sub), colnames(comm_sub)]

cat("Running ses.mpd / ses.mntd (picante), null.model = 'richness', runs = 199...\n")
t0 <- Sys.time()
mpd_res  <- ses.mpd(comm_sub, phylo_dist, null.model = "richness",
                     abundance.weighted = FALSE, runs = 199)
mntd_res <- ses.mntd(comm_sub, phylo_dist, null.model = "richness",
                      abundance.weighted = FALSE, runs = 199)
cat(sprintf("Done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

nri_nti <- data.frame(
  ix = cell_xy_sub[, "ix"], iy = cell_xy_sub[, "iy"],
  x_km = cell_xy_sub[, "x_km"], y_km = cell_xy_sub[, "y_km"],
  lon = cell_xy_sub[, "lon"], lat = cell_xy_sub[, "lat"],
  Richness_247 = richness_247[keep],
  NRI = -1 * mpd_res$mpd.obs.z,
  NTI = -1 * mntd_res$mntd.obs.z
)
write_csv(nri_nti, "output/phylo_export/figure4d_NRI_NTI_by_cell.csv")

# Covariates aggregated to the same equal-area cells
cat("Loading + aggregating richness_grid.csv to equal-area cells...\n")
env_fine <- read_csv("updatedata/richness_grid.csv", show_col_types = FALSE)

HALF_W <- 17367530.445161372  # Behrmann x at 180 degrees
CELL <- 2 * HALF_W / 316        # 316 columns: ~109.9 km square cells
xy <- project(cbind(env_fine$lon, env_fine$lat), from = "EPSG:4326", to = "ESRI:54017")
env_fine <- env_fine %>%
  mutate(ix = floor((xy[, 1] + HALF_W) / CELL), iy = floor(xy[, 2] / CELL))

landuse_mode <- function(x) names(sort(table(x), decreasing = TRUE))[1]

env_coarse <- env_fine %>%
  group_by(ix, iy) %>%
  summarise(
    Richness_full = sum(Richness, na.rm = TRUE),
    Bio_4         = mean(Bio_4, na.rm = TRUE),
    Bio_15        = mean(Bio_15, na.rm = TRUE),
    Elevation     = mean(Elevation, na.rm = TRUE),
    HII           = mean(HII, na.rm = TRUE),
    Landuse       = landuse_mode(Landuse),
    n_fine_cells  = n(),
    .groups = "drop"
  )

model_df <- nri_nti %>%
  inner_join(env_coarse, by = c("ix", "iy")) %>%
  mutate(Landuse = factor(Landuse)) %>%
  filter(!is.na(Bio_4), !is.na(Bio_15), !is.na(Elevation), !is.na(HII))
model_df$Landuse <- relevel(model_df$Landuse, ref = "cropland")

cat(sprintf("Final modeling dataset: %d equal-area grid cells\n", nrow(model_df)))
write_csv(model_df, "output/phylo_export/figure4d_phylo_sensitivity_model_data.csv")
