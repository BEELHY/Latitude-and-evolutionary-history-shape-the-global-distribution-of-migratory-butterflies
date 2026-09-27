# Richness GAMs with NRI/NTI phylogenetic sensitivity.

suppressMessages({
  library(terra)
  library(ape)
  library(picante)
  library(mgcv)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
})

set.seed(1)
dir.create("output/phylo_export", showWarnings = FALSE, recursive = TRUE)
dir.create("output/tmp", showWarnings = FALSE, recursive = TRUE)

tree <- readRDS("output/checkpoints/tree_final_fixed.rds")
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

AGG_FACT <- 24
raster_dir <- "data/SuitabilityMaps_MigratorySpecies"
files <- list.files(raster_dir, pattern = "\\.tif$", full.names = TRUE)
file_sp <- sub("^Binary_S[0-9]+", "", sub("\\.tif$", "", basename(files)))

files_tree <- files[file_sp %in% tree_sp]
sp_of_file <- file_sp[file_sp %in% tree_sp]
sp_list <- sort(unique(sp_of_file))
cat(sprintf("Rasters to aggregate (tree species x season): %d for %d species\n",
            length(files_tree), length(sp_list)))

ref_template <- rast(xmin = -180, xmax = 180, ymin = -60, ymax = 85,
                      resolution = 1, crs = "EPSG:4326")
n_cell <- ncell(ref_template)
cell_xy <- xyFromCell(ref_template, 1:n_cell)
colnames(cell_xy) <- c("lon1", "lat1")

agg_presence <- function(f) {
  r <- rast(f)
  r_agg <- terra::aggregate(r, fact = AGG_FACT, fun = "max", na.rm = TRUE)
  r_ext <- terra::extend(r_agg, ref_template)
  r_ext <- terra::resample(r_ext, ref_template, method = "near")
  v <- terra::values(r_ext)[, 1]
  v[is.na(v)] <- 0
  as.integer(v > 0)
}

comm_tree <- matrix(0L, nrow = n_cell, ncol = length(sp_list),
                     dimnames = list(NULL, sp_list))

cat("Aggregating rasters to 1 degree (progress every 25 species)...\n")
t0 <- Sys.time()
for (i in seq_along(sp_list)) {
  sp <- sp_list[i]
  sp_files <- files_tree[sp_of_file == sp]
  sp_presence <- rep(0L, n_cell)
  for (f in sp_files) sp_presence <- pmax(sp_presence, agg_presence(f))
  comm_tree[, sp] <- sp_presence
  if (i %% 25 == 0 || i == length(sp_list)) {
    cat(sprintf("  [%d/%d] %s | elapsed %.1f min\n", i, length(sp_list), sp,
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
}

richness_247 <- rowSums(comm_tree)
saveRDS(list(comm_tree = comm_tree, cell_xy = cell_xy, richness_247 = richness_247),
        "output/tmp/comm_tree_1deg.rds")
cat(sprintf("Cells with >=1 tree-covered species present: %d / %d\n",
            sum(richness_247 > 0), n_cell))

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
  lon1 = cell_xy_sub[, "lon1"],
  lat1 = cell_xy_sub[, "lat1"],
  Richness_247 = richness_247[keep],
  NRI = -1 * mpd_res$mpd.obs.z,
  NTI = -1 * mntd_res$mntd.obs.z
)
write_csv(nri_nti, "output/phylo_export/figure4d_NRI_NTI_by_cell.csv")

cat("Loading + aggregating cleaned_data_for_SAR.csv to 1 degree...\n")
env_fine <- read_csv("output/BAM/output/cleaned_data_for_SAR.csv", show_col_types = FALSE)

env_fine <- env_fine %>%
  mutate(lon1 = floor(lon) + 0.5, lat1 = floor(lat) + 0.5)

landuse_mode <- function(x) names(sort(table(x), decreasing = TRUE))[1]

env_coarse <- env_fine %>%
  group_by(lon1, lat1) %>%
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
  inner_join(env_coarse, by = c("lon1", "lat1")) %>%
  mutate(Landuse = factor(Landuse)) %>%
  filter(!is.na(Bio_4), !is.na(Bio_15), !is.na(Elevation), !is.na(HII))
model_df$Landuse <- relevel(model_df$Landuse, ref = "cropland")

cat(sprintf("Final modeling dataset: %d grid cells (1-degree)\n", nrow(model_df)))
write_csv(model_df, "output/phylo_export/figure4d_phylo_sensitivity_model_data.csv")

fit_gam <- function(formula, data) {
  bam(formula, data = data, family = nb(), method = "fREML")
}

form_base <- function(resp) {
  as.formula(sprintf(
    "%s ~ s(Bio_4) + s(Bio_15) + s(Elevation) + HII + Landuse + s(lon1, lat1, k=100)",
    resp))
}

cat("Fitting model A (full species pool, no phylogeny)...\n")
m_full   <- fit_gam(form_base("Richness_full"), model_df)

cat("Fitting model B (247-species subset, no phylogeny)...\n")
m_subset <- fit_gam(form_base("Richness_247"), model_df)

cat("Fitting model C (247-species subset + NRI + NTI)...\n")
m_phylo  <- fit_gam(
  Richness_247 ~ s(Bio_4) + s(Bio_15) + s(Elevation) +
    HII + Landuse + NRI + NTI + s(lon1, lat1, k=100),
  data = model_df)

saveRDS(list(m_full = m_full, m_subset = m_subset, m_phylo = m_phylo, data = model_df),
        "output/phylo_export/figure4d_phylo_sensitivity_models.rds")

extract_terms <- function(model, label) {
  s <- summary(model)
  para <- as.data.frame(s$p.table) %>% rownames_to_column("term") %>% mutate(part = "parametric")
  smooth <- as.data.frame(s$s.table) %>% rownames_to_column("term") %>% mutate(part = "smooth")
  bind_rows(para, smooth) %>% mutate(model = label)
}

comparison <- bind_rows(
  extract_terms(m_full,   "A_full_400species_no_phylo"),
  extract_terms(m_subset, "B_subset_247species_no_phylo"),
  extract_terms(m_phylo,  "C_subset_247species_plus_NRI_NTI")
)
write_csv(comparison, "output/phylo_export/figure4d_phylo_sensitivity_term_comparison.csv")

cat("\n=== Sensitivity test complete ===\n")
cat("Outputs:\n")
cat("  output/phylo_export/figure4d_NRI_NTI_by_cell.csv\n")
cat("  output/phylo_export/figure4d_phylo_sensitivity_model_data.csv\n")
cat("  output/phylo_export/figure4d_phylo_sensitivity_models.rds\n")
cat("  output/phylo_export/figure4d_phylo_sensitivity_term_comparison.csv\n")
print(comparison)
