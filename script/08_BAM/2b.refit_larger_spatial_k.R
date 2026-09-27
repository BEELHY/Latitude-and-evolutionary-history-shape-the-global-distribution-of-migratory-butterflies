# Refit sensitivity GAMs with larger spatial basis.

suppressMessages({
  library(mgcv)
  library(dplyr)
  library(readr)
  library(tibble)
})

model_df <- read_csv("output/phylo_export/figure4d_phylo_sensitivity_model_data.csv", show_col_types = FALSE)
model_df$Landuse <- factor(model_df$Landuse)
model_df$Landuse <- relevel(model_df$Landuse, ref = "cropland")

cat(sprintf("Modeling dataset: %d grid cells\n", nrow(model_df)))

SPATIAL_K <- 800

fit_gam <- function(formula, data) {
  bam(formula, data = data, family = nb(), method = "fREML")
}

form_base <- function(resp) {
  as.formula(sprintf(
    "%s ~ s(Bio_4) + s(Bio_15) + s(Elevation) + HII + Landuse + s(lon1, lat1, k=%d)",
    resp, SPATIAL_K))
}

cat(sprintf("Fitting model A (full species pool, no phylogeny), spatial k=%d...\n", SPATIAL_K))
m_full   <- fit_gam(form_base("Richness_full"), model_df)

cat(sprintf("Fitting model B (247-species subset, no phylogeny), spatial k=%d...\n", SPATIAL_K))
m_subset <- fit_gam(form_base("Richness_247"), model_df)

cat(sprintf("Fitting model C (247-species subset + NRI + NTI), spatial k=%d...\n", SPATIAL_K))
m_phylo  <- fit_gam(
  as.formula(sprintf(
    "Richness_247 ~ s(Bio_4) + s(Bio_15) + s(Elevation) + HII + Landuse + NRI + NTI + s(lon1, lat1, k=%d)",
    SPATIAL_K)),
  data = model_df)

for (mod in list(A = m_full, B = m_subset, C = m_phylo)) {
  sp_edf <- summary(mod)$s.table["s(lon1,lat1)", "edf"]
  cat(sprintf("  spatial smooth edf = %.1f (k = %d)\n", sp_edf, SPATIAL_K))
}

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

cat("\n=== Refit with larger spatial k complete ===\n")
print(comparison)
