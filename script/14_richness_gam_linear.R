# Linear GAM, NRI/NTI appendix.

suppressMessages({
  library(dplyr)
  library(readr)
  library(mgcv)
  library(broom)
  library(spdep)
})

model_df <- read_csv("output/phylo_export/figure4d_phylo_sensitivity_model_data.csv", show_col_types = FALSE)
model_df$Landuse <- factor(model_df$Landuse)
model_df$Landuse <- relevel(model_df$Landuse, ref = "cropland")
cat(sprintf("Modelling dataset: %d grid cells\n", nrow(model_df)))

SPATIAL_K <- 800
fit_gam <- function(formula, data) bam(formula, data = data, family = nb(), method = "fREML")

coords <- as.matrix(model_df[, c("lon1", "lat1")])
nb <- dnearneigh(coords, d1 = 0, d2 = 1.5)
listw <- nb2listw(nb, style = "W", zero.policy = TRUE)

moran_report <- function(mod, label) {
  mt <- moran.test(residuals(mod, type = "deviance"), listw, zero.policy = TRUE)
  p_str <- if (mt$p.value < 0.001) "< 0.001" else sprintf("= %.3f", mt$p.value)
  cat(sprintf("[%s] Moran's I = %.3f, p %s\n", label, unname(mt$estimate["Moran I statistic"]), p_str))
}

# Main: global richness
cat("\n=== MAIN: A_lin_sm, global richness (Richness_full) ===\n")
m_A <- fit_gam(
  Richness_full ~ Bio_4 + Bio_15 + Elevation + HII + Landuse + s(lon1, lat1, k = SPATIAL_K),
  model_df)
print(summary(m_A))
moran_report(m_A, "A_lin_sm (global, main)")

# Appendix: 247 species, NRI/NTI
cat("\n=== SUPP: B_lin_sm, 247-sp subset, no NRI/NTI ===\n")
m_B <- fit_gam(
  Richness_247 ~ Bio_4 + Bio_15 + Elevation + HII + Landuse + s(lon1, lat1, k = SPATIAL_K),
  model_df)
print(summary(m_B))
moran_report(m_B, "B_lin_sm (247-sp, no NRI/NTI)")

cat("\n=== SUPP: C_lin_sm, 247-sp subset, with NRI/NTI ===\n")
m_C <- fit_gam(
  Richness_247 ~ Bio_4 + Bio_15 + Elevation + HII + Landuse + NRI + NTI + s(lon1, lat1, k = SPATIAL_K),
  model_df)
print(summary(m_C))
moran_report(m_C, "C_lin_sm (247-sp, with NRI/NTI)")

cat(sprintf("\nAIC: A=%.1f | B=%.1f | C=%.1f\n", AIC(m_A), AIC(m_B), AIC(m_C)))

saveRDS(list(m_A = m_A, m_B = m_B, m_C = m_C, data = model_df, listw = listw),
        "output/phylo_export/linear_covariates_with_spatial_smooth_models.rds")

extract_terms <- function(model, label) {
  s <- summary(model)
  para <- as.data.frame(s$p.table); para$term <- rownames(para); para$part <- "parametric"
  smooth <- as.data.frame(s$s.table); smooth$term <- rownames(smooth); smooth$part <- "smooth"
  rbind(para[, c("term","Estimate","Std. Error","Pr(>|z|)","part")] %>%
          setNames(c("term","estimate","std_error","p_value","part")),
        data.frame(term = smooth$term, estimate = NA, std_error = NA,
                   p_value = smooth[["p-value"]], part = "smooth")) %>%
    mutate(model = label)
}
comparison <- bind_rows(extract_terms(m_A, "A_global_main"),
                         extract_terms(m_B, "B_247sp_no_NRI_NTI"),
                         extract_terms(m_C, "C_247sp_with_NRI_NTI"))
write_csv(comparison, "output/phylo_export/linear_covariates_with_spatial_smooth_comparison.csv")
cat("\nSaved comparison table.\n")
