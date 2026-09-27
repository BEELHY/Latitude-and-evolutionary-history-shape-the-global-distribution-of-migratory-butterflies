# Linear richness GLMs with and without NRI/NTI.

suppressMessages({
  library(dplyr)
  library(readr)
  library(MASS)
  library(broom)
  library(spdep)
})

model_df <- read_csv("output/phylo_export/figure4d_phylo_sensitivity_model_data.csv", show_col_types = FALSE)
model_df$Landuse <- factor(model_df$Landuse)
model_df$Landuse <- relevel(model_df$Landuse, ref = "cropland")

cat(sprintf("Modelling dataset: %d grid cells\n", nrow(model_df)))

m_B_lin <- glm.nb(Richness_247 ~ Bio_4 + Bio_15 + Elevation + HII + Landuse, data = model_df)
m_C_lin <- glm.nb(Richness_247 ~ Bio_4 + Bio_15 + Elevation + HII + Landuse + NRI + NTI, data = model_df)

cat("\n=== B_lin: no NRI/NTI ===\n"); print(summary(m_B_lin))
cat("\n=== C_lin: with NRI/NTI ===\n"); print(summary(m_C_lin))

coords <- as.matrix(model_df[, c("lon1", "lat1")])
nb <- dnearneigh(coords, d1 = 0, d2 = 1.5)
listw <- nb2listw(nb, style = "W", zero.policy = TRUE)

moran_report <- function(mod, label) {
  mt <- moran.test(residuals(mod, type = "deviance"), listw, zero.policy = TRUE)
  cat(sprintf("\n[%s] Moran's I = %.3f, p %s\n", label,
              unname(mt$estimate["Moran I statistic"]),
              if (mt$p.value < 0.001) "< 0.001" else sprintf("= %.3f", mt$p.value)))
}
moran_report(m_B_lin, "B_lin (no NRI/NTI)")
moran_report(m_C_lin, "C_lin (with NRI/NTI)")

tidy_rr <- function(mod, label) {
  tidy(mod, conf.int = TRUE) %>%
    filter(term != "(Intercept)") %>%
    mutate(model = label,
           rate_ratio = exp(estimate),
           rr_lo = exp(conf.low), rr_hi = exp(conf.high))
}
comp <- bind_rows(tidy_rr(m_B_lin, "B_lin_no_NRI_NTI"), tidy_rr(m_C_lin, "C_lin_with_NRI_NTI"))
write_csv(comp, "output/phylo_export/linear_NRI_NTI_coefficient_comparison.csv")

cat("\n=== Coefficient comparison (rate ratios) ===\n")
print(comp %>% dplyr::select(model, term, estimate, std.error, p.value, rate_ratio, rr_lo, rr_hi), n = 30)

cat("\n=== AIC comparison ===\n")
cat(sprintf("B_lin AIC = %.1f | C_lin AIC = %.1f | delta = %.1f\n",
            AIC(m_B_lin), AIC(m_C_lin), AIC(m_B_lin) - AIC(m_C_lin)))

saveRDS(list(m_B_lin = m_B_lin, m_C_lin = m_C_lin, data = model_df, listw = listw),
        "output/phylo_export/linear_NRI_NTI_models.rds")
