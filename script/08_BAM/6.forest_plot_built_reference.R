# Forest plot of main GAM, built-up reference.

suppressMessages({
  library(mgcv)
  library(dplyr)
  library(ggplot2)
  library(ragg)
})

bundle <- readRDS("output/phylo_export/linear_covariates_with_spatial_smooth_models.rds")
mod <- bundle$m_A

cf <- coef(mod)
V  <- vcov(mod)

env_vars <- c(Bio_4 = "Bio_4", Bio_15 = "Bio_15", Elevation = "Elevation", HII = "HII")
env_df <- tibble::tibble(
  variable = names(env_vars),
  estimate = cf[unname(env_vars)],
  se       = sqrt(diag(V)[unname(env_vars)]),
  type     = "Environmental"
)

b_built <- cf["Landusebuilt"]; v_built <- V["Landusebuilt", "Landusebuilt"]
landuse_contrast <- function(name, coef_name) {
  if (is.na(coef_name)) {
    est <- -b_built; se <- sqrt(v_built)
  } else {
    est <- cf[coef_name] - b_built
    se  <- sqrt(V[coef_name, coef_name] + v_built - 2 * V[coef_name, "Landusebuilt"])
  }
  tibble::tibble(variable = name, estimate = unname(est), se = unname(se), type = "Land use (vs. built)")
}
landuse_df <- dplyr::bind_rows(
  landuse_contrast("Cropland",  NA),
  landuse_contrast("Grassland", "Landusegrassland"),
  landuse_contrast("Shrubs",    "Landuseshrubs"),
  landuse_contrast("Trees",     "Landusetrees")
)

plot_df <- dplyr::bind_rows(env_df, landuse_df) %>%
  mutate(
    ci_lo = estimate - 1.96 * se,
    ci_hi = estimate + 1.96 * se,
    p_value = 2 * pnorm(-abs(estimate / se)),
    sig = ifelse(ci_lo > 0 | ci_hi < 0, "Significant", "Non-significant"),
    variable = factor(variable, levels = c("Trees", "Shrubs", "Grassland", "Cropland",
                                            "HII", "Elevation", "Bio_15", "Bio_4"))
  )

readr::write_csv(plot_df, "output/phylo_export/main_model_forest_data_built_ref.csv")
print(plot_df %>% dplyr::select(variable, type, estimate, ci_lo, ci_hi, p_value, sig))

PUB_FONT <- "Helvetica"

p <- ggplot(plot_df, aes(x = estimate, y = variable, color = type, shape = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.22, linewidth = 0.7) +
  geom_point(size = 3) +
  scale_color_manual(values = c("Environmental" = "#1b7837", "Land use (vs. built)" = "#762a83")) +
  scale_shape_manual(values = c("Significant" = 16, "Non-significant" = 1)) +
  labs(
    title = "Drivers of global migratory butterfly richness",
    subtitle = "Linear GAM (spatial smooth retained), built-up land as reference category",
    x = "Coefficient estimate (95% CI)", y = NULL, color = NULL, shape = NULL
  ) +
  theme_minimal(base_size = 11, base_family = PUB_FONT) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(face = "bold", color = "black"),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, colour = "grey30")
  )

dir.create("output/BAM/figures", showWarnings = FALSE, recursive = TRUE)
ggsave("output/BAM/figures/Figure2d_forest_main_model_built_ref.png", p, width = 7, height = 5,
       dpi = 600, bg = "white", device = ragg::agg_png)
ggsave("output/phylo_export/Figure2d_forest_main_model_built_ref.pdf", p, width = 7, height = 5, device = "pdf")

cat("\nSaved: output/BAM/figures/Figure2d_forest_main_model_built_ref.png\n")
cat("Saved: output/phylo_export/Figure2d_forest_main_model_built_ref.pdf\n")
