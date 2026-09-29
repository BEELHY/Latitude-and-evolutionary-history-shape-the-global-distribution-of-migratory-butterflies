# Forest plot data.

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

# Environmental predictors
env_vars <- c(Bio_4 = "Bio_4", Bio_15 = "Bio_15", Elevation = "Elevation", HII = "HII")
env_df <- tibble::tibble(
  variable = names(env_vars),
  estimate = cf[unname(env_vars)],
  se       = sqrt(diag(V)[unname(env_vars)]),
  type     = "Environmental"
)

# Land use vs built
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

# Table S3 (GAM columns) and Fig. S5: the three equal-area models, land cover vs built-up
contrasts <- function(mod, label) {
  cf <- coef(mod); V <- vcov(mod)
  terms <- intersect(c("Bio_4", "Bio_15", "Elevation", "HII", "NRI", "NTI"), names(cf))
  env <- tibble::tibble(term = terms, est = cf[terms], se = sqrt(diag(V)[terms]))
  bb <- "Landusebuilt"
  lu <- dplyr::bind_rows(
    tibble::tibble(term = "Cropland", est = -cf[bb], se = sqrt(V[bb, bb])),
    lapply(c(Grassland = "Landusegrassland", Shrubs = "Landuseshrubs", Trees = "Landusetrees"), function(k)
      tibble::tibble(est = cf[k] - cf[bb], se = sqrt(V[k, k] + V[bb, bb] - 2 * V[k, bb]))) %>%
      dplyr::bind_rows(.id = "term"))
  dplyr::bind_rows(env, lu) %>%
    mutate(model = label, lo = est - 1.96 * se, hi = est + 1.96 * se, p = 2 * pnorm(-abs(est / se)))
}
gam_rows <- dplyr::bind_rows(contrasts(bundle$m_A, "main"), contrasts(bundle$m_B, "sp247"),
                             contrasts(bundle$m_C, "sp247_NRI_NTI"))
si_file <- "output/SI/richness_models_builtref.csv"
old <- readr::read_csv(si_file, show_col_types = FALSE) %>%
  dplyr::filter(!model %in% c("main", "sp247", "sp247_NRI_NTI"))
readr::write_csv(dplyr::bind_rows(gam_rows %>% dplyr::select(model, term, est, se, lo, hi, p), old), si_file)

lab <- c(Bio_4 = "Temperature\nseasonality", Bio_15 = "Precipitation\nseasonality", Elevation = "Elevation",
         HII = "Human influence", Cropland = "Cropland", Grassland = "Grassland", Shrubs = "Shrubland",
         Trees = "Forest", NRI = "Net relatedness\nindex", NTI = "Nearest taxon\nindex")
panel <- c(main = "a  All species (main model)", sp247 = "b  247 phylogenetically matched species",
           sp247_NRI_NTI = "c  247 species + community phylogenetic structure")
s5 <- gam_rows %>%
  mutate(group = dplyr::case_when(term %in% c("NRI", "NTI") ~ "Phylogenetic structure",
                                  term %in% c("Cropland", "Grassland", "Shrubs", "Trees") ~ "Land cover (vs built-up)",
                                  TRUE ~ "Environment"),
         sig = ifelse(lo > 0 | hi < 0, "P < 0.05", "n.s."),
         term = factor(lab[term], levels = rev(lab)), model = factor(panel[model], levels = panel))
p_s5 <- ggplot(s5, aes(x = est, y = term, colour = group, shape = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.25, linewidth = 0.6, orientation = "y") +
  geom_point(size = 2.2, fill = "white", stroke = 0.8) +
  facet_wrap(~ model, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = c("Environment" = "#1b7837", "Land cover (vs built-up)" = "#762a83",
                                 "Phylogenetic structure" = "#E08214"), name = NULL) +
  scale_shape_manual(values = c("P < 0.05" = 16, "n.s." = 21), name = NULL) +
  labs(x = "Coefficient estimate (95% confidence interval, log scale)", y = NULL) +
  theme_classic(base_size = 11, base_family = "Arial") +
  theme(legend.position = "top", strip.background = element_blank(),
        strip.text = element_text(face = "bold", hjust = 0), axis.text = element_text(colour = "black"))
ggsave("output/SI/FigS_NRI_NTI.png", p_s5, width = 11, height = 5.2, dpi = 400, bg = "white", device = ragg::agg_png)
cat("Saved: output/SI/richness_models_builtref.csv (GAM rows) and output/SI/FigS_NRI_NTI.png\n")
