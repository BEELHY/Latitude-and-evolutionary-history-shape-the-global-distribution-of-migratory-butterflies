# Figure 2 panels.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(ggplot2)
  library(lme4)
  library(lmerTest)
  library(ggeffects)
  library(patchwork)
  library(viridis)
  library(mgcv)
})

out_dir <- "output/Manuscript/reproducibility_code"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Load data
trait_range  <- read_csv("output/trait_range.csv", show_col_types = FALSE)
tree_sp      <- ape::read.tree("updatedata/phylogeny_matched.tre")$tip.label
# BIO4 colour: phylogeny subset
mean_bio4_ss <- read_csv("updatedata/species_season_metrics.csv", show_col_types = FALSE) %>%
  dplyr::filter(species %in% tree_sp) %>%
  dplyr::select(species, season, mean_bio4)
df_lat_sp    <- read_csv("output/df_lat_sp.csv", show_col_types = FALSE)

final_df <- trait_range %>%
  mutate(species = str_replace_all(species, " ", "_")) %>%
  group_by(species) %>%
  mutate(prop_mean = mean(prop_tropics, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(prop_within = prop_tropics - prop_mean) %>%
  left_join(mean_bio4_ss %>% mutate(species = str_replace_all(species, " ", "_")),
            by = c("species", "season"))

final_df_lat <- final_df %>%
  left_join(df_lat_sp %>% dplyr::select(species, abs_lat), by = "species")

BASE_SIZE <- 12
LEGEND_POS <- c(0.82, 0.86)  # Top-right legend
theme_fig2 <- theme_classic(base_size = BASE_SIZE) +
  theme(
    plot.title = element_text(face = "bold", size = BASE_SIZE),
    legend.background = element_rect(fill = alpha("white", 0.75), color = NA),
    legend.key.size = unit(0.3, "cm"),
    legend.text = element_text(size = BASE_SIZE - 2),
    legend.title = element_text(size = BASE_SIZE - 2),
    legend.position = "inside", legend.position.inside = LEGEND_POS
  )

# Panel a
m_pattern <- lmer(log10(range_km2) ~ prop_mean + prop_within + season + (1 | species),
                   data = final_df_lat)

pred_pattern <- predict_response(m_pattern, terms = "prop_mean", back_transform = FALSE) %>%
  as_tibble() %>% rename(prop_mean = x, log10_range = predicted)

p2a <- ggplot() +
  geom_ribbon(data = pred_pattern, aes(x = prop_mean, ymin = conf.low, ymax = conf.high),
              fill = "#D55E00", alpha = 0.2) +
  geom_point(data = final_df, aes(x = prop_mean, y = log10(range_km2), color = mean_bio4),
             alpha = 0.3, size = 1) +
  scale_color_viridis_c(option = "viridis") +
  # Blank margin for legend
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.20))) +
  geom_line(data = pred_pattern, aes(x = prop_mean, y = log10_range),
            color = "#D55E00", linewidth = 1.2) +
  labs(x = "Mean proportion in tropical area",
       y = expression(log[10] * " Range size (" * km^2 * ")"),
       color = "Temp.\nseasonality\n(BIO4)") +
  theme_fig2 +
  # Anchor legend top-left
  theme(legend.position.inside = c(0.85, 0.97),
        legend.justification.inside = c(0, 1))

# Panel b
m_lat_L <- lmer(log10(WS_L) ~ abs_lat + (1 | Family),
                 data = df_lat_sp %>% dplyr::filter(is.finite(WS_L), WS_L > 0))
m_lat_U <- lmer(log10(WS_U) ~ abs_lat + (1 | Family),
                 data = df_lat_sp %>% dplyr::filter(is.finite(WS_U), WS_U > 0))

pred_L <- predict_response(m_lat_L, terms = "abs_lat", back_transform = FALSE) %>%
  as_tibble() %>% mutate(WS_type = "Lower")
pred_U <- predict_response(m_lat_U, terms = "abs_lat", back_transform = FALSE) %>%
  as_tibble() %>% mutate(WS_type = "Upper")
df_preds2b <- bind_rows(pred_L, pred_U) %>% rename(abs_lat = x, log10_WS = predicted)

df_plot2b <- df_lat_sp %>%
  dplyr::filter(is.finite(WS_L), WS_L > 0, is.finite(WS_U), WS_U > 0) %>%
  pivot_longer(cols = c(WS_L, WS_U), names_to = "WS_type", values_to = "WS_value") %>%
  mutate(WS_type = recode(WS_type, WS_U = "Upper", WS_L = "Lower"))

p2b <- ggplot() +
  geom_ribbon(data = df_preds2b, aes(x = abs_lat, ymin = conf.low, ymax = conf.high, fill = WS_type),
              alpha = 0.2) +
  geom_point(data = df_plot2b, aes(x = abs_lat, y = log10(WS_value), color = WS_type),
             alpha = 0.4, size = 1.1, shape = 16) +
  geom_line(data = df_preds2b, aes(x = abs_lat, y = log10_WS, color = WS_type), linewidth = 1.2) +
  scale_color_manual(values = c("Lower" = "#56B4E9", "Upper" = "#D55E00")) +
  scale_fill_manual(values = c("Lower" = "#56B4E9", "Upper" = "#D55E00")) +
  labs(x = "Absolute latitude (°)", y = expression(log[10] * " Wingspan (cm)"),
       color = "Wingspan", fill = "Wingspan") +
  theme_fig2

# Panel c
m_interact <- lmer(log10(range_km2) ~ abs_lat * log10(WS_L) + prop_within + season + (1 | species),
                    data = final_df_lat)

lat_vals <- c(0, 20, 40, 60)
pred_c <- predict_response(m_interact, terms = c("WS_L [n=50]", "abs_lat [0,20,40,60]"),
                            back_transform = FALSE) %>%
  as_tibble() %>%
  rename(WS_L = x, fit = predicted) %>%
  mutate(lat_group = factor(paste0(group, "°"), levels = paste0(lat_vals, "°")))

p2c <- ggplot(pred_c, aes(x = WS_L, y = fit, color = lat_group, fill = lat_group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), color = NA, alpha = 0.15) +
  geom_line(linewidth = 1.1) +
  scale_x_log10() +
  scale_color_viridis_d(option = "viridis") +
  scale_fill_viridis_d(option = "viridis") +
  labs(x = "Wingspan (log10 cm)",
       y = expression(log[10] * " Range size (" * km^2 * ")"),
       color = "Absolute\nlatitude", fill = "Absolute\nlatitude") +
  theme_fig2

# Panel d: main GAM
forest_df <- read_csv("output/phylo_export/main_model_forest_data_built_ref.csv", show_col_types = FALSE) %>%
  mutate(variable = recode(variable, "Bio_4" = "Temperature\nseasonality", "Bio_15" = "Precipitation\nseasonality"),
         # Row order, top to bottom
         variable = factor(variable, levels = c("Trees", "Shrubs", "Grassland", "Cropland",
                                                 "HII", "Precipitation\nseasonality",
                                                 "Temperature\nseasonality", "Elevation")))

p2d <- ggplot(forest_df, aes(x = estimate, y = variable, color = type, shape = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.22, linewidth = 0.7) +
  geom_point(size = 2.6) +
  scale_color_manual(values = c("Environmental" = "#1b7837", "Land use (vs. built)" = "#762a83")) +
  scale_shape_manual(values = c("Significant" = 16, "Non-significant" = 1), guide = "none") +
  labs(x = "Coefficient estimate (95% CI)", y = NULL, color = NULL) +
  theme_fig2

# Export panels separately
panel_dir <- file.path(out_dir, "panels")
dir.create(panel_dir, showWarnings = FALSE, recursive = TRUE)

add_tag <- function(p, tag) {
  p + labs(tag = tag) +
    theme(plot.tag = element_text(face = "bold", size = BASE_SIZE + 2),
          plot.tag.position = c(0.02, 0.98))
}

PANEL_W <- 5.5; PANEL_H <- 4.25; PANEL_DPI <- 400

save_panel <- function(p, name) {
  ggsave(file.path(panel_dir, paste0(name, ".png")), p, width = PANEL_W, height = PANEL_H,
         dpi = PANEL_DPI, bg = "white")
  ggsave(file.path(panel_dir, paste0(name, ".pdf")), p, width = PANEL_W, height = PANEL_H)
}

save_panel(add_tag(p2a, "a"), "panel_a")
save_panel(add_tag(p2b, "b"), "panel_b")
save_panel(add_tag(p2c, "c"), "panel_c")
save_panel(add_tag(p2d, "d"), "panel_d")

cat(sprintf("Saved 4 panels (PNG @ %d dpi + PDF) to: %s\n", PANEL_DPI, panel_dir))

# Fig. S3: phylogenetic version of c
pp <- read_csv("output/phylo_export/interaction_phylo_predictions.csv", show_col_types = FALSE) %>%
  mutate(lat_group = factor(paste0(abs_lat, "°"), levels = paste0(c(0, 20, 40, 60), "°")))
p_s3 <- ggplot(pp, aes(x = WS_U, y = fit, colour = lat_group, fill = lat_group)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), colour = NA, alpha = 0.15) +
  geom_line(linewidth = 1.1) +
  scale_x_log10() +
  scale_colour_viridis_d(option = "viridis") + scale_fill_viridis_d(option = "viridis") +
  labs(x = "Upper wingspan (log10 cm)", y = expression(log[10] * " Range size (" * km^2 * ")"),
       colour = "Absolute\nlatitude", fill = "Absolute\nlatitude") +
  theme_fig2
dir.create("output/SI", showWarnings = FALSE)
ggsave("output/SI/FigS_interaction_phylo.png", p_s3, width = PANEL_W, height = PANEL_H, dpi = PANEL_DPI, bg = "white")
cat("Next: run script/16_figure2_compose.py to lay them out into the final Figure 2 grid.\n")
