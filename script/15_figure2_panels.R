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

BASE_SIZE <- 12  # text size once the figure is placed at 6.5 in width (manuscript body text)
theme_fig2 <- theme_classic(base_size = BASE_SIZE, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.3),
    axis.text = element_text(colour = "black", size = BASE_SIZE),
    axis.title = element_text(size = BASE_SIZE),
    legend.background = element_blank(),
    legend.key.width = unit(4, "mm"), legend.key.height = unit(3, "mm"),
    legend.key.spacing = unit(1, "mm"), legend.spacing = unit(0, "mm"),
    axis.title.y = element_text(margin = margin(r = 3)), axis.title.x = element_text(margin = margin(t = 3)),
    legend.text = element_text(size = BASE_SIZE),
    legend.title = element_text(size = BASE_SIZE),
    legend.position = "top", legend.justification = "left", legend.direction = "horizontal",
    legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(1, "mm"),
    plot.margin = margin(3, 4, 3, 3)
  )

# Horizontal colour bar with the minimum and maximum written at its two ends (as in Fig. 1)
end_bar <- function(lo, hi, title) {
  ggplot(data.frame(x = seq(0, 1, length.out = 256))) +
    geom_raster(aes(x = x, y = 0, fill = x)) +
    scale_fill_viridis_c(guide = "none") +
    annotate("rect", xmin = 0, xmax = 1, ymin = -0.5, ymax = 0.5, fill = NA, colour = "black", linewidth = 0.2) +
    annotate("text", x = -0.03, y = 0, label = lo, hjust = 1, size = BASE_SIZE / .pt, family = "Arial") +
    annotate("text", x = 1.03, y = 0, label = hi, hjust = 0, size = BASE_SIZE / .pt, family = "Arial") +
    annotate("text", x = 0.5, y = 1.1, label = title, vjust = 0, size = BASE_SIZE / .pt, family = "Arial") +
    coord_cartesian(xlim = c(0, 1), ylim = c(-0.5, 0.5), expand = FALSE, clip = "off") +
    theme_void()
}
lab_range <- expression("Range size (" * log[10] * ", " * km^2 * ")")
lab_ws <- function(w) bquote(.(w) * " (" * log[10] * ", cm)")

# Panel a
m_pattern <- lmer(log10(range_km2) ~ prop_mean + prop_within + season + (1 | species),
                   data = final_df_lat)

pred_pattern <- predict_response(m_pattern, terms = "prop_mean", back_transform = FALSE) %>%
  as_tibble() %>% rename(prop_mean = x, log10_range = predicted)

p2a <- ggplot() +
  geom_ribbon(data = pred_pattern, aes(x = prop_mean, ymin = conf.low, ymax = conf.high),
              fill = "grey40", alpha = 0.25) +
  geom_point(data = final_df, aes(x = prop_mean, y = log10(range_km2), color = mean_bio4),
             alpha = 0.5, size = 0.8, stroke = 0) +
  scale_color_viridis_c(option = "viridis", na.value = "grey75") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(limits = c(1.4, NA)) +  # empty band at the bottom for the colour bar
  geom_line(data = pred_pattern, aes(x = prop_mean, y = log10_range),
            color = "black", linewidth = 0.6) +
  labs(x = "Mean proportion of range in the tropics", y = lab_range,
       color = "Temperature\nseasonality") +
  guides(color = "none") +
  theme_fig2 +
  theme(plot.margin = margin(3, 4, 3, 3))
bio4_rng <- range(final_df$mean_bio4, na.rm = TRUE)
p2a <- p2a + inset_element(end_bar(round(bio4_rng[1]), round(bio4_rng[2]), "Temperature seasonality"),
                           left = 0.30, right = 0.72, bottom = 0.03, top = 0.07, align_to = "panel",
                           clip = FALSE, ignore_tag = TRUE)

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
             alpha = 0.55, size = 0.8, stroke = 0) +
  geom_line(data = df_preds2b, aes(x = abs_lat, y = log10_WS, color = WS_type), linewidth = 0.6) +
  scale_color_manual(values = c("Lower" = viridisLite::viridis(5)[4], "Upper" = viridisLite::viridis(5)[1])) +
  scale_fill_manual(values = c("Lower" = viridisLite::viridis(5)[4], "Upper" = viridisLite::viridis(5)[1])) +
  labs(x = "Absolute latitude (°)", y = lab_ws("Wingspan"),
       color = "Wingspan", fill = "Wingspan") +
  theme_fig2 +
  theme(legend.position = "inside", legend.position.inside = c(0.99, 0.99), legend.justification.inside = c(1, 1),
        legend.direction = "vertical")

# Panel c
m_interact <- lmer(log10(range_km2) ~ abs_lat * log10(WS_L) + prop_within + season + (1 | species),
                    data = final_df_lat)

lat_vals <- c(0, 20, 40, 60)
pred_c <- predict_response(m_interact, terms = c("WS_L [n=50]", "abs_lat [0,20,40,60]"),
                            back_transform = FALSE) %>%
  as_tibble() %>%
  rename(WS_L = x, fit = predicted) %>%
  mutate(WS_L = log10(WS_L)) %>%
  mutate(lat_group = factor(paste0(group, "°"), levels = paste0(lat_vals, "°")))

p2c <- ggplot(pred_c, aes(x = WS_L, y = fit, color = lat_group, fill = lat_group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), color = NA, alpha = 0.15) +
  geom_line(linewidth = 0.6) +
  scale_color_viridis_d(option = "viridis", end = 0.85) +
  scale_fill_viridis_d(option = "viridis", end = 0.85) +
  labs(x = lab_ws("Lower wingspan"), y = lab_range,
       color = "Absolute latitude", fill = "Absolute latitude") +
  guides(color = guide_legend(nrow = 2), fill = guide_legend(nrow = 2)) +
  theme_fig2 +
  theme(legend.position = "inside", legend.position.inside = c(0.01, 0.99), legend.justification.inside = c(0, 1),
        legend.direction = "horizontal", legend.title.position = "top")

# Panel d: main GAM (black only; filled = 95% CI excludes zero)
forest_df <- read_csv("output/phylo_export/main_model_forest_data_built_ref.csv", show_col_types = FALSE) %>%
  mutate(variable = recode(variable, "Bio_4" = "Temperature seasonality", "Bio_15" = "Precipitation seasonality",
                           "HII" = "Human influence", "Trees" = "Forest", "Shrubs" = "Shrubland"),
         group = factor(ifelse(type == "Environmental", "Environment", "Land cover\n(vs built-up)"),
                        levels = c("Environment", "Land cover\n(vs built-up)")),
         variable = factor(variable, levels = rev(c("Temperature seasonality", "Precipitation seasonality",
                                                    "Elevation", "Human influence",
                                                    "Cropland", "Grassland", "Shrubland", "Forest"))))

p2d <- ggplot(forest_df, aes(x = estimate, y = variable, shape = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.25, linewidth = 0.5, colour = "black") +
  geom_point(size = 2.2, colour = "black", fill = "white", stroke = 0.6) +
  scale_shape_manual(values = c("Significant" = 16, "Non-significant" = 21), name = NULL,
                     breaks = c("Significant", "Non-significant"), labels = c("P < 0.05", "P ≥ 0.05")) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Coefficient estimate", y = NULL) +
  theme_fig2 +
  theme(strip.background = element_blank(), strip.text = element_blank(), panel.spacing.y = unit(3, "mm"),
        legend.position = "inside", legend.position.inside = c(0.02, 0.02), legend.justification.inside = c(0, 0),
        legend.direction = "vertical")

# Compose Figure 2 (6.5 in wide, the text width of the manuscript)
fig2 <- (p2a + p2b) / (p2c + p2d + plot_layout(widths = c(0.9, 1.1))) + plot_layout(heights = c(1, 1.05)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = BASE_SIZE, family = "Arial"))
ggsave(file.path(out_dir, "Figure2_final.png"), fig2, width = 6.5, height = 5.3, units = "in",
       dpi = 600, device = ragg::agg_png, bg = "white")
ggsave(file.path(out_dir, "Figure2_final.pdf"), fig2, width = 6.5, height = 5.3, units = "in",
       device = cairo_pdf)
cat("Saved Figure 2 to", out_dir, "\n")

PANEL_W <- 5.5; PANEL_H <- 4.25; PANEL_DPI <- 400

# Fig. S3: phylogenetic version of c (styled as in the Supporting Information)
theme_s3 <- theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    legend.background = element_rect(fill = alpha("white", 0.75), color = NA),
    legend.key.size = unit(0.3, "cm"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.position = "inside", legend.position.inside = c(0.82, 0.86)
  )
pp <- read_csv("output/phylo_export/interaction_phylo_predictions.csv", show_col_types = FALSE) %>%
  mutate(lat_group = factor(paste0(abs_lat, "°"), levels = paste0(c(0, 20, 40, 60), "°")))
p_s3 <- ggplot(pp, aes(x = WS_U, y = fit, colour = lat_group, fill = lat_group)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), colour = NA, alpha = 0.15) +
  geom_line(linewidth = 1.1) +
  scale_x_log10() +
  scale_colour_viridis_d(option = "viridis") + scale_fill_viridis_d(option = "viridis") +
  labs(x = "Upper wingspan (log10 cm)", y = expression(log[10] * " Range size (" * km^2 * ")"),
       colour = "Absolute\nlatitude", fill = "Absolute\nlatitude") +
  theme_s3
dir.create("output/SI", showWarnings = FALSE)
ggsave("output/SI/FigS_interaction_phylo.png", p_s3, width = PANEL_W, height = PANEL_H, dpi = PANEL_DPI, bg = "white")
