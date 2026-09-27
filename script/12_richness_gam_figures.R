# Richness GAM diagnostic figures.

library(mgcv)
library(ggplot2)
library(gratia)
library(patchwork)
library(dplyr)

cat("\n--- [1/4] 环境就绪，正在加载模型与数据... ---\n")

model_final <- readRDS("/Users/hlii0385/Desktop/Phd_Haiyu_LI/Rapoport-s-rule-and-Bergmann-s-rule-of-Migratory-butterflies/output/BAM/output/butterfly_bam_model.rds")

file_path <- "/Users/hlii0385/Desktop/Phd_Haiyu_LI/Rapoport-s-rule-and-Bergmann-s-rule-of-Migratory-butterflies/output/BAM/output/cleaned_data_for_SAR.csv"
df_aggressive <- read.csv(file_path)
df_aggressive$Landuse <- factor(df_aggressive$Landuse)
df_aggressive$Landuse <- relevel(df_aggressive$Landuse, ref = "cropland")

map_data <- read.csv("/Users/hlii0385/Desktop/Phd_Haiyu_LI/Rapoport-s-rule-and-Bergmann-s-rule-of-Migratory-butterflies/output/BAM/output/spatial_residuals_map_data.csv")

boot_data <- read.csv("/Users/hlii0385/Desktop/Phd_Haiyu_LI/Rapoport-s-rule-and-Bergmann-s-rule-of-Migratory-butterflies/output/BAM/output/moran_bootstrap_distribution.csv")

cat("✅ 模型与数据加载成功！\n")

cat("\n--- [2/4] 正在生成 Main Figure (聚光灯截距平移图)... ---\n")

get_smooth_curve <- function(mod, term_name, var_name) {
  dat <- smooth_estimates(mod, select = term_name)
  x_val <- if (".value" %in% names(dat)) dat$.value else dat[[var_name]]
  y_val <- if (".estimate" %in% names(dat)) dat$.estimate else dat$partial
  x_scaled <- (x_val - min(x_val)) / (max(x_val) - min(x_val))
  data.frame(Variable = var_name, X_Relative = x_scaled, Effect = y_val)
}

curve_bio4  <- get_smooth_curve(model_final, "s(Bio_4)", "Bio_4")
curve_bio15 <- get_smooth_curve(model_final, "s(Bio_15)", "Bio_15")
curve_elev  <- get_smooth_curve(model_final, "s(Elevation)", "Elevation")

hii_seq <- seq(min(df_aggressive$HII, na.rm = TRUE),
               max(df_aggressive$HII, na.rm = TRUE), length.out = 100)

pred_df_hii <- df_aggressive[1:100, ]
pred_df_hii$HII <- hii_seq

hii_terms <- predict(model_final, newdata = pred_df_hii, type = "terms")

curve_hii <- data.frame(
  Variable = "HII",
  X_Relative = seq(0, 1, length.out = 100),
  Effect = hii_terms[, "HII"]
)

all_curves <- bind_rows(curve_bio4, curve_bio15, curve_elev, curve_hii)

lu_levels <- levels(df_aggressive$Landuse)
pred_df_lu <- df_aggressive[1:length(lu_levels), ]
pred_df_lu$Landuse <- factor(lu_levels, levels = lu_levels)

lu_terms <- predict(model_final, newdata = pred_df_lu, type = "terms")

landuse_eff <- data.frame(
  Landuse = lu_levels,
  Shift   = lu_terms[, "Landuse"]
)

ref_shift <- landuse_eff$Shift[landuse_eff$Landuse == "cropland"]
landuse_eff$Shift <- landuse_eff$Shift - ref_shift

plot_data_shifted <- cross_join(all_curves, landuse_eff) |>
  mutate(
    Shifted_Effect = Effect + Shift,

    Target_Group = "built",

    Is_Target = (Landuse == Target_Group)
  )

main_figure <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +

  geom_line(
    data = plot_data_shifted |> filter(!Is_Target),
    aes(x = X_Relative, y = Shifted_Effect, group = interaction(Variable, Landuse)),
    color = "grey60", alpha = 0.25, linewidth = 0.8
  ) +

  geom_line(
    data = plot_data_shifted |> filter(Is_Target),
    aes(x = X_Relative, y = Shifted_Effect, color = Variable),
    linewidth = 1.5
  ) +

  scale_color_manual(values = c("Bio_4" = "#d73027",
                                "Bio_15" = "#4575b4",
                                "Elevation" = "#1b7837",
                                "HII" = "#762a83")) +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Relative Drivers of Global Butterfly Richness",
    subtitle = paste0("Highlighted group: [", unique(plot_data_shifted$Target_Group),
                      "] | Background: Other land use categories"),
    x = "Relative Gradient (0 to 1)",
    y = "Partial Effect on Richness (Intercept Shifted)"
  )
ggsave("/Users/hlii0385/Desktop/Phd_Haiyu_LI/Rapoport-s-rule-and-Bergmann-s-rule-of-Migratory-butterflies/Figure_1_Main_Spotlight.pdf", main_figure, width = 10, height = 7, dpi = 300)
ggsave("/Users/hlii0385/Desktop/Phd_Haiyu_LI/Rapoport-s-rule-and-Bergmann-s-rule-of-Migratory-butterflies/Figure_1_Main_Spotlight.png", main_figure, width = 10, height = 7, dpi = 300, bg = "white")

cat("\n--- [3/4] 正在生成 Support Figures (模型诊断与稳健性)... ---\n")

sf1_map <- ggplot(map_data, aes(x = lon, y = lat, color = res)) +
  geom_point(alpha = 0.6, size = 0.5) +
  scale_color_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", midpoint = 0,
                        name = "Pearson\nResidual") +
  theme_minimal(base_size = 14) +
  labs(title = "Spatial Distribution of Model Residuals",
       subtitle = "Residuals exhibit a near-random spatial pattern after incorporating the s(lon, lat) smooth term",
       x = "Longitude", y = "Latitude") +
  theme(legend.position = "right",
        panel.background = element_rect(fill = "aliceblue", color = NA))

ggsave("/Users/hlii0385/Desktop/Phd_Haiyu_LI/Rapoport-s-rule-and-Bergmann-s-rule-of-Migratory-butterflies/Support_Figure_1_Residual_Map.png", sf1_map, width = 12, height = 6, dpi = 300, bg = "white")

mean_I <- mean(boot_data$Moran_I)
ci_lower <- quantile(boot_data$Moran_I, 0.025)
ci_upper <- quantile(boot_data$Moran_I, 0.975)

sf2_moran <- ggplot(boot_data, aes(x = Moran_I)) +
  geom_histogram(aes(y = after_stat(density)), bins = 15, fill = "#a6bddb", color = "black", alpha = 0.7) +
  geom_density(color = "#045a8d", linewidth = 1) +
  geom_vline(xintercept = mean_I, color = "#cb181d", linetype = "dashed", linewidth = 1.2) +
  geom_vline(xintercept = c(ci_lower, ci_upper), color = "black", linetype = "dotted", linewidth = 1) +
  theme_minimal(base_size = 14) +
  labs(title = "Robustness Check: Bootstrapped Residual Moran's I",
       subtitle = paste0("100 Resamples (N=10000 each) | Mean: ", round(mean_I, 4),
                         " | 95% CI: [", round(ci_lower, 4), ", ", round(ci_upper, 4), "]"),
       x = "Global Moran's I of Model Residuals", y = "Density") +
  annotate("text", x = mean_I, y = max(density(boot_data$Moran_I)$y) * 1.05,
           label = "Mean", vjust = -0.5, color = "#cb181d") +
  coord_cartesian(clip = "off")

ggsave("/Users/hlii0385/Desktop/Phd_Haiyu_LI/Rapoport-s-rule-and-Bergmann-s-rule-of-Migratory-butterflies/Support_Figure_2_Moran_Distribution.png", sf2_moran, width = 8, height = 6, dpi = 300, bg = "white")

sf3_spatial_smooth <- draw(model_final, select = "s(lon,lat)") +
  theme_minimal(base_size = 14) +
  scale_fill_viridis_c(option = "mako", name = "Spatial\nEffect") +
  labs(title = "Extracted Spatial Smooth Effect",
       subtitle = "Based on s(lon, lat, bs='ds', k=1500)",
       x = "Longitude", y = "Latitude")

ggsave("/Users/hlii0385/Desktop/Phd_Haiyu_LI/Rapoport-s-rule-and-Bergmann-s-rule-of-Migratory-butterflies/Support_Figure_3_Spatial_Smooth.png", sf3_spatial_smooth, width = 10, height = 6, dpi = 300, bg = "white")

cat("\n--- [4/4] 恭喜！所有高清图表 (PDF/PNG) 均已成功生成并保存在指定目录中！ ---\n")
