# ==============================================================================
# 全球蝴蝶多样性 BAM 模型 - 最终版可视化脚本
# 包含：Main Figure (生态学驱动力) & Support Figures (统计学稳健性检验)
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. 加载必备扩展包
# ------------------------------------------------------------------------------
install.packages(c("mgcv", "ggplot2", "gratia", "patchwork"))
library(mgcv)       # BAM 模型基础包
library(ggplot2)    # 绘图底层逻辑
library(gratia)     # GAM/BAM 模型专属高颜值可视化提取包
library(patchwork)  # 优雅的图表拼图神器

cat("\n--- [1/4] 环境就绪，正在加载模型与数据... ---\n")

# ------------------------------------------------------------------------------
# 2. 读取持久化的模型与诊断数据
# ------------------------------------------------------------------------------
# 读取核心 BAM 模型 (加载时间视内存大小而定，请耐心等待)
model_final <- readRDS("butterfly_bam_model.rds")

# 读取空间残差诊断数据
map_data <- read.csv("spatial_residuals_map_data.csv")

# 读取 Moran's I 的 100 次抽样诊断数据
boot_data <- read.csv("moran_bootstrap_distribution.csv")

cat("✅ 模型与数据加载成功！\n")

# ==============================================================================
# 3. 绘制 Main Figure：驱动全球蝴蝶多样性的核心机制
# ==============================================================================
cat("\n--- [2/4] 正在生成 Main Figure (生态学驱动力核心图)... ---\n")

# A. 提取自然环境因子 (非线性平滑项 Smooth effects)
plot_bio4 <- draw(model_final, select = "s(Bio_4)") + 
  theme_minimal(base_size = 14) + 
  labs(title = "Temperature Seasonality", x = "Bio 4", y = "Partial Effect on Richness")

plot_bio15 <- draw(model_final, select = "s(Bio_15)") + 
  theme_minimal(base_size = 14) + 
  labs(title = "Precipitation Seasonality", x = "Bio 15", y = "")

plot_elev <- draw(model_final, select = "s(Elevation)") + 
  theme_minimal(base_size = 14) + 
  labs(title = "Elevation", x = "Elevation (m)", y = "")

# B. 提取人为影响因子 (线性参数项 Parametric effects)
# HII (连续变量线性项)
plot_hii <- parametric_effects(model_final, term = "HII") |> 
  draw() + 
  theme_minimal(base_size = 14) + 
  labs(title = "Human Influence Index", x = "HII", y = "Partial Effect on Richness")

# Landuse (分类变量参数项)
plot_landuse <- parametric_effects(model_final, term = "Landuse") |> 
  draw() + 
  theme_minimal(base_size = 14) + 
  labs(title = "Land Use Type", x = "Landuse", y = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

# C. 完美拼图 (上排自然因子，下排人为因子)
main_figure <- (plot_bio4 | plot_bio15 | plot_elev) / 
               (plot_hii | plot_landuse) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 18, face = "bold"))

# 保存高分辨率 Main Figure
ggsave("Figure_1_Main_Drivers.pdf", main_figure, width = 16, height = 10, dpi = 300)
ggsave("Figure_1_Main_Drivers.png", main_figure, width = 16, height = 10, dpi = 300, bg = "white")

# ==============================================================================
# 4. 绘制 Support Figures：严谨的统计学稳健性证明
# ==============================================================================
cat("\n--- [3/4] 正在生成 Support Figures (模型诊断与稳健性)... ---\n")

# ----------------------------------------------------------------
# Support Figure 1: 空间残差分布地图 (证明偏误已被吸收)
# ----------------------------------------------------------------
sf1_map <- ggplot(map_data, aes(x = lon, y = lat, color = res)) +
  geom_point(alpha = 0.6, size = 0.5) +
  scale_color_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", midpoint = 0,
                        name = "Pearson\nResidual") +
  theme_minimal(base_size = 14) +
  labs(title = "Spatial Distribution of Model Residuals",
       subtitle = "Residuals exhibit a near-random spatial pattern after incorporating the s(lon, lat) smooth term",
       x = "Longitude", y = "Latitude") +
  theme(legend.position = "right", 
        panel.background = element_rect(fill = "aliceblue", color = NA)) # 加一点淡淡的海洋蓝背景

ggsave("Support_Figure_1_Residual_Map.png", sf1_map, width = 12, height = 6, dpi = 300, bg = "white")

# ----------------------------------------------------------------
# Support Figure 2: Moran's I 重抽样分布 (证明残差空间自相关显著降低)
# ----------------------------------------------------------------
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

ggsave("Support_Figure_2_Moran_Distribution.png", sf2_moran, width = 8, height = 6, dpi = 300, bg = "white")

# ----------------------------------------------------------------
# Support Figure 3: s(lon, lat) 提取的空间平滑效应基底
# ----------------------------------------------------------------
sf3_spatial_smooth <- draw(model_final, select = "s(lon,lat)") + 
  theme_minimal(base_size = 14) +
  scale_fill_viridis_c(option = "mako", name = "Spatial\nEffect") +
  labs(title = "Extracted Spatial Smooth Effect",
       subtitle = "Based on s(lon, lat, bs='ds', k=1500)",
       x = "Longitude", y = "Latitude")

ggsave("Support_Figure_3_Spatial_Smooth.png", sf3_spatial_smooth, width = 10, height = 6, dpi = 300, bg = "white")

cat("\n--- [4/4] 恭喜！所有高清图表 (PDF/PNG) 均已成功生成并保存在当前工作目录中！ ---\n")