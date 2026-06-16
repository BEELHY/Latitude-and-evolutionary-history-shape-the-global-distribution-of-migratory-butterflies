# ==============================================================================
# 全球蝴蝶多样性 BAM 模型 - 最终版可视化脚本 (聚光灯高亮版)
# 包含：Main Figure (生态学驱动力同框对比) & Support Figures (统计学稳健性检验)
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. 加载必备扩展包
# ------------------------------------------------------------------------------
library(mgcv)       # BAM 模型基础包
library(ggplot2)    # 绘图底层逻辑
library(gratia)     # GAM/BAM 模型专属高颜值可视化提取包
library(patchwork)  # 优雅的图表拼图神器
library(dplyr)      # 数据清洗与合并必备包

cat("\n--- [1/4] 环境就绪，正在加载模型与数据... ---\n")

# ------------------------------------------------------------------------------
# 2. 读取持久化的模型与诊断数据
# ------------------------------------------------------------------------------
# 读取核心 BAM 模型 (加载时间视内存大小而定，请耐心等待)
model_final <- readRDS("/Users/hlii0385/Desktop/Phd_Haiyu_LI/Rapoport-s-rule-and-Bergmann-s-rule-of-Migratory-butterflies/output/BAM/output/butterfly_bam_model.rds")

# 读取原始清洗数据并处理 Landuse 因子
file_path <- "/Users/hlii0385/Desktop/Phd_Haiyu_LI/Rapoport-s-rule-and-Bergmann-s-rule-of-Migratory-butterflies/output/BAM/output/cleaned_data_for_SAR.csv" 
df_aggressive <- read.csv(file_path)
df_aggressive$Landuse <- factor(df_aggressive$Landuse)
df_aggressive$Landuse <- relevel(df_aggressive$Landuse, ref = "cropland")

# 读取空间残差诊断数据
map_data <- read.csv("/Users/hlii0385/Desktop/Phd_Haiyu_LI/Rapoport-s-rule-and-Bergmann-s-rule-of-Migratory-butterflies/output/BAM/output/spatial_residuals_map_data.csv")

# 读取 Moran's I 的 100 次抽样诊断数据
boot_data <- read.csv("/Users/hlii0385/Desktop/Phd_Haiyu_LI/Rapoport-s-rule-and-Bergmann-s-rule-of-Migratory-butterflies/output/BAM/output/moran_bootstrap_distribution.csv")

cat("✅ 模型与数据加载成功！\n")


# ==============================================================================
# 3. 绘制 Main Figure：聚光灯平移图 (驱动全球蝴蝶多样性的核心机制)
# ==============================================================================
cat("\n--- [2/4] 正在生成 Main Figure (聚光灯截距平移图)... ---\n")

# ------------------------------------------------------------------------------
# A. 提取非线性平滑项 (Bio_4, Bio_15, Elevation) - 继续用 gratia，这部分它很稳定
# ------------------------------------------------------------------------------
get_smooth_curve <- function(mod, term_name, var_name) {
  dat <- smooth_estimates(mod, select = term_name)
  # 兼容不同版本的列名
  x_val <- if (".value" %in% names(dat)) dat$.value else dat[[var_name]]
  y_val <- if (".estimate" %in% names(dat)) dat$.estimate else dat$partial
  # X 轴缩放
  x_scaled <- (x_val - min(x_val)) / (max(x_val) - min(x_val))
  data.frame(Variable = var_name, X_Relative = x_scaled, Effect = y_val)
}

curve_bio4  <- get_smooth_curve(model_final, "s(Bio_4)", "Bio_4")
curve_bio15 <- get_smooth_curve(model_final, "s(Bio_15)", "Bio_15")
curve_elev  <- get_smooth_curve(model_final, "s(Elevation)", "Elevation")

# ------------------------------------------------------------------------------
# B. 提取连续线性项 (HII) - 使用最稳健的原生 mgcv predict
# ------------------------------------------------------------------------------
# 创建 100 个 HII 的等距序列 (从你的真实数据最小值到最大值)
hii_seq <- seq(min(df_aggressive$HII, na.rm = TRUE), 
               max(df_aggressive$HII, na.rm = TRUE), length.out = 100)

# 构造一个预测专用的伪数据框
pred_df_hii <- df_aggressive[1:100, ] 
pred_df_hii$HII <- hii_seq

# type = "terms" 是精髓：它会自动剥离出单一变量的纯偏效应
hii_terms <- predict(model_final, newdata = pred_df_hii, type = "terms")

curve_hii <- data.frame(
  Variable = "HII",
  X_Relative = seq(0, 1, length.out = 100), # 因为是线性，相对梯度直接就是 0 到 1
  Effect = hii_terms[, "HII"]               # 精准提取 HII 的列
)

# 合并所有连续变量的曲线
all_curves <- bind_rows(curve_bio4, curve_bio15, curve_elev, curve_hii)

# ------------------------------------------------------------------------------
# C. 提取分类线性项 (Landuse) - 同理，使用稳健的原生 predict
# ------------------------------------------------------------------------------
lu_levels <- levels(df_aggressive$Landuse)
pred_df_lu <- df_aggressive[1:length(lu_levels), ]
pred_df_lu$Landuse <- factor(lu_levels, levels = lu_levels)

# 提取分类变量的偏效应
lu_terms <- predict(model_final, newdata = pred_df_lu, type = "terms")

landuse_eff <- data.frame(
  Landuse = lu_levels,
  Shift   = lu_terms[, "Landuse"]
)

# 强制将 reference level (cropland) 的效应归 0，以此为平移基准
ref_shift <- landuse_eff$Shift[landuse_eff$Landuse == "cropland"]
landuse_eff$Shift <- landuse_eff$Shift - ref_shift

# ------------------------------------------------------------------------------
# D. 笛卡尔积连接计算所有分组的平行移动
# ------------------------------------------------------------------------------
plot_data_shifted <- cross_join(all_curves, landuse_eff) |>
  mutate(
    Shifted_Effect = Effect + Shift,
    
    # 💡 已经为你修改为高亮 urban 组
    Target_Group = "built", 
    
    Is_Target = (Landuse == Target_Group)
  )

# ------------------------------------------------------------------------------
# E. 绘制最终的主图
# ------------------------------------------------------------------------------
main_figure <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  
  # 背景层：其他所有 Landuse 的淡化线条（此时 forest, cropland 等都会变成背景灰色）
  geom_line(
    data = plot_data_shifted |> filter(!Is_Target),
    aes(x = X_Relative, y = Shifted_Effect, group = interaction(Variable, Landuse)),
    color = "grey60", alpha = 0.25, linewidth = 0.8
  ) +
  
  # 顶层高亮：只有 urban 组对应的 4 条鲜艳曲线
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
# 导出主图至你的工作目录 (基于你在运行脚本时的当前 setwd 目录)
ggsave("/Users/hlii0385/Desktop/Phd_Haiyu_LI/Rapoport-s-rule-and-Bergmann-s-rule-of-Migratory-butterflies/Figure_1_Main_Spotlight.pdf", main_figure, width = 10, height = 7, dpi = 300)
ggsave("/Users/hlii0385/Desktop/Phd_Haiyu_LI/Rapoport-s-rule-and-Bergmann-s-rule-of-Migratory-butterflies/Figure_1_Main_Spotlight.png", main_figure, width = 10, height = 7, dpi = 300, bg = "white")

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
        panel.background = element_rect(fill = "aliceblue", color = NA)) 

ggsave("/Users/hlii0385/Desktop/Phd_Haiyu_LI/Rapoport-s-rule-and-Bergmann-s-rule-of-Migratory-butterflies/Support_Figure_1_Residual_Map.png", sf1_map, width = 12, height = 6, dpi = 300, bg = "white")

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

ggsave("/Users/hlii0385/Desktop/Phd_Haiyu_LI/Rapoport-s-rule-and-Bergmann-s-rule-of-Migratory-butterflies/Support_Figure_2_Moran_Distribution.png", sf2_moran, width = 8, height = 6, dpi = 300, bg = "white")

# ----------------------------------------------------------------
# Support Figure 3: s(lon, lat) 提取的空间平滑效应基底
# ----------------------------------------------------------------
sf3_spatial_smooth <- draw(model_final, select = "s(lon,lat)") + 
  theme_minimal(base_size = 14) +
  scale_fill_viridis_c(option = "mako", name = "Spatial\nEffect") +
  labs(title = "Extracted Spatial Smooth Effect",
       subtitle = "Based on s(lon, lat, bs='ds', k=1500)",
       x = "Longitude", y = "Latitude")

ggsave("/Users/hlii0385/Desktop/Phd_Haiyu_LI/Rapoport-s-rule-and-Bergmann-s-rule-of-Migratory-butterflies/Support_Figure_3_Spatial_Smooth.png", sf3_spatial_smooth, width = 10, height = 6, dpi = 300, bg = "white")

cat("\n--- [4/4] 恭喜！所有高清图表 (PDF/PNG) 均已成功生成并保存在指定目录中！ ---\n")