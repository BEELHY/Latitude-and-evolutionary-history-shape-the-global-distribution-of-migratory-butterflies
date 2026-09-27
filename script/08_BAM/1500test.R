# Spatial richness model test with Moran bootstrap.

library(mgcv)
library(parallel)
library(spdep)
library(car)
library(ggplot2)

file_path <- "/home/hli/im26/Distribution_model/data/cleaned_data_for_SAR.csv"

df_aggressive <- read.csv(file_path)

df_aggressive$Landuse <- factor(df_aggressive$Landuse)
df_aggressive$Landuse <- relevel(df_aggressive$Landuse, ref = "cropland")

target_vars <- c("Richness", "Bio_4", "Bio_15", "Elevation", "HII", "Landuse", "lon", "lat")
df_aggressive <- na.omit(df_aggressive[, target_vars])

n_cores <- max(1, detectCores() - 1)

cat("\n--- [1/5] 正在拟合稳健版空间模型 (k=1500, bs='ds') --- \n")
cat("数据集大小:", nrow(df_aggressive), "条记录\n")

model_final_v2 <- bam(
  Richness ~ s(Bio_4, k=10) + s(Bio_15, k=10) + s(Elevation, k=10) +
    HII + Landuse +
    s(lon, lat, bs = "ds", m = c(1, 0.5), k = 1500),
  data = df_aggressive,
  family = nb(),
  method = "fREML",
  discrete = TRUE,
  nthreads = n_cores
)

summary(model_final_v2)

spatial_edf <- summary(model_final_v2)$s.table[4, "edf"]
cat("\n空间项的有效自由度 (edf) 为:", spatial_edf, "\n")
if(spatial_edf > 1450) {
  cat("警告：edf 接近 k (1500)，可能仍存在未吸收的空间自相关。\n")
} else {
  cat("提示：edf 远小于 k，说明空间平滑项捕捉效果良好。\n")
}

cat("\n--- [2/5] 正在通过代理模型评估共线性 (内存安全版) --- \n")

set.seed(123)
df_proxy <- df_aggressive[sample(1:nrow(df_aggressive), 50000), ]

m_proxy <- gam(
  Richness ~ s(Bio_4, k=10) + s(Bio_15, k=10) + s(Elevation, k=10) +
    HII + Landuse +
    s(lon, lat, bs = "ds", m = c(1, 0.5), k = 1500),
  data = df_proxy,
  family = nb()
)

con_check <- concurvity(m_proxy, full = TRUE)
cat("\n--- 共线性诊断结果 (Full Concurvity) ---\n")
print(round(con_check, 3))

cat("\n--- [3/5] 正在执行残差空间自相关检验 (100次重复抽样自助法) --- \n")

df_aggressive$all_res <- residuals(model_final_v2, type = "pearson")

n_boot <- 100
moran_estimates <- numeric(n_boot)
p_values <- numeric(n_boot)

set.seed(123)

for(i in 1:n_boot) {
  sub_idx <- sample(1:nrow(df_aggressive), 10000)
  df_sub_loop <- df_aggressive[sub_idx, ]

  coords_loop <- cbind(df_sub_loop$lon, df_sub_loop$lat)
  nb_loop <- knn2nb(knearneigh(coords_loop, k = 12, longlat = TRUE))
  lw_loop <- nb2listw(nb_loop, style = "W", zero.policy = TRUE)

  test_loop <- moran.test(df_sub_loop$all_res, lw_loop, zero.policy = TRUE)
  moran_estimates[i] <- test_loop$estimate[1]
  p_values[i] <- test_loop$p.value

  if(i %% 10 == 0) cat(paste0("进度: [", i, "/", n_boot, "] 轮计算完成...\n"))
}

mean_moran <- mean(moran_estimates)
sd_moran   <- sd(moran_estimates)
ci_moran   <- quantile(moran_estimates, probs = c(0.025, 0.975))
prop_significant <- mean(p_values < 0.05) * 100

cat("\n======================================================\n")
cat("          Moran's I 多次抽样诊断报告 (N=100)          \n")
cat("======================================================\n")
cat("Moran's I 均值       :", round(mean_moran, 4), "\n")
cat("Moran's I 标准差     :", round(sd_moran, 4), "\n")
cat("95% 置信区间 (CI)    : [", round(ci_moran[1], 4), ",", round(ci_moran[2], 4), "]\n")
cat("P值显著(<0.05)的比例 :", prop_significant, "%\n")
cat("======================================================\n")

cat("\n[数据保存] 正在将模型对象与诊断结果导出至当前工作目录...\n")

saveRDS(model_final_v2, "butterfly_bam_model.rds")

boot_results_df <- data.frame(
  Repetition = 1:n_boot,
  Moran_I = moran_estimates,
  P_Value = p_values
)
write.csv(boot_results_df, "moran_bootstrap_distribution.csv", row.names = FALSE)

df_sub <- df_sub_loop
df_sub$res <- df_sub$all_res
map_data_df <- df_sub[, c("lon", "lat", "res", "Richness")]
write.csv(map_data_df, "spatial_residuals_map_data.csv", row.names = FALSE)

cat("✅ 成功保存以下文件：\n")
cat("   - butterfly_bam_model.rds          (完整的 BAM 模型二进制对象)\n")
cat("   - moran_bootstrap_distribution.csv (100次抽样的 Moran's I 值)\n")
cat("   - spatial_residuals_map_data.csv   (用于画残差地图的经纬度点)\n")

cat("\n--- [4/5] 正在生成默认诊断图表 --- \n")

p_res <- ggplot(map_data_df, aes(x = lon, y = lat, color = res)) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  theme_minimal() +
  labs(title = "Pearson Residuals Spatial Distribution",
       subtitle = paste0("Mean Moran's I: ", round(mean_moran, 4),
                         " | 95% CI: [", round(ci_moran[1], 4), ", ", round(ci_moran[2], 4), "]"))

try(gam.check(model_final_v2))

ggsave("residuals_map.png", p_res, width=12, height=6)
cat("\n--- [5/5] 脚本运行全部结束，模型与可视化基础数据已安全落盘！ ---\n")
