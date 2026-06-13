# ==============================================================================
# 最终整合版：全球蝴蝶多样性 BAM 模型拟合与全套诊断
# 针对 27 万行大数据集优化 (k=1500, bs='ds')
# ==============================================================================

# 1. 环境准备
library(mgcv)
library(parallel)
library(spdep)  # 空间诊断
library(car)    # VIF诊断
library(ggplot2) # 绘图

# 2. 数据加载与预处理
file_path <- "/home/hli/im26/Distribution_model/data/cleaned_data_for_SAR.csv" 

df_aggressive <- read.csv(file_path)

# 预处理：因子化与缺失值处理 (严格对照你原本的代码)
df_aggressive$Landuse <- factor(df_aggressive$Landuse)
df_aggressive$Landuse <- relevel(df_aggressive$Landuse, ref = "cropland")

# 确保所有使用的变量都没有 NA
target_vars <- c("Richness", "Bio_4", "Bio_15", "Elevation", "HII", "Landuse", "lon", "lat")
df_aggressive <- na.omit(df_aggressive[, target_vars])

# 3. 配置多核运算
n_cores <- max(1, detectCores() - 1)

# ------------------------------------------------------------------------------
# 4. 拟合稳健版空间模型 (k=1500, bs='ds')
# ------------------------------------------------------------------------------
cat("\n--- [1/5] 正在拟合稳健版空间模型 (k=1500, bs='ds') --- \n")
cat("数据集大小:", nrow(df_aggressive), "条记录\n")

# 使用 bam 配合 fREML 和 discrete=TRUE 是处理 27万个点的唯一可行方案
model_final_v2 <- bam(
  Richness ~ s(Bio_4, k=10) + s(Bio_15, k=10) + s(Elevation, k=10) + 
    HII + Landuse + 
    s(lon, lat, bs = "ds", m = c(1, 0.5), k = 1500), 
  data = df_aggressive, 
  family = nb(),       
  method = "fREML",    
  discrete = TRUE,     # 加速大数据计算的关键
  nthreads = n_cores   # 多核并行
)

# 输出基础模型结果
summary(model_final_v2)

# 检查 edf
spatial_edf <- summary(model_final_v2)$s.table[4, "edf"]
cat("\n空间项的有效自由度 (edf) 为:", spatial_edf, "\n")
if(spatial_edf > 1450) {
  cat("警告：edf 接近 k (1500)，可能仍存在未吸收的空间自相关。\n")
} else {
  cat("提示：edf 远小于 k，说明空间平滑项捕捉效果良好。\n")
}


# ------------------------------------------------------------------------------
# 5. 内存避雷版：共线性诊断 (使用代理模型)
# ------------------------------------------------------------------------------
cat("\n--- [2/5] 正在通过代理模型评估共线性 (内存安全版) --- \n")

set.seed(123)
# 抽取 50,000 个点作为代表
df_proxy <- df_aggressive[sample(1:nrow(df_aggressive), 50000), ]

# 在子集上拟合一个一模一样的简易模型 (用 gam 即可，速度快)
m_proxy <- gam(
  Richness ~ s(Bio_4, k=10) + s(Bio_15, k=10) + s(Elevation, k=10) + 
    HII + Landuse + 
    s(lon, lat, bs = "ds", m = c(1, 0.5), k = 1500), 
  data = df_proxy, 
  family = nb()
)

# 现在可以安全地计算共线性了
con_check <- concurvity(m_proxy, full = TRUE)
cat("\n--- 共线性诊断结果 (Full Concurvity) ---\n")
print(round(con_check, 3))

# ------------------------------------------------------------------------------
# 6. 空间自相关诊断 (Moran's I 多次抽样) 与模型持久化保存 (策略一核心)
# ------------------------------------------------------------------------------
cat("\n--- [3/5] 正在执行残差空间自相关检验 (100次重复抽样自助法) --- \n")

# 1. 一次性提取全量模型的 Pearson 残差
df_aggressive$all_res <- residuals(model_final_v2, type = "pearson")

# 2. 配置重复抽样参数
n_boot <- 100         
moran_estimates <- numeric(n_boot)  
p_values <- numeric(n_boot)         

set.seed(123) # 保证抽样结果可重复

# 3. 执行循环抽样诊断
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

# 4. 计算统计学指标
mean_moran <- mean(moran_estimates)
sd_moran   <- sd(moran_estimates)
ci_moran   <- quantile(moran_estimates, probs = c(0.025, 0.975)) # 95% 置信区间
prop_significant <- mean(p_values < 0.05) * 100

cat("\n======================================================\n")
cat("          Moran's I 多次抽样诊断报告 (N=100)          \n")
cat("======================================================\n")
cat("Moran's I 均值       :", round(mean_moran, 4), "\n")
cat("Moran's I 标准差     :", round(sd_moran, 4), "\n")
cat("95% 置信区间 (CI)    : [", round(ci_moran[1], 4), ",", round(ci_moran[2], 4), "]\n")
cat("P值显著(<0.05)的比例 :", prop_significant, "%\n")
cat("======================================================\n")


# ==============================================================================
# 【策略一核心】持久化保存：模型对象与诊断数据一同导出
# ==============================================================================
cat("\n[数据保存] 正在将模型对象与诊断结果导出至当前工作目录...\n")

# A. 保存完整模型对象（后续随时随地加载、提取任意因子效应并绘图）
saveRDS(model_final_v2, "butterfly_bam_model.rds")

# B. 保存 Moran's I 分布数据（用于后续绘制置信区间直方图/密度图）
boot_results_df <- data.frame(
  Repetition = 1:n_boot,
  Moran_I = moran_estimates,
  P_Value = p_values
)
write.csv(boot_results_df, "moran_bootstrap_distribution.csv", row.names = FALSE)

# C. 保存地理空间残差数据（用于后续绘制精细的全球残差分布图）
df_sub <- df_sub_loop
df_sub$res <- df_sub$all_res
map_data_df <- df_sub[, c("lon", "lat", "res", "Richness")]
write.csv(map_data_df, "spatial_residuals_map_data.csv", row.names = FALSE)

cat("✅ 成功保存以下文件：\n")
cat("   - butterfly_bam_model.rds          (完整的 BAM 模型二进制对象)\n")
cat("   - moran_bootstrap_distribution.csv (100次抽样的 Moran's I 值)\n")
cat("   - spatial_residuals_map_data.csv   (用于画残差地图的经纬度点)\n")


# ------------------------------------------------------------------------------
# 7. 绘图诊断 (基于当前导出的轻量数据做图)
# ------------------------------------------------------------------------------
cat("\n--- [4/5] 正在生成默认诊断图表 --- \n")

# 空间残差图 (副标题动态展示 100 次抽样后的均值与置信区间)
p_res <- ggplot(map_data_df, aes(x = lon, y = lat, color = res)) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  theme_minimal() +
  labs(title = "Pearson Residuals Spatial Distribution",
       subtitle = paste0("Mean Moran's I: ", round(mean_moran, 4), 
                         " | 95% CI: [", round(ci_moran[1], 4), ", ", round(ci_moran[2], 4), "]"))

# 全量模型的标准诊断（如果内存不够，try 会自动跳过，保证程序能顺利走完）
try(gam.check(model_final_v2))

ggsave("residuals_map.png", p_res, width=12, height=6)
cat("\n--- [5/5] 脚本运行全部结束，模型与可视化基础数据已安全落盘！ ---\n")