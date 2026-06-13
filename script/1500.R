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
# 6. 空间自相关诊断 (Moran's I) - 修正版
# ------------------------------------------------------------------------------
cat("\n--- [3/5] 正在执行残差空间自相关检验 (Moran's I) --- \n")

# 注意：残差要从全量模型 model_final_v2 中提取，保证真实性
set.seed(123)
sub_idx <- sample(1:nrow(df_aggressive), 10000)
df_sub <- df_aggressive[sub_idx, ]
# 提取全量模型的残差
df_sub$res <- residuals(model_final_v2, type = "pearson")[sub_idx]

# 构建权重矩阵
coords <- cbind(df_sub$lon, df_sub$lat)
nb <- knn2nb(knearneigh(coords, k = 12, longlat = TRUE))
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

moran_test <- moran.test(df_sub$res, lw, zero.policy = TRUE)
print(moran_test)

# ------------------------------------------------------------------------------
# 7. 绘图诊断
# ------------------------------------------------------------------------------
cat("\n--- [4/5] 正在生成诊断图表 --- \n")

# A. 空间残差图
p_res <- ggplot(df_sub, aes(x = lon, y = lat, color = res)) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  theme_minimal() +
  labs(title = "Pearson Residuals (Full Model Resampled)",
       subtitle = paste("Moran's I:", round(moran_test$estimate[1], 4)))

# B. 标准诊断图 (直接针对全量模型)
# 如果 gam.check 也内存报错，就跳过它，看 summary 就够了
try(gam.check(model_final_v2))

ggsave("residuals_map.png", p_res, width=12, height=6)