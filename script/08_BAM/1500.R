# Spatial richness model on HPC.

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

cat("\n--- [3/5] 正在执行残差空间自相关检验 (Moran's I) --- \n")

set.seed(123)
sub_idx <- sample(1:nrow(df_aggressive), 10000)
df_sub <- df_aggressive[sub_idx, ]
df_sub$res <- residuals(model_final_v2, type = "pearson")[sub_idx]

coords <- cbind(df_sub$lon, df_sub$lat)
nb <- knn2nb(knearneigh(coords, k = 12, longlat = TRUE))
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

moran_test <- moran.test(df_sub$res, lw, zero.policy = TRUE)
print(moran_test)

cat("\n--- [4/5] 正在生成诊断图表 --- \n")

p_res <- ggplot(df_sub, aes(x = lon, y = lat, color = res)) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  theme_minimal() +
  labs(title = "Pearson Residuals (Full Model Resampled)",
       subtitle = paste("Moran's I:", round(moran_test$estimate[1], 4)))

try(gam.check(model_final_v2))

ggsave("residuals_map.png", p_res, width=12, height=6)
