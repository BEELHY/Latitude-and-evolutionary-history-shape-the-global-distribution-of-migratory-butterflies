# Spatial richness GAM.

# Libraries
library(mgcv)
library(parallel)
library(spdep)  # Spatial diagnostics
library(car)  # VIF
library(ggplot2)  # Plots

# Load and clean data
file_path <- "updatedata/richness_grid.csv"
out_dir <- "output/BAM/output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

df_aggressive <- read.csv(file_path)

# Factors, reference level
df_aggressive$Landuse <- factor(df_aggressive$Landuse)
df_aggressive$Landuse <- relevel(df_aggressive$Landuse, ref = "cropland")

# Drop missing values
target_vars <- c("Richness", "Bio_4", "Bio_15", "Elevation", "HII", "Landuse", "lon", "lat")
df_aggressive <- na.omit(df_aggressive[, target_vars])

# Parallel cores
n_cores <- max(1, detectCores() - 1)

# Fit spatial GAM
cat("\n--- [1/5] Fitting spatial GAM --- \n")
cat("Rows:", nrow(df_aggressive), "\n")

# bam with discrete fREML
model_final_v2 <- bam(
  Richness ~ s(Bio_4, k=10) + s(Bio_15, k=10) + s(Elevation, k=10) +
    HII + Landuse +
    s(lon, lat, bs = "ds", m = c(1, 0.5), k = 1500),
  data = df_aggressive,
  family = nb(),
  method = "fREML",
  discrete = TRUE,  # Speeds large data
  nthreads = n_cores  # Parallel threads
)

# Model summary
summary(model_final_v2)

# Check edf
spatial_edf <- summary(model_final_v2)$s.table[4, "edf"]
cat("\nSpatial smooth edf:", spatial_edf, "\n")
if(spatial_edf > 1450) {
  cat("Warning: edf close to k.\n")
} else {
  cat("edf well below k.\n")
}

# Collinearity on subsample
cat("\n--- [2/5] Concurvity check --- \n")

set.seed(123)
# Sample 50,000 cells
df_proxy <- df_aggressive[sample(1:nrow(df_aggressive), 50000), ]

# Proxy GAM on subset
m_proxy <- gam(
  Richness ~ s(Bio_4, k=10) + s(Bio_15, k=10) + s(Elevation, k=10) +
    HII + Landuse +
    s(lon, lat, bs = "ds", m = c(1, 0.5), k = 1500),
  data = df_proxy,
  family = nb()
)

# VIF
con_check <- concurvity(m_proxy, full = TRUE)
cat("\n--- Full concurvity ---\n")
print(round(con_check, 3))

# Moran's I, save outputs
cat("\n--- [3/5] Residual Moran's I (100 subsamples) --- \n")

# Pearson residuals
df_aggressive$all_res <- residuals(model_final_v2, type = "pearson")

# Resampling settings
n_boot <- 100
moran_estimates <- numeric(n_boot)
p_values <- numeric(n_boot)

set.seed(123)  # Reproducible sampling

# Repeated Moran's I
for(i in 1:n_boot) {
  sub_idx <- sample(1:nrow(df_aggressive), 10000)
  df_sub_loop <- df_aggressive[sub_idx, ]

  coords_loop <- cbind(df_sub_loop$lon, df_sub_loop$lat)
  nb_loop <- knn2nb(knearneigh(coords_loop, k = 12, longlat = TRUE))
  lw_loop <- nb2listw(nb_loop, style = "W", zero.policy = TRUE)

  test_loop <- moran.test(df_sub_loop$all_res, lw_loop, zero.policy = TRUE)
  moran_estimates[i] <- test_loop$estimate[1]
  p_values[i] <- test_loop$p.value

  if(i %% 10 == 0) cat(paste0("Progress: [", i, "/", n_boot, "]\n"))
}

# Summary statistics
mean_moran <- mean(moran_estimates)
sd_moran   <- sd(moran_estimates)
ci_moran   <- quantile(moran_estimates, probs = c(0.025, 0.975))  # 95% interval
prop_significant <- mean(p_values < 0.05) * 100

cat("\n======================================================\n")
cat("Moran's I over 100 subsamples\n")
cat("======================================================\n")
cat("Mean Moran's I :", round(mean_moran, 4), "\n")
cat("SD Moran's I   :", round(sd_moran, 4), "\n")
cat("95% CI         : [", round(ci_moran[1], 4), ",", round(ci_moran[2], 4), "]\n")
cat("% with P < 0.05:", prop_significant, "%\n")
cat("======================================================\n")

# Save model and diagnostics
cat("\nSaving model and diagnostics...\n")

# Model object
saveRDS(model_final_v2, file.path(out_dir, "butterfly_bam_model.rds"))

# Moran's I draws
boot_results_df <- data.frame(
  Repetition = 1:n_boot,
  Moran_I = moran_estimates,
  P_Value = p_values
)
write.csv(boot_results_df, file.path(out_dir, "moran_bootstrap_distribution.csv"), row.names = FALSE)

# Residual map data
df_sub <- df_sub_loop
df_sub$res <- df_sub$all_res
map_data_df <- df_sub[, c("lon", "lat", "res", "Richness")]
write.csv(map_data_df, file.path(out_dir, "spatial_residuals_map_data.csv"), row.names = FALSE)

cat("Saved:\n")
cat("   - butterfly_bam_model.rds\n")
cat("   - moran_bootstrap_distribution.csv\n")
cat("   - spatial_residuals_map_data.csv\n")

# Diagnostic plots
cat("\n--- [4/5] Diagnostic plots --- \n")

# Residual map
p_res <- ggplot(map_data_df, aes(x = lon, y = lat, color = res)) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  theme_minimal() +
  labs(title = "Pearson Residuals Spatial Distribution",
       subtitle = paste0("Mean Moran's I: ", round(mean_moran, 4),
                         " | 95% CI: [", round(ci_moran[1], 4), ", ", round(ci_moran[2], 4), "]"))

# Standard GAM checks
try(gam.check(model_final_v2))

ggsave(file.path(out_dir, "residuals_map.png"), p_res, width=12, height=6)
cat("\n--- [5/5] Done ---\n")
