# Phylogenetic mixed models.

library(tidybayes)
library(bayesplot)
library(broom.mixed)
library(GGally)

path <- "data/climate/wc2.1_2.5m_bio_15.tif"

bio15_ras <- rast(path)

path <- "data/climate/wc2.1_2.5m_elev.tif"

elev_ras <- rast(path)

extract_mechanism_metrics <- function(f, env_layer) {
  r <- rast(f)

  if (!is.lonlat(r)) r <- project(r, "EPSG:4326", method = "near")

  env_resampled <- resample(env_layer, r, method = "bilinear")

  if (global(r == 1, "sum", na.rm = TRUE)[1, 1] == 0) return(NA_real_)

  occ_env <- mask(env_resampled, r, maskvalues = 0)

  res <- global(occ_env, "mean", na.rm = TRUE)[1, 1]
  return(res)
}

mechanism_df <- ras_meta %>%
  mutate(
    mean_elev  = map_dbl(file, ~extract_mechanism_metrics(.x, elev_ras)),
    mean_bio15 = map_dbl(file, ~extract_mechanism_metrics(.x, bio15_ras)),
    mean_bio4  = map_dbl(file, ~extract_mechanism_metrics(.x, bio4_ras))
  )

final_df <- df_wb %>%
  left_join(mechanism_df %>%
              dplyr::select(species, season, mean_bio4, mean_bio15, mean_elev),
            by = c("species", "season"))

df_lat_sp_env <- final_df %>%
  left_join(df_lat_sp , by = "species")

m_combined_cli<- lmer(log10(range_km2) ~ prop_mean+prop_within+season+mean_bio15+mean_elev+ mean_bio4 + (1 | species), data = df_lat_sp_env)
summary(m_combined_cli)
final_df_scaled <- df_lat_sp_env
vars_to_scale <- c("prop_mean", "abs_lat", "prop_within", "mean_bio15", "mean_elev", "mean_bio4")

final_df_scaled[vars_to_scale] <- lapply(final_df_scaled[vars_to_scale], function(x) {
  as.numeric(datawizard::standardize(x))
})

m_combined_cli_scaled <- lmer(log10(range_km2) ~ prop_mean + abs_lat + prop_within +
                                season + mean_bio15 + mean_elev + mean_bio4 +
                                (1 | species),
                              data = final_df_scaled)

summary(m_combined_cli_scaled)
performance::check_model(m_combined_cli_scaled)

m_combined_cli_WS_U<- lmer(log10(WS_U) ~ abs_lat+mean_bio15+mean_elev+ (1 | Family), data = final_df_scaled)
m_combined_cli_WS_U1<- lmer(log10(WS_U) ~ mean_bio4+mean_bio15+mean_elev+ (1 | Family), data = final_df_scaled)
summary(m_combined_cli_WS_U1)
summary(m_combined_cli_WS_U)
m_combined_cli_WS_L<- lmer(log10(WS_L) ~ abs_lat+mean_bio15+mean_elev+ (1 | Family), data = final_df_scaled)
summary(m_combined_cli_WS_L)
m_combined_cli_WS_L1<- lmer(log10(WS_L) ~ mean_bio4+mean_bio15+mean_elev+ (1 | Family), data = final_df_scaled)
summary(m_combined_cli_WS_L1)

df_combined<- df_combined %>%
  left_join(mechanism_df %>%
              dplyr::select(species, season, mean_elev, mean_bio15),
            by = c("species", "season"))

df_combined <- df_combined %>%
  mutate(
    mean_elev_z = as.numeric(scale(mean_elev)),
    mean_bio15_z = as.numeric(scale(mean_bio15))
  )

df_model_WS_L <- df_combined %>%
  filter(!is.na(WS_L), WS_L > 0, !is.na(mean_elev), !is.na(mean_bio15)) %>%
  mutate(log_WS_L_z = as.numeric(scale(log10(WS_L))))

df_model_WS_U <- df_combined %>%
  filter(!is.na(WS_U), WS_U > 0, !is.na(mean_elev), !is.na(mean_bio15)) %>%
  mutate(log_WS_U_z = as.numeric(scale(log10(WS_U))))

m_combined_phylo_L <- brm(
  log10(range_km2) ~ mean_bio4_z + mean_bio15_z + mean_elev_z +
    prop_within_z + log_WS_L_z + season +
    (1 | gr(species_phylo, dist = "gaussian")),
  data = df_model_WS_L,
  data2 = list(species_phylo = A),
  family = gaussian(),
  prior = c(
    prior(normal(0, 1), class = "b"),
    prior(exponential(1), class = "sd"),
    prior(exponential(1), class = "sigma")
  ),
  chains = 4, iter = 6000, warmup = 2000, cores = 4,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

m_combined_phylo_U <- brm(
  log10(range_km2) ~ mean_bio4_z + mean_bio15_z + mean_elev_z +
    prop_within_z + log_WS_U_z + season +
    (1 | gr(species_phylo, dist = "gaussian")),
  data = df_model_WS_U,
  data2 = list(species_phylo = A),
  family = gaussian(),
  prior = c(
    prior(normal(0, 1), class = "b"),
    prior(exponential(1), class = "sd"),
    prior(exponential(1), class = "sigma")
  ),
  chains = 4, iter = 6000, warmup = 2000, cores = 4,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

m_combined_phylo <- brm(
  log10(range_km2) ~ mean_bio4_z + mean_bio15_z + mean_elev_z +
    prop_within_z + season +
    (1 | gr(species_phylo, dist = "gaussian")),
  data = df_model_WS_L,
  data2 = list(species_phylo = A),
  family = gaussian(),
  prior = c(
    prior(normal(0, 1), class = "b"),
    prior(exponential(1), class = "sd"),
    prior(exponential(1), class = "sigma")
  ),
  chains = 4, iter = 6000, warmup = 2000, cores = 4,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

summary(m_combined_phylo_L)
summary(m_combined_phylo_U)
summary(m_combined_phylo)

posterior <- as.array(m_combined_phylo)
mcmc_pairs(posterior, pars = c("b_Intercept", "b_mean_bio4_z", "b_mean_bio15_z", "sigma"))

m_subset_glmm <- brm(
  log10(range_km2) ~ mean_bio4_z + mean_bio15_z + mean_elev_z +
    prop_within_z + season +
    (1 | species_phylo),
  data = df_model_WS_L,
  family = gaussian(),
  prior = c(
    prior(normal(0, 1), class = "b"),
    prior(exponential(1), class = "sd"),
    prior(exponential(1), class = "sigma")
  ),
  chains = 4, iter = 6000, warmup = 2000, cores = 4,
  control = list(adapt_delta = 0.99)
)

summary(m_subset_glmm)

m_full_glmm <- lmer(log10(range_km2) ~ prop_within +
                                season + mean_bio15 + mean_elev + mean_bio4 +
                                (1 | species),
                              data = final_df_scaled)
d1 <- tidy(m_full_glmm, conf.int = TRUE) %>% mutate(model = "Full GLMM (N=1377)")
d2 <- tidy(m_subset_glmm, conf.int = TRUE) %>% mutate(model = "Subset GLMM (N=783)")
d3 <- tidy(m_combined_phylo, conf.int = TRUE) %>% mutate(model = "Subset PGLMM (N=783)")

process_model_data <- function(df) {
  df %>%
    mutate(term = str_remove(term, "^b_")) %>%
    mutate(term = str_remove(term, "_z$")) %>%
    filter(!term %in% c("(Intercept)", "Intercept")) %>%
    mutate(
      clean_term = case_when(
        term == "mean_bio4"   ~ "Temp. Seasonality (Bio4)",
        term == "mean_bio15"  ~ "Precip. Seasonality (Bio15)",
        term == "mean_elev"   ~ "Mean Elevation",
        term == "prop_within" ~ "Tropical Proportion",
        term == "prop_mean"   ~ "Mean Habitat Prop.",
        term == "seasonS2"    ~ "Season S2",
        term == "seasonS3"    ~ "Season S3",
        term == "seasonS4"    ~ "Season S4",
        TRUE ~ term
      ),
      is_significant = ifelse(conf.low > 0 | conf.high < 0, "Significant", "Non-significant")
    )
}

plot_df <- bind_rows(d1, d2, d3) %>% process_model_data()

ggplot(plot_df, aes(y = reorder(clean_term, estimate), x = estimate,
                    color = model, group = model)) +
  geom_vline(xintercept = 0, color = "gray50", linetype = "dashed", size = 0.6) +

  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 position = position_dodge(width = 0.7),
                 height = 0.25, size = 0.8) +

  geom_point(aes(shape = is_significant),
             position = position_dodge(width = 0.7),
             size = 3.5) +

  scale_color_manual(values = c(
    "Full GLMM (N=1377)" = "#16A085",
    "Subset GLMM (N=783)" = "#2E86C1",
    "Subset PGLMM (N=783)" = "#E67E22"
  )) +
  scale_shape_manual(values = c("Significant" = 16, "Non-significant" = 1)) +

  labs(
    title = "Sensitivity Analysis: Environmental Drivers of Range Size",
    subtitle = "Comparing Full Data vs. Phylogenetic Subset Models",
    x = "Standardized Coefficient Estimate (95% CI/CrI)",
    y = NULL,
    color = "Model Version",
    shape = "Significance"
  ) +

  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    axis.text.y = element_text(face = "bold", color = "black"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14)
  )

df_full_corr <- final_df_scaled %>%
  dplyr::select(mean_bio4, mean_bio15, mean_elev, prop_within)

df_subset_corr <- df_model_WS_L %>%
  dplyr::select(mean_bio4_z, mean_bio15_z, mean_elev_z, prop_within_z) %>%
  rename_with(~stringr::str_remove(., "_z$"))

ggpairs(df_full_corr,
        title = "Correlation Matrix: Full Dataset (N=1377)",
        upper = list(continuous = wrap("cor", size = 4, color = "black")),
        diag = list(continuous = wrap("densityDiag", fill = "#16A085", alpha = 0.5)),
        lower = list(continuous = wrap("smooth", alpha = 0.1, size = 0.1, color = "#16A085"))) +
  theme_bw()

ggpairs(df_subset_corr,
        title = "Correlation Matrix: Subset (N=783)",
        upper = list(continuous = wrap("cor", size = 4, color = "black")),
        diag = list(continuous = wrap("densityDiag", fill = "#E67E22", alpha = 0.5)),
        lower = list(continuous = wrap("smooth", alpha = 0.1, size = 0.1, color = "#E67E22"))) +
  theme_bw()

df <- data.frame(
  Variable = factor(c('Grassland', 'Shrubs', 'Trees/Forest', 'Cropland', 'Bio_4', 'Bio_15', 'Elevation', 'HII'),
                    levels = rev(c('Grassland', 'Shrubs', 'Trees/Forest', 'Cropland', 'Bio_4', 'Bio_15', 'Elevation', 'HII'))),
  Model = rep(c("OLS", "GLM"), each = 8),
  Estimate = c(-3.311, -3.117, -2.408, -1.753, -0.689, -0.310, 0.036, 0.136,
               -1.786, -1.188, -0.892, -0.396, -0.648, -0.454, 0.053, 0.137),
  Lower = c(-3.463, -3.272, -2.567, -1.909, -0.696, -0.321, 0.026, 0.130,
            -1.895, -1.317, -1.000, -0.525, -0.655, -0.466, 0.039, 0.129),
  Upper = c(-3.111, -2.925, -2.201, -1.561, -0.682, -0.299, 0.045, 0.141,
            -1.668, -1.035, -0.766, -0.259, -0.641, -0.442, 0.066, 0.146)
)

df <- data.frame(
  Variable = factor(c('Grassland', 'Shrubs', 'Trees/Forest', 'Cropland', 'Bio_4', 'Bio_15', 'Elevation', 'HII'),
                    levels = rev(c('Grassland', 'Shrubs', 'Trees/Forest', 'Cropland', 'Bio_4', 'Bio_15', 'Elevation', 'HII'))),
  Model = rep(c("OLS", "GLM"), each = 8),
  Estimate = c(-3.311, -3.117, -2.408, -1.753, -0.689, -0.310, 0.036, 0.136,
               -1.786, -1.188, -0.892, -0.396, -0.648, -0.454, 0.053, 0.137),
  Lower = c(-3.463, -3.272, -2.567, -1.909, -0.696, -0.321, 0.026, 0.130,
            -1.895, -1.317, -1.000, -0.525, -0.655, -0.466, 0.039, 0.129),
  Upper = c(-3.111, -2.925, -2.201, -1.561, -0.682, -0.299, 0.045, 0.141,
            -1.668, -1.035, -0.766, -0.259, -0.641, -0.442, 0.066, 0.146)
)

pd <- position_dodge(width = 0.6)

ggplot(df, aes(x = Estimate, y = Variable, color = Model, fill = Model)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +

  geom_point(aes(shape = Model),
             size = 2.5,
             position = pd,
             alpha = 0.7) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper),
                 height = 0.4,
                 linewidth = 0.8,
                 position = pd) +

  scale_color_manual(values = c("OLS" = "#56B4E9", "GLM" = "#D55E00")) +
  theme_bw() +
  theme(legend.position = "top",
        panel.grid.minor = element_blank()) +
  labs(title = "Model Comparison with Visible Error Bars",
       x = "Estimate (95% CI)")

m_no_phylo <- m_no_phylo_bayes <- brm(
  formula = log10(range_km2) ~ mean_bio4_z + mean_bio15_z + mean_elev_z + prop_within_z + season,
  data = df_model_WS_L,
  family = gaussian(),
  chains = 4,
  iter = 2000,
  cores = 4
)

summary(m_no_phylo)

m_no_phylo_all <- m_no_phylo_bayes <- brm(
  formula = log10(range_km2) ~ mean_bio4 + mean_bio15 + mean_elev + prop_within + season,
  data = final_df,
  family = gaussian(),
  chains = 4,
  iter = 2000,
  cores = 4
)

summary(m_no_phylo_all)

loo_compare(loo(m_combined_phylo), loo(m_no_phylo))
