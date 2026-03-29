#use elevation and precipitation to see new model

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

m_combined_cli<- lmer(log10(range_km2) ~ prop_mean+abs_lat+prop_within+season+mean_bio15+mean_elev+ mean_bio4 + (1 | species), data = df_lat_sp_env)
summary(m_combined_cli)

final_df_scaled <- df_lat_sp_env
vars_to_scale <- c("prop_mean", "abs_lat", "prop_within", "mean_bio15", "mean_elev", "mean_bio4")

final_df_scaled[vars_to_scale] <- lapply(final_df_scaled[vars_to_scale], scale)

m_combined_cli_scaled <- lmer(log10(range_km2) ~ prop_mean + abs_lat + prop_within + 
                                season + mean_bio15 + mean_elev + mean_bio4 + 
                                (1 | species), 
                              data = final_df_scaled)

summary(m_combined_cli_scaled)

m_combined_cli_WS_U<- lmer(log10(WS_U) ~ prop_mean+abs_lat+prop_within+season+mean_bio15+mean_elev+ mean_bio4 + (1 | Family), data = final_df_scaled)
summary(m_combined_cli_WS_U)
#lmer(log10(WS_U) ~ mean_bio15+mean_elev + (1 | Family), data = final_df_scaled)
m_combined_cli_WS_L<- lmer(log10(WS_L) ~ prop_mean+abs_lat+prop_within+season+mean_bio15+mean_elev+ mean_bio4 + (1 | Family), data = final_df_scaled)
summary(m_combined_cli_WS_L)
#lmer(log10(WS_L) ~ mean_elev + (1 | Family), data = final_df_scaled)

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

summary(m_combined_phylo_L)
summary(m_combined_phylo_U)


