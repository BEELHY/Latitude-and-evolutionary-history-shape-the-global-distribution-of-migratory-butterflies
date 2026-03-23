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



df_lat_sp <- final_df %>%
  left_join(df_lat_sp , by = "species")

m_combined_cli<- lmer(log10(range_km2) ~ prop_mean+abs_lat+prop_within+season+mean_bio15+mean_elev+ mean_bio4 + (1 | species), data = final_df_lat)
summary(m_combined_cli)

final_df_scaled <- final_df_lat
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
#lmer(log10(WS_U) ~ mean_elev + (1 | Family), data = final_df_scaled)


path <- "data/climate/landuse.tif"
landuse_ras <- rast(path)
summary(landuse_ras)
