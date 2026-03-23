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


path <- "data/climate/landuse.tiff"

landuse_ras <- rast(path)


