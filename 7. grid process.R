path <- "data/climate/landuse.tif"
landuse_ras <- rast(path)
summary(landuse_ras)

lc_mapping <- c(
  "0"   = "No_Data",
  "111" = "Closed_Forest_Evergreen_Needle",
  "112" = "Closed_Forest_Evergreen_Broad",
  "113" = "Closed_Forest_Deciduous_Needle",
  "114" = "Closed_Forest_Deciduous_Broad",
  "115" = "Closed_Forest_Mixed",
  "116" = "Closed_Forest_Unknown",
  "121" = "Open_Forest_Evergreen_Needle",
  "122" = "Open_Forest_Evergreen_Broad",
  "123" = "Open_Forest_Deciduous_Needle",
  "124" = "Open_Forest_Deciduous_Broad",
  "125" = "Open_Forest_Mixed",
  "126" = "Open_Forest_Unknown",
  "20"  = "Shrubs",
  "30"  = "Herbaceous_Vegetation",
  "90"  = "Herbaceous_Wetland",
  "100" = "Moss_and_Lichen",
  "60"  = "Bare_Sparse_Vegetation",
  "40"  = "Cultivated_Agriculture",
  "50"  = "Urban_Built_up",
  "70"  = "Snow_and_Ice",
  "80"  = "Permanent_Water",
  "200" = "Open_Sea"
)

env_cont <- c(elev_ras, bio15_ras, bio4_ras)
names(env_cont) <- c("elevation", "bio15", "bio4")
names(landuse_ras) <- "landuse"

extract_pixel_data_all <- function(f, sp_name, season_name, env_cont, landuse_ras) {
  tryCatch({
    r <- rast(f) 
    if (!any(is.lonlat(r))) {
      r <- project(r, "EPSG:4326", method = "near")
    }
    
    if (global(r == 1, "sum", na.rm = TRUE)[1, 1] == 0) {
      return(tibble())
    }
    
    env_cont_proj <- project(env_cont, r, method = "bilinear") 
    
    landuse_proj <- project(landuse_ras, r, method = "near")
    
    env_all_proj <- c(env_cont_proj, landuse_proj)
  
    r_mask <- r
    r_mask[r_mask != 1] <- NA 
    env_masked <- mask(env_all_proj, r_mask)
    
    pixel_df <- as.data.frame(env_masked, xy = TRUE, na.rm = TRUE)
    pixel_df <- pixel_df %>%
      mutate(
        species = sp_name,
        season = season_name
      )
    
    return(pixel_df)
    
  }, error = function(e) {
    message("\n❌ 提取失败文件: ", f, " | 错误: ", conditionMessage(e))
    return(tibble()) 
  })
}

#pixel_level_df <- pmap_dfr(
#  list(ras_meta$file, ras_meta$species, ras_meta$season),
#  ~ extract_pixel_data_all(
#    f = ..1, 
#    sp_name = ..2, 
#    season_name = ..3, 
#    env_cont = env_cont, 
#    landuse_ras = landuse_ras
#  )
#)

#lc_df <- data.frame(
#  landuse = as.numeric(names(lc_mapping)), 
#  landuse_name = unname(lc_mapping)        
#)

#pixel_level_df <- pixel_level_df %>%
#  left_join(lc_df, by = "landuse") %>%
#  mutate(landuse = landuse_name) %>% 
#  dplyr::select(-landuse_name)

#write_csv(pixel_level_df, "output/pixel_level_env_data.csv")

pixel_level_env_data <- read_csv("output/pixel_level_env_data.csv")

pixel_level_df <- pixel_level_df %>%
  filter(!landuse %in% c("No_Data", "Open_Sea")) %>%
  filter(!is.na(landuse))

table(pixel_level_df$landuse)






