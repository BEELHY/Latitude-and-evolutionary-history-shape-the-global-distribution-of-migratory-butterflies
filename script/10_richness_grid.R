# Build richness grid data.

library(ranger)

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

pixel_level_env_data <- read_csv("output/pixel_level_env_data.csv")
pixel_level_env_data_clean <- pixel_level_env_data %>%
  filter(!landuse %in% c("No_Data", "Open_Sea", "Permanent_Water")) %>%
  mutate(landuse_grouped = case_when(
    grepl("Closed_Forest", landuse) ~ "Closed_Forest",
    grepl("Open_Forest", landuse)  ~ "Open_Forest",
    landuse %in% c("Shrubs", "Herbaceous_Vegetation") ~ "Shrub_Grass",
    landuse == "Cultivated_Agriculture" ~ "Agriculture",
    landuse == "Urban_Built_up" ~ "Urban",
    TRUE ~ "Other_Natural"
  )) %>%
  mutate(landuse_grouped = as.factor(landuse_grouped))

table(pixel_level_env_data_clean$landuse_grouped)

present_data <- pixel_level_env_data_clean %>% mutate(occ = 1)

bg_mask <- test_mask <- rast("output/masks/landuse_mask_no_water.tif")

n_samples <- nrow(present_data)

bg_pool <- spatSample(bg_mask, size = 20000000, method = "random", na.rm = TRUE, xy = TRUE) %>% as_tibble()
gc()
message("📡 正在生成背景候选池...")
bg_pool <- bg_pool %>%
  mutate(xr = round(x, 5), yr = round(y, 5))

get_specific_absent <- function(sp, se, p_data, pool) {
  pres_coords <- p_data %>%
    filter(species == sp, season == se) %>%
    dplyr::select(x, y) %>%
    mutate(x = round(x, 5), y = round(y, 5)) %>%
    distinct()

  specific_absent <- pool %>%
    anti_join(pres_coords, by = c("xr" = "x", "yr" = "y")) %>%
    slice_sample(n = nrow(pres_coords)) %>%
    mutate(species = sp, season = se, occ = 0) %>%
    dplyr::select(x, y, landuse, species, season, occ)

  return(specific_absent)
}
message("🧬 正在生成专属不在点...")
all_comb_final <- present_data %>% distinct(species, season)

absent_final_raw <- map2_dfr(
  all_comb_final$species,
  all_comb_final$season,
  ~get_specific_absent(.x, .y, present_data, bg_pool)
)

message("🔍 正在点对点提取环境因子...")
env_values_final <- terra::extract(env_cont, absent_final_raw[, c("x", "y")])

present_final <- present_data %>%
  mutate(occ = 1) %>%
  dplyr::select(x, y, elevation, bio15, bio4,
                landuse = landuse_grouped, species, season, occ)

absent_final_ready <- bind_cols(absent_final_raw, env_values_final) %>%
  dplyr::select(-any_of("ID")) %>%
  mutate(
    landuse_grouped = case_when(
      landuse %in% c(111, 112, 113, 114, 115, 116) ~ "Closed_Forest",
      landuse %in% c(121, 122, 123, 124, 125, 126) ~ "Open_Forest",
      landuse %in% c(20, 30) ~ "Shrub_Grass",
      landuse == 40 ~ "Agriculture",
      landuse == 50 ~ "Urban",
      TRUE ~ "Other_Natural"
    )
  ) %>%
  mutate(occ = 0) %>%
  dplyr::select(x, y, elevation, bio15, bio4,
                landuse = landuse_grouped, species, season, occ)

full_model_data <- bind_rows(present_final, absent_final_ready) %>%
  mutate(across(c(elevation, bio15, bio4), ~as.numeric(scale(.))))
dir.create("output/final", recursive = TRUE, showWarnings = FALSE)
saveRDS(full_model_data, "output/final/full_model_data_88m.rds")
