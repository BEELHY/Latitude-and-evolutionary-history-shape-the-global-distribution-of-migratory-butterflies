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