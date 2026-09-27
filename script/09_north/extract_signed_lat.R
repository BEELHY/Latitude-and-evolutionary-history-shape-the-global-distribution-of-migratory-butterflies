# Extract signed centroid latitude per species.

setwd("/Users/hlii0385/Desktop/Phd_Haiyu_LI/Rapoport-s-rule-and-Bergmann-s-rule-of-Migratory-butterflies")
suppressPackageStartupMessages({ library(terra); library(dplyr); library(stringr) })

cat("Extracting signed centroid latitude from rasters (script 04 logic)...\n")

path <- "data/SuitabilityMaps_MigratorySpecies"
ras_files <- list.files(path, pattern = "^Binary_S[1-4].*[.]tif$", full.names = TRUE)
cat("Total raster files:", length(ras_files), "\n")

ras_meta <- tibble(file = ras_files) %>%
  mutate(
    stem    = tools::file_path_sans_ext(basename(file)),
    season  = str_extract(stem, "S[1-4]"),
    species = str_replace(stem, "^Binary_S[1-4]", "")
  )

df_lat_sp <- read.csv("data/df_lat_sp.csv") %>%
  mutate(species = str_replace_all(species, " ", "_"))
cat("Species in df_lat_sp:", nrow(df_lat_sp), "\n")

one_per_sp <- ras_meta %>%
  arrange(species, season) %>%
  group_by(species) %>% slice(1) %>% ungroup() %>%
  filter(species %in% df_lat_sp$species)
cat("Species with raster:", nrow(one_per_sp), "\n")

get_signed_lat <- function(f) {
  tryCatch({
    r <- rast(f)
    if (!is.lonlat(r)) r <- project(r, "EPSG:4326", method = "near")
    lat    <- init(r, "y")
    masked <- mask(lat, r, maskvalues = 0)
    as.numeric(global(masked, "mean", na.rm = TRUE)[1, 1])
  }, error = function(e) NA_real_)
}

cat("Computing signed lat...\n")
signed_lats <- numeric(nrow(one_per_sp))
for (i in seq_len(nrow(one_per_sp))) {
  if (i %% 50 == 0 || i == 1) cat(" progress:", i, "/", nrow(one_per_sp), "\n")
  signed_lats[i] <- get_signed_lat(one_per_sp$file[i])
}

result <- one_per_sp %>%
  mutate(centroid_latitude = signed_lats) %>%
  dplyr::select(species, centroid_latitude)

write.csv(result, "output/north/signed_lat_per_species.csv", row.names = FALSE)
cat("Saved: output/north/signed_lat_per_species.csv\n")
cat("N with signed lat:", sum(!is.na(result$centroid_latitude)), "\n")
cat("Preview:\n")
print(head(result, 10))
