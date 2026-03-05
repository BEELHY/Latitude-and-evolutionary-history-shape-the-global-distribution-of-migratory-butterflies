# Load libraries
library(terra)
library(dplyr)
library(stringr)
library(tidyr)
library(purrr)
library(ggplot2)
library(lme4)
library(lmerTest)

# Import data
trait_range <- read_csv("output/trait_range.csv")

# Parse raster filenames
path <- "data/SuitabilityMaps_MigratorySpecies"

ras_files <- list.files(path, pattern = "^Binary_S[1-4].*", full.names = TRUE, recursive = TRUE)

ras_meta <- tibble(file = ras_files) %>%
  mutate(stem = tools::file_path_sans_ext(basename(file))) %>%
  mutate(
    season  = str_extract(stem, "S[1-4]"),
    species = str_replace(stem, "^Binary_S[1-4]", "")
  ) %>%
  dplyr::select(file, species, season)

# Latitude metrics from one raster
lat_metrics_one <- function(f) {
  r <- rast(f)
  
  # Ensure lon/lat so init(r, "y") is latitude in degrees
  if (!is.lonlat(r)) r <- project(r, "EPSG:4326", method = "near")
  
  occ <- (r == 1)
  if (global(occ, "sum", na.rm = TRUE)[1, 1] == 0) {
    return(tibble(mean_lat = NA_real_, mean_abs_lat = NA_real_,
                  lat_min = NA_real_, lat_max = NA_real_, lat_span = NA_real_))
  }
  
  lat <- init(r, "y")
  a   <- cellSize(r, unit = "km")  # area weights
  
  # Total occupied area
  area_total <- global(mask(a, occ, maskvalues = 0), "sum", na.rm = TRUE)[1, 1]
  
  # Area-weighted means
  mean_lat     <- global(mask(lat * a, occ, maskvalues = 0), "sum", na.rm = TRUE)[1, 1] / area_total
  mean_abs_lat <- global(mask(abs(lat) * a, occ, maskvalues = 0), "sum", na.rm = TRUE)[1, 1] / area_total
  
  # Range in latitude
  lat_min <- global(mask(lat, occ, maskvalues = 0), "min", na.rm = TRUE)[1, 1]
  lat_max <- global(mask(lat, occ, maskvalues = 0), "max", na.rm = TRUE)[1, 1]
  
  tibble(
    mean_lat     = mean_lat,
    mean_abs_lat = mean_abs_lat,
    lat_min      = lat_min,
    lat_max      = lat_max,
    lat_span     = lat_max - lat_min
  )
}

# Compute latitude metrics for all rasters + join to trait table
lat_df <- ras_meta %>%
  mutate(metrics = map(file, lat_metrics_one)) %>%
  unnest(metrics) %>%
  dplyr::select(species, season, mean_lat, mean_abs_lat, lat_min, lat_max, lat_span)

trait_range2 <- trait_range %>%
  mutate(species_key = str_replace_all(species, " ", "_")) %>%
  left_join(
    lat_df %>% rename(species_key = species),
    by = c("species_key", "season")
  ) %>%
  dplyr::select(-species_key)

# Species-level latitude (abs) + traits (WS_L, WS_U)
df_lat_sp <- trait_range2 %>%
  group_by(species, Family) %>%
  summarise(
    abs_lat = mean(mean_abs_lat, na.rm = TRUE),
    WS_L    = first(WS_L[!is.na(WS_L)]),
    WS_U    = first(WS_U[!is.na(WS_U)]),
    .groups = "drop"
  ) %>%
  filter(is.finite(abs_lat))

# Fit mixed models (Family random intercept)
m_lat_L <- lmer(log10(WS_L) ~ abs_lat + (1 | Family),
                data = df_lat_sp %>% filter(is.finite(WS_L), WS_L > 0))

m_lat_U <- lmer(log10(WS_U) ~ abs_lat + (1 | Family),
                data = df_lat_sp %>% filter(is.finite(WS_U), WS_U > 0))

summary(m_lat_L)
summary(m_lat_U)

# Plot (WS_U shown; swap WS_L if needed)
p_lat_U <- ggplot(df_lat_sp %>% filter(is.finite(WS_U), WS_U > 0),
                  aes(abs_lat, log10(WS_U))) +
  geom_point(aes(shape = Family), alpha = 0.5, size = 1.6) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_classic() +
  labs(x = "Mean absolute latitude of occupied cells (°)",
       y = expression(log[10]*"(WS_U)"))

p_lat_U

# Effect size helper: % change per delta degrees
effect_per_deg <- function(mod, term = "abs_lat", delta = 10) {
  b  <- fixef(mod)[term]
  se <- summary(mod)$coefficients[term, "Std. Error"]
  
  b_ci <- b + c(-1, 1) * 1.96 * se
  
  mult    <- 10^(b * delta)
  mult_ci <- 10^(b_ci * delta)
  
  pct    <- (mult - 1) * 100
  pct_ci <- (mult_ci - 1) * 100
  
  tibble(
    term = term,
    delta_deg = delta,
    multiplier = mult,
    pct_change = pct,
    pct_lwr = min(pct_ci),
    pct_upr = max(pct_ci)
  )
}

effect_per_deg(m_lat_L, delta = 10)
effect_per_deg(m_lat_U, delta = 10)

#### add in rapport
df_lat_sp <- df_lat_sp %>%
  mutate(species = str_replace_all(species, " ", "_"))

df_combined <- df_phylo_scaled %>%
  left_join(df_lat_sp %>% dplyr::select(species, WS_L, WS_U), by = "species")

df_model_WS_L <- df_combined %>%
  filter(!is.na(WS_L), WS_L > 0) %>%
  mutate(
    log_WS_L_z = as.numeric(scale(log10(WS_L))) 
  )

df_model_WS_U <- df_combined %>%
  filter(!is.na(WS_U), WS_U > 0) %>%
  mutate(
    log_WS_U_z = as.numeric(scale(log10(WS_U))) 
  )

cat("WS_L:", length(unique(df_model_WS_L$species)), "\n")
cat("WS_U:", length(unique(df_model_WS_U$species)), "\n")

m_rapoport_with_size_L <- brm(
  log10(range_km2) ~ mean_bio4_z + prop_within_z + log_WS_L_z + season + 
    (1 | gr(species_phylo, dist = "gaussian")),
  data = df_model_WS_L,
  data2 = list(species_phylo = A),
  family = gaussian(),
  prior = c(
    prior(normal(0, 1), class = "b"),         
    prior(exponential(1), class = "sd"),     
    prior(exponential(1), class = "sigma") 
  ),
  chains = 4, 
  iter = 6000,     
  warmup = 2000,   
  cores = 4,
  control = list(
    adapt_delta = 0.99,       
    max_treedepth = 15        
  )
)

m_rapoport_with_size_U <- brm(
  log10(range_km2) ~ mean_bio4_z + prop_within_z + log_WS_U_z + season + 
    (1 | gr(species_phylo, dist = "gaussian")),
  data = df_model_WS_U,
  data2 = list(species_phylo = A),
  family = gaussian(),
  prior = c(
    prior(normal(0, 1), class = "b"),         
    prior(exponential(1), class = "sd"),     
    prior(exponential(1), class = "sigma") 
  ),
  chains = 4, 
  iter = 6000,     
  warmup = 2000,   
  cores = 4,
  control = list(
    adapt_delta = 0.99,       
    max_treedepth = 15        
  )
)

summary(m_rapoport_with_size_L)
summary(m_rapoport_with_size_U)

