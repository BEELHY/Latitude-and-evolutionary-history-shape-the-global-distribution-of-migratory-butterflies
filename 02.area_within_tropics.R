# Load libraries
library(terra)
library(dplyr)
library(stringr)
library(tidyverse)
library(ggplot2)
library(lme4)
library(lmerTest)
library(MuMIn) 
library(dismo)
library(geodata)
library(raster)

# Path
path <- "data/SuitabilityMaps_MigratorySpecies"

# List all files
files <- list.files(path,
                    pattern = "^Binary_S[1-4].*",
                    full.names = TRUE,
                    recursive = TRUE)

# Parse season & species
meta <- tibble(file = files) %>%
  mutate(stem = tools::file_path_sans_ext(basename(file))) %>%
  mutate(
    season  = str_match(stem, "^Binary_(S[1-4])")[,2],
    species = str_match(stem, "^Binary_(S[1-4])(.*)$")[,3]
  ) %>%
  select(file, species, season)

tropic_lat <- 23.4366

calc_metrics <- function(f) {
  r <- rast(f)
  
  # ensure lon/lat so latitude is meaningful
  if (!is.lonlat(r)) {
    r <- project(r, "EPSG:4326", method = "near")  # keeps 0/1
  }
  
  occ  <- (r == 1)
  lat  <- init(r, "y")
  trop <- (lat >= -tropic_lat) & (lat <= tropic_lat)
  
  a <- cellSize(r, unit = "km")  # per-cell area in km^2
  
  # total range size (km^2)
  range_km2 <- global(mask(a, occ, maskvalues = 0), "sum", na.rm = TRUE)[1,1]
  
  # tropical part (km^2)
  trop_km2  <- global(mask(a, occ & trop, maskvalues = 0), "sum", na.rm = TRUE)[1,1]
  
  prop_tropics <- if (is.na(range_km2) || range_km2 == 0) NA_real_ else trop_km2 / range_km2
  
  tibble(prop_tropics = prop_tropics,
         range_km2     = range_km2)
}

out <- meta %>%
  rowwise() %>%
  mutate(tmp = list(calc_metrics(file))) %>%
  tidyr::unnest(tmp) %>%
  ungroup() %>%
  select(species, season, prop_tropics, range_km2) %>%
  arrange(species, season)

# Export output
write_csv(out, "output/range_tropics.csv")

# Import data
out <- read_csv("output/range_tropics.csv")

# Plot
ggplot(out, aes(range_km2, prop_tropics)) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "lm")

# Linear model
# Considering seasonal distribution
# Within and between species
df_ss <- out %>% 
  filter(is.finite(prop_tropics), is.finite(range_km2), range_km2 > 0) %>% 
  mutate(season = factor(season))

df_wb <- df_ss %>%
  group_by(species) %>%
  mutate(
    prop_mean   = mean(prop_tropics, na.rm = TRUE),      # between-species tropicality
    prop_within = prop_tropics - prop_mean              # within-species seasonal deviation
  ) %>%
  ungroup()

m_wb <- lmer(log10(range_km2) ~ prop_mean + prop_within + season + (1 | species),
             data = df_wb)

summary(m_wb)

m_wb_p <- lmer(log10(range_km2) ~ prop_mean + prop_within + season + (1 | species),
               data = df_wb)
summary(m_wb_p)

# Effect size
# 95% CI
b  <- fixef(m_wb_p)["prop_mean"]
se <- summary(m_wb_p)$coefficients["prop_mean","Std. Error"]

# 95% Wald CI for beta
b_ci <- b + c(-1, 1) * 1.96 * se

delta <- 0.1

mult_ci <- 10^(b_ci * delta)
pct_ci  <- (1 - mult_ci) * 100

mult_ci
pct_ci

# Plot
# One point per species (between-species pattern)
df_between <- df_wb %>%
  group_by(species) %>%
  summarise(
    prop_mean = first(prop_mean),
    range_km2 = mean(range_km2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(is.finite(prop_mean), is.finite(range_km2), range_km2 > 0) %>%
  mutate(log10_range = log10(range_km2))

# Prediction line + 95% CI from the mixed model (fixed effects only)
b <- fixef(m_wb_p)
V <- vcov(m_wb_p)

grid <- data.frame(
  prop_mean   = seq(0, 1, length.out = 200),
  prop_within = 0,
  season      = factor("S1", levels = levels(df_wb$season))  # reference season for display
)

X <- model.matrix(~ prop_mean + prop_within + season, grid)
grid$fit <- as.numeric(X %*% b)
grid$se  <- sqrt(diag(X %*% V %*% t(X)))
grid$lwr <- grid$fit - 1.96 * grid$se
grid$upr <- grid$fit + 1.96 * grid$se

# Plot
ggplot(df_between, aes(x = prop_mean, y = log10_range)) +
  geom_point(alpha = 0.4, size = 1.6) +
  geom_ribbon(
    data = grid,
    aes(x = prop_mean, ymin = lwr, ymax = upr),
    inherit.aes = FALSE,
    alpha = 0.2
  ) +
  geom_line(
    data = grid,
    aes(x = prop_mean, y = fit),
    inherit.aes = FALSE,
    linewidth = 1
  ) +
  theme_classic() +
  labs(
    x = "Mean proportion of range in the tropics (|lat| ≤ 23.4366°)",
    y = expression(log[10]*"(range size, km"^2*")")
  )

######################################
# Coefficient plot
# extract betas + Wald CIs
terms <- c("prop_mean", "prop_within")
se_fix <- sqrt(diag(vcov(m_wb_p)))[terms]
beta   <- fixef(m_wb_p)[terms]

coef_df <- data.frame(
  term = c("Between-species tropicality (prop_mean)",
           "Within-species seasonal deviation (prop_within)"),
  beta = as.numeric(beta),
  lwr  = as.numeric(beta - 1.96 * se_fix),
  upr  = as.numeric(beta + 1.96 * se_fix)
)

ggplot(coef_df, aes(x = beta, y = term)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2) +
  geom_point(size = 2) +
  theme_classic() +
  labs(x = "Effect on log10(range_km2)", y = NULL)












# try to add temperature seasonality data
path <- "data/climate/wc2.1_10m_bio_4.tif"

bio4_ras <- rast(path)

#get data+calculate mean
extract_mechanism_metrics <- function(f, bio_layer) {
  r <- rast(f)
  if (!is.lonlat(r)) r <- project(r, "EPSG:4326", method = "near")
  
  bio_layer <- resample(bio_layer, r, method="bilinear") 
  
  occ <- (r == 1)
  if (global(occ, "sum", na.rm = TRUE)[1, 1] == 0) return(NA_real_)
  
  occ_bio <- mask(bio_layer, r, maskvalues = 0)
  
  mean_bio4 <- global(occ_bio, "mean", na.rm = TRUE)[1, 1]
  return(mean_bio4)
}

# calculation
mechanism_df <- ras_meta %>%
  mutate(mean_bio4 = map_dbl(file, ~extract_mechanism_metrics(.x, bio4_ras)))

#rbind
final_df <- df_wb %>%
  left_join(mechanism_df %>% dplyr::select(species, season, mean_bio4), 
            by = c("species", "season"))
#pattern model
m_pattern<- lmer(log10(range_km2) ~ prop_mean+prop_within + season + (1 | species), data = final_df)

# meca model
m_mechanism <- lmer(log10(range_km2) ~ mean_bio4 + season + (1 | species), data = final_df)

# both model
m_combined <- lmer(log10(range_km2) ~ prop_mean+prop_within+season+ mean_bio4 + (1 | species), data = final_df)

summary(m_mechanism)
anova(m_pattern, m_combined)





