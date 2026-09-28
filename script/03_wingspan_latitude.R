# Wingspan-latitude models.
# Run after 01-02, same session

# Libraries
library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(purrr)
library(ggplot2)
library(lme4)
library(lmerTest)

# Import data
trait_range <- read_csv("output/trait_range.csv")

# Latitude per species-season
lat_df <- read_csv("updatedata/species_season_metrics.csv") %>%
  dplyr::select(species, season, mean_lat, mean_abs_lat, lat_min, lat_max, lat_span)

trait_range2 <- trait_range %>%
  mutate(species_key = str_replace_all(species, " ", "_")) %>%
  left_join(
    lat_df %>% rename(species_key = species),
    by = c("species_key", "season")
  ) %>%
  dplyr::select(-species_key)

# Species-level latitude, wingspan
df_lat_sp <- trait_range2 %>%
  group_by(species, Family) %>%
  summarise(
    abs_lat = mean(mean_abs_lat, na.rm = TRUE),
    WS_L    = first(WS_L[!is.na(WS_L)]),
    WS_U    = first(WS_U[!is.na(WS_U)]),
    .groups = "drop"
  ) %>%
  filter(is.finite(abs_lat))

# Save species-level table
write_csv(df_lat_sp %>% mutate(species = str_replace_all(species, " ", "_")),
          "output/df_lat_sp.csv")

# Mixed models, family intercept
m_lat_L <- lmer(log10(WS_L) ~ abs_lat + (1 | Family),
                data = df_lat_sp %>% filter(is.finite(WS_L), WS_L > 0))

m_lat_U <- lmer(log10(WS_U) ~ abs_lat + (1 | Family),
                data = df_lat_sp %>% filter(is.finite(WS_U), WS_U > 0))

summary(m_lat_L)
summary(m_lat_U)

# Plot
p_lat_U <- ggplot(df_lat_sp %>% filter(is.finite(WS_U), WS_U > 0),
                  aes(abs_lat, log10(WS_U))) +
  geom_point(aes(shape = Family), alpha = 0.5, size = 1.6) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_classic() +
  labs(x = "Mean absolute latitude of occupied cells (°)",
       y = expression(log[10]*"(WS_U)"))

p_lat_U

# Percent change per degree
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

# Add to range models
df_lat_sp <- df_lat_sp %>%
  mutate(species = str_replace_all(species, " ", "_"))

df_combined <- df_phylo_scaled %>%
  left_join(df_lat_sp %>% dplyr::select(abs_lat, species, WS_L, WS_U), by = "species")

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
    (1 | gr(species_phylo, cov = A)) + (1 | species),
  data = df_model_WS_L,
  data2 = list(A = A),
  family = gaussian(),
  prior = c(
    prior(normal(0, 1), class = "b"),
    prior(exponential(1), class = "sd"),
    prior(exponential(1), class = "sigma")
  ),
  chains = 4,
  iter = 6000,
  warmup = 2000,
  cores = 4, seed = 1,
  control = list(
    adapt_delta = 0.99,
    max_treedepth = 15
  )
)

m_rapoport_with_size_U <- brm(
  log10(range_km2) ~ mean_bio4_z + prop_within_z + log_WS_U_z + season +
    (1 | gr(species_phylo, cov = A)) + (1 | species),
  data = df_model_WS_U,
  data2 = list(A = A),
  family = gaussian(),
  prior = c(
    prior(normal(0, 1), class = "b"),
    prior(exponential(1), class = "sd"),
    prior(exponential(1), class = "sigma")
  ),
  chains = 4,
  iter = 6000,
  warmup = 2000,
  cores = 4, seed = 1,
  control = list(
    adapt_delta = 0.99,
    max_treedepth = 15
  )
)

summary(m_rapoport_with_size_L)
summary(m_rapoport_with_size_U)

# Wingspan and range
final_df_lat <- final_df %>%
 left_join(df_lat_sp , by = "species")

m_pattern_lat<- lmer(log10(range_km2) ~ abs_lat+prop_within + season + (1 | species), data = final_df_lat)
summary(m_pattern_lat)

m_pattern<- lmer(log10(range_km2) ~ prop_mean+prop_within + season + (1 | species), data = final_df_lat)

anova(m_pattern,m_pattern_lat)

# m_pattern fits better

# Additive models, no interaction
m_add_L <- lmer(log10(range_km2) ~ abs_lat + log10(WS_L) + prop_within + season + (1 | species),
                data = final_df_lat %>% filter(is.finite(WS_L), WS_L > 0))
m_add_U <- lmer(log10(range_km2) ~ abs_lat + log10(WS_U) + prop_within + season + (1 | species),
                data = final_df_lat %>% filter(is.finite(WS_U), WS_U > 0))
summary(m_add_L)
summary(m_add_U)

# Range per wingspan ratio

final_df_lat <- final_df_lat %>%
  mutate(expansion_efficiency = log10(range_km2) - log10(WS_L))

m_interact<- lmer(log10(range_km2) ~ abs_lat*log10(WS_L) + prop_within + season + (1 | species),
                  data = final_df_lat)
summary(m_interact)
# Latitude-dependent wingspan effect

m_efficiency <- lmer(expansion_efficiency ~ abs_lat + prop_within + season + (1 | species),
                     data = final_df_lat)

summary(m_efficiency)

# Interaction with phylogeny

final_df_lat_z <- final_df_lat %>%
  filter(!is.na(WS_U), WS_U > 0, is.finite(abs_lat)) %>%
  mutate(
    log_WS_U_z   = as.numeric(scale(log10(WS_U))),
    abs_lat_z    = as.numeric(scale(abs_lat)),
    prop_within_z = as.numeric(scale(prop_within))
  )

m_interact_full_z <- lmer(log10(range_km2) ~ abs_lat_z * log_WS_U_z + prop_within_z + season + (1 | species),
                           data = final_df_lat_z)

df_model_WS_U_lat <- df_model_WS_U %>%
  filter(is.finite(abs_lat)) %>%
  mutate(abs_lat_z = as.numeric(scale(abs_lat)))

m_interact_subset_z <- lmer(log10(range_km2) ~ abs_lat_z * log_WS_U_z + prop_within_z + season + (1 | species),
                             data = df_model_WS_U_lat)

m_interact_phylo <- brm(
  log10(range_km2) ~ abs_lat_z * log_WS_U_z + prop_within_z + season +
    (1 | gr(species_phylo, cov = A)) + (1 | species),
  data = df_model_WS_U_lat,
  data2 = list(A = A),
  family = gaussian(),
  prior = c(
    prior(normal(0, 1), class = "b"),
    prior(student_t(3, 0, 1), class = "sd"),
    prior(student_t(3, 0, 1), class = "sigma")
  ),
  chains = 4,
  iter = 6000,
  warmup = 2000,
  cores = 4, seed = 1,
  control = list(
    adapt_delta = 0.99,
    max_treedepth = 15
  )
)

summary(m_interact_full_z)
summary(m_interact_subset_z)
summary(m_interact_phylo)

# Phylo interaction predictions (Fig. S3)
ws_mu <- mean(log10(df_model_WS_U$WS_U)); ws_sd <- sd(log10(df_model_WS_U$WS_U))
lat_mu <- mean(df_model_WS_U_lat$abs_lat); lat_sd <- sd(df_model_WS_U_lat$abs_lat)
ws_grid <- seq(min(df_model_WS_U_lat$log_WS_U_z), max(df_model_WS_U_lat$log_WS_U_z), length.out = 50)
nd <- expand.grid(log_WS_U_z = ws_grid, abs_lat = c(0, 20, 40, 60)) %>%
  mutate(abs_lat_z = (abs_lat - lat_mu) / lat_sd, prop_within_z = 0,
         season = "S1")  # Reference season
ep <- posterior_epred(m_interact_phylo, newdata = nd, re_formula = NA)
pred_phylo <- nd %>%
  mutate(WS_U = 10^(log_WS_U_z * ws_sd + ws_mu), fit = colMeans(ep),
         lo = apply(ep, 2, quantile, 0.025), hi = apply(ep, 2, quantile, 0.975))
dir.create("output/phylo_export", showWarnings = FALSE, recursive = TRUE)
write_csv(pred_phylo, "output/phylo_export/interaction_phylo_predictions.csv")

# Species-level phylogenetic checks
sp_prop <- df_wb %>% group_by(species) %>% summarise(prop_mean = first(prop_mean), .groups = "drop")
tree_sp <- drop.tip(tree_final, which(duplicated(tree_final$tip.label)))
wing_phylo <- list()
for (w in c("WS_U", "WS_L")) {
  d <- df_lat_sp %>% left_join(sp_prop, by = "species") %>%
    filter(is.finite(.data[[w]]), .data[[w]] > 0, species %in% tree_sp$tip.label) %>%
    mutate(lw = log10(.data[[w]]), species_phylo = species)
  A_sp <- vcv.phylo(keep.tip(tree_sp, d$species))
  pri <- c(prior(normal(0, 1), class = "b"), prior(exponential(1), class = "sd"), prior(exponential(1), class = "sigma"))
  m_lat <- brm(lw ~ abs_lat + (1 | gr(species_phylo, cov = A_sp)), data = d, data2 = list(A_sp = A_sp),
               prior = pri, chains = 4, iter = 6000, warmup = 2000, cores = 4, seed = 1, control = list(adapt_delta = 0.99))
  m_trop <- brm(prop_mean ~ lw + (1 | gr(species_phylo, cov = A_sp)), data = d, data2 = list(A_sp = A_sp),
                prior = pri, chains = 4, iter = 6000, warmup = 2000, cores = 4, seed = 1, control = list(adapt_delta = 0.99))
  wing_phylo[[w]] <- rbind(data.frame(model = "wingspan ~ latitude", term = "abs_lat", n = nrow(d), t(fixef(m_lat)["abs_lat", ])),
                           data.frame(model = "tropical ~ wingspan", term = "lw", n = nrow(d), t(fixef(m_trop)["lw", ])))
  wing_phylo[[w]]$wingspan <- w
}
write_csv(bind_rows(wing_phylo), "output/phylo_export/wingspan_phylo_models.csv")
