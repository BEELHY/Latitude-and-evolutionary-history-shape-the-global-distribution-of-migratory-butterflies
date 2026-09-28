# Join wingspan and range.

# Libraries
library(tidyverse)
library(dplyr)
library(stringr)
library(lme4)
library(lmerTest)
library(ggplot2)

# Import data
migr_trait <- read_csv("updatedata/wingspan_traits.csv")
range <- read_csv("updatedata/species_season_metrics.csv") %>%
  dplyr::select(species, season, prop_tropics, range_km2)

# Wrangle data
colnames(migr_trait)[1] <- "species"
migr_trait <- migr_trait %>%
  arrange(species, is.na(WS_L) & is.na(WS_U)) %>%
  distinct(species, .keep_all = TRUE)  # One row per species
range <- range %>%
  mutate(species = str_replace_all(species, "_", " "))

# Combine tables
trait_range <- dplyr::left_join(range, migr_trait, by = c("species"))
trait_range <- unique(trait_range)

# Export
dir.create("output", showWarnings = FALSE)
write_csv(trait_range, "output/trait_range.csv")

# Species-level data
df_sp_trait <- trait_range %>%
  group_by(species, Family) %>%
  summarise(
    prop_mean = mean(prop_tropics, na.rm = TRUE),
    range_km2 = mean(range_km2, na.rm = TRUE),  # Optional covariate
    WS_L = dplyr::first(WS_L[!is.na(WS_L)]),
    WS_U = dplyr::first(WS_U[!is.na(WS_U)]),
    .groups = "drop"
  ) %>%
  filter(is.finite(prop_mean))

# Model
df_L <- df_sp_trait %>% filter(is.finite(WS_L), WS_L > 0)
df_U <- df_sp_trait %>% filter(is.finite(WS_U), WS_U > 0)

m_L <- lmer(log10(WS_L) ~ prop_mean + (1 | Family), data = df_L)
m_U <- lmer(log10(WS_U) ~ prop_mean + (1 | Family), data = df_U)

summary(m_L)
summary(m_U)

# Sensitivity
m_L2 <- lmer(log10(WS_L) ~ prop_mean + log10(range_km2) + (1 | Family),
             data = df_L %>% filter(is.finite(range_km2), range_km2 > 0))
m_U2 <- lmer(log10(WS_U) ~ prop_mean + log10(range_km2) + (1 | Family),
             data = df_U %>% filter(is.finite(range_km2), range_km2 > 0))

summary(m_L2)
summary(m_U2)

# prop_mean not significant

# Trend
effect_per_0.1 <- function(mod, term = "prop_mean", delta = 0.1) {
  coefs <- summary(mod)$coefficients
  b  <- coefs[term, "Estimate"]
  se <- coefs[term, "Std. Error"]

  b_ci <- b + c(-1, 1) * 1.96 * se

  mult    <- 10^(b * delta)
  mult_ci <- 10^(b_ci * delta)

  pct    <- (mult - 1) * 100
  pct_ci <- (mult_ci - 1) * 100

  data.frame(
    multiplier = mult,
    pct_change = pct,
    pct_lwr = min(pct_ci),
    pct_upr = max(pct_ci)
  )
}

effect_per_0.1(m_L)  # Lower wingspan
effect_per_0.1(m_U)  # Upper wingspan

df_sp_trait <- df_sp_trait %>%
  mutate(WS_mid = (WS_L + WS_U) / 2)

# Plot
ggplot(df_sp_trait, aes(prop_mean, log10(WS_mid))) +
  geom_point(alpha = 0.5, size = 1.6) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_classic() +
  labs(
    x = "Mean proportion of range in the tropics",
    y = expression(log[10]*"(wingspan metric)")
  )
