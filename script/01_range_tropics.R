# Range size and tropicality.

# Libraries
library(dplyr)
library(stringr)
library(tidyverse)
library(ggplot2)
library(lme4)
library(lmerTest)
library(MuMIn)
library(ape)
library(brms)

# Species-season metrics
metrics <- read_csv("updatedata/species_season_metrics.csv")
out <- metrics %>%
  dplyr::select(species, season, prop_tropics, range_km2)

# Plot
ggplot(out, aes(range_km2, prop_tropics)) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "lm")

# Within/between-species model
df_ss <- out %>%
  filter(is.finite(prop_tropics), is.finite(range_km2), range_km2 > 0) %>%
  mutate(season = factor(season))

df_wb <- df_ss %>%
  group_by(species) %>%
  mutate(
    prop_mean   = mean(prop_tropics, na.rm = TRUE),  # Between-species tropicality
    prop_within = prop_tropics - prop_mean  # Within-species seasonal deviation
  ) %>%
  ungroup()

m_wb <- lmer(log10(range_km2) ~ prop_mean + prop_within + season + (1 | species),
             data = df_wb)

summary(m_wb)

m_wb_p <- lmer(log10(range_km2) ~ prop_mean + prop_within + season + (1 | species),
               data = df_wb)
summary(m_wb_p)

# Effect size, 95% CI
b  <- fixef(m_wb_p)["prop_mean"]
se <- summary(m_wb_p)$coefficients["prop_mean","Std. Error"]

# Wald 95% CI
b_ci <- b + c(-1, 1) * 1.96 * se

delta <- 0.1

mult_ci <- 10^(b_ci * delta)
pct_ci  <- (1 - mult_ci) * 100

mult_ci
pct_ci

# Species-level plot
df_between <- df_wb %>%
  group_by(species) %>%
  summarise(
    prop_mean = first(prop_mean),
    range_km2 = mean(range_km2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(is.finite(prop_mean), is.finite(range_km2), range_km2 > 0) %>%
  mutate(log10_range = log10(range_km2))

# Fixed-effect prediction, 95% CI
b <- fixef(m_wb_p)
V <- vcov(m_wb_p)

grid <- data.frame(
  prop_mean   = seq(0, 1, length.out = 200),
  prop_within = 0,
  season      = factor("S1", levels = levels(df_wb$season))  # Reference season
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

# Coefficient plot
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

# Add temperature seasonality
final_df <- df_wb %>%
  left_join(metrics %>% dplyr::select(species, season, mean_bio4),
            by = c("species", "season"))
# Pattern model
m_pattern<- lmer(log10(range_km2) ~ prop_mean+prop_within + season + (1 | species), data = final_df)

# Mechanism model
m_mechanism <- lmer(log10(range_km2) ~ mean_bio4 +prop_within+ season + (1 | species), data = final_df)

# Combined model
m_combined <- lmer(log10(range_km2) ~ prop_mean+prop_within+season+ mean_bio4 + (1 | species), data = final_df)

summary(m_pattern)
summary(m_combined)
anova(m_pattern, m_mechanism)

# Season interaction check
m_pattern_season<- lmer(log10(range_km2) ~ (prop_within +prop_mean)*season + (1 | species), data = final_df)
summary(m_pattern_season)

m_mechanism_season<- lmer(log10(range_km2) ~ (mean_bio4 +prop_within)* season + (1 | species), data = final_df)
summary(m_mechanism_season)
# No seasonal interaction

# Phylogenetic signal
tree_final <- read.tree("updatedata/phylogeny_matched.tre")

df_phylo <- final_df %>%
  filter(species %in% tree_final$tip.label)
df_phylo <-df_phylo %>%
  filter(species %in% tree_final$tip.label) %>%
  mutate(species = factor(species))
df_phylo$species_phylo <- df_phylo$species

df_phylo_scaled <- df_phylo %>%
  mutate(
    mean_bio4_z = as.numeric(scale(mean_bio4)),
    prop_mean_z = as.numeric(scale(prop_mean)),
    prop_within_z = as.numeric(scale(prop_within))
  )

print(paste("final species", length(unique(df_phylo$species))))

# 247 species matched

# Phylogenetic covariance matrix
if(!is.ultrametric(tree_final)) tree_final <- phytools::force.ultrametric(tree_final)
A <- vcv.phylo(tree_final)

# Phylogenetic model
m_rapoport_optimized <- brm(
  log10(range_km2) ~ prop_mean+prop_within_z + season +
    (1 | gr(species_phylo, dist = "gaussian")),
  data = df_phylo_scaled,
  data2 = list(species_phylo = A),
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

summary(m_rapoport_optimized)
plot(m_rapoport_optimized)
bayes_R2(m_rapoport_optimized)

draws <- as.data.frame(m_rapoport_optimized)
draws <- draws %>%
  mutate(
    var_phylo = sd_species_phylo__Intercept^2,
    var_res = sigma^2,
    lambda = var_phylo / (var_phylo + var_res)
  )
lambda_summary <- draws %>%
  summarise(
    mean = mean(lambda),
    median = median(lambda),
    lower_95 = quantile(lambda, 0.025),
    upper_95 = quantile(lambda, 0.975)
  )

print("Phylogenetic signal (lambda/H2):")
print(lambda_summary)

ggplot(draws, aes(x = lambda)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = lambda_summary$mean, linetype = "dashed", color = "red") +
  labs(
    title = "Posterior Distribution of Phylogenetic Signal (Lambda)",
    x = "Lambda (H2)",
    y = "Density"
  ) +
  theme_minimal()

# Bio4 model with phylogeny

final_df_z <- final_df %>%
  mutate(
    mean_bio4_z   = as.numeric(scale(mean_bio4)),
    prop_within_z = as.numeric(scale(prop_within))
  )

m_mechanism_full_z <- lmer(log10(range_km2) ~ mean_bio4_z + prop_within_z + season + (1 | species),
                            data = final_df_z)

m_mechanism_subset_z <- lmer(log10(range_km2) ~ mean_bio4_z + prop_within_z + season + (1 | species),
                              data = df_phylo_scaled)

m_mechanism_phylo <- brm(
  log10(range_km2) ~ mean_bio4_z + prop_within_z + season +
    (1 | gr(species_phylo, dist = "gaussian")),
  data = df_phylo_scaled,
  data2 = list(species_phylo = A),
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

summary(m_mechanism_full_z)
summary(m_mechanism_subset_z)
summary(m_mechanism_phylo)

m_sensitive <- brm(
  log10(range_km2) ~  season +
    (1 | gr(species_phylo, dist = "gaussian")),
  data = df_phylo_scaled,
  data2 = list(species_phylo = A),
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

summary(m_sensitive)

# Retained vs excluded species
df_wb$in_bpmm <- df_wb$species %in% df_phylo$species

t.test(log10(range_km2) ~ in_bpmm, data = df_wb)

# Intercept-only model
m_intercept <- brm(
  formula = log10(range_km2) ~ 1 + (1 | gr(species_phylo, dist = "gaussian")),
  data = df_phylo_scaled,
  data2 = list(species_phylo = A),
  family = gaussian(),
  prior = c(
    prior(normal(0, 1), class = "Intercept"),
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

summary(m_intercept)
