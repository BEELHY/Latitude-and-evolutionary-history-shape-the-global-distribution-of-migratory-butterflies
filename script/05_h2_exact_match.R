# H2 exact-match sensitivity.

suppressMessages({
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(readr)
  library(ape)
  library(brms)
})

set.seed(1)
dir.create("output/h2_sensitivity", showWarnings = FALSE, recursive = TRUE)
log <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), sprintf(...), "\n", sep = "")

# Base table with covariates
log("Loading secondary data...")
metrics <- read_csv("updatedata/species_season_metrics.csv", show_col_types = FALSE)
tree <- read.tree("updatedata/phylogeny_matched.tre")

df_env <- metrics %>%
  filter(is.finite(prop_tropics), is.finite(range_km2), range_km2 > 0) %>%
  group_by(species) %>%
  mutate(prop_mean = mean(prop_tropics, na.rm = TRUE),
         prop_within = prop_tropics - prop_mean) %>%
  ungroup() %>%
  filter(species %in% tree$tip.label)

log("Base table: %d obs, %d species; missing bio15=%d, missing elev=%d",
    nrow(df_env), length(unique(df_env$species)),
    sum(is.na(df_env$mean_bio15)), sum(is.na(df_env$mean_elev)))

# Keep species with wingspan
ws_tbl <- read_csv("output/df_lat_sp.csv", show_col_types = FALSE) %>%
  dplyr::select(species, WS_L)

df_ws <- df_env %>%
  left_join(ws_tbl, by = "species") %>%
  filter(!is.na(WS_L), WS_L > 0, !is.na(mean_bio15), !is.na(mean_elev)) %>%
  mutate(species = factor(species), season = factor(season))

log("After WS_L + bio15/elev completeness filter: %d obs, %d species (manuscript: 783 obs, 207 species)",
    nrow(df_ws), length(unique(df_ws$species)))

# Match type and tree
matching <- read_csv("updatedata/phylogeny_matching.csv", show_col_types = FALSE)
log("Match types: %s", paste(capture.output(print(table(matching$match_type))), collapse = " | "))

exact_species <- matching %>% filter(match_type == "exact") %>% pull(species)
log("Exact-match species available: %d", length(exact_species))

dup_idx <- which(duplicated(tree$tip.label))
if (length(dup_idx) > 0) {
  log("Dropping %d duplicated tip label(s): %s", length(dup_idx),
      paste(tree$tip.label[dup_idx], collapse = ", "))
  tree <- drop.tip(tree, dup_idx)
}
stopifnot(!any(duplicated(tree$tip.label)))

build_A <- function(species_subset) {
  t <- keep.tip(tree, intersect(tree$tip.label, species_subset))
  if (!is.ultrametric(t)) t <- phytools::force.ultrametric(t)
  v <- vcv.phylo(t)
  v / max(diag(v))  # Correlation scale
}

fit_bpmm <- function(data, A, label) {
  data <- data %>%
    mutate(
      mean_bio4_z   = as.numeric(scale(mean_bio4)),
      mean_bio15_z  = as.numeric(scale(mean_bio15)),
      mean_elev_z   = as.numeric(scale(mean_elev)),
      prop_within_z = as.numeric(scale(prop_within)),
      species_phylo = species
    ) %>%
    filter(species_phylo %in% rownames(A)) %>%
    mutate(species_phylo = factor(species_phylo), species = factor(species))

  log("[%s] final modelling data: %d obs, %d species", label, nrow(data), length(unique(data$species)))

  mod <- brm(
    log10(range_km2) ~ mean_bio4_z + mean_bio15_z + mean_elev_z + prop_within_z + season +
      (1 | gr(species_phylo, cov = A)) + (1 | species),
    data = data,
    data2 = list(A = A),
    family = gaussian(),
    prior = c(
      prior(normal(0, 1), class = "b"),
      prior(exponential(1), class = "sd"),
      prior(exponential(1), class = "sigma")
    ),
    chains = 4, iter = 6000, warmup = 2000, cores = 4,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    seed = 1
  )

  draws <- as_draws_df(mod)
  lambda <- draws$sd_species_phylo__Intercept^2 /
    (draws$sd_species_phylo__Intercept^2 + draws$sd_species__Intercept^2 + draws$sigma^2)

  list(
    label = label,
    n_obs = nrow(data),
    n_species = length(unique(data$species)),
    model = mod,
    h2_mean   = mean(lambda),
    h2_median = median(lambda),
    h2_lo95   = quantile(lambda, 0.025),
    h2_hi95   = quantile(lambda, 0.975),
    fixef = fixef(mod)
  )
}

# Model A: full pool
log("=== Fitting Model A: full reconstructed pool ===")
A_full <- build_A(unique(df_ws$species))
res_full <- fit_bpmm(df_ws, A_full, "A_full_reconstructed")

# Model B: exact matches
log("=== Fitting Model B: exact-match-only subset ===")
df_exact <- df_ws %>% filter(species %in% exact_species)
A_exact  <- build_A(unique(as.character(df_exact$species)))
res_exact <- fit_bpmm(df_exact, A_exact, "B_exact_match_only")

# Report
summary_tbl <- tibble(
  model     = c(res_full$label, res_exact$label),
  n_species = c(res_full$n_species, res_exact$n_species),
  n_obs     = c(res_full$n_obs, res_exact$n_obs),
  H2_mean   = c(res_full$h2_mean, res_exact$h2_mean),
  H2_lo95   = c(res_full$h2_lo95, res_exact$h2_lo95),
  H2_hi95   = c(res_full$h2_hi95, res_exact$h2_hi95)
)

write_csv(summary_tbl, "output/h2_sensitivity/H2_exact_match_sensitivity_summary.csv")
saveRDS(list(res_full = res_full, res_exact = res_exact),
        "output/h2_sensitivity/H2_exact_match_sensitivity_models.rds")

log("=== DONE ===")
print(summary_tbl)
cat("\nFixed effects, full reconstructed pool:\n"); print(res_full$fixef)
cat("\nFixed effects, exact-match-only subset:\n"); print(res_exact$fixef)
