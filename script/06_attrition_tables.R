# Sample attrition tables.

suppressMessages({
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(purrr)
  library(lme4)
  library(lmerTest)
  library(ape)
  library(effsize)
})

set.seed(1)  # Reproducible simulated P
t_start <- Sys.time()

log <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), sprintf(...), "\n", sep = "")

# 426 species with maps
log("PHASE 1a: sp_426 from species-season metrics")
metrics <- read_csv("updatedata/species_season_metrics.csv", show_col_types = FALSE)
out <- metrics %>% dplyr::select(species, season, prop_tropics, range_km2)
sp_426 <- sort(unique(out$species))
log("sp_426: %d species (manuscript target: 426)", length(sp_426))

df_ss <- out %>%
  filter(is.finite(prop_tropics), is.finite(range_km2), range_km2 > 0) %>%
  mutate(season = factor(season))

df_wb <- df_ss %>%
  group_by(species) %>%
  mutate(
    prop_mean   = mean(prop_tropics, na.rm = TRUE),
    prop_within = prop_tropics - prop_mean
  ) %>%
  ungroup()

# Species with wingspan data
log("PHASE 1b: species with wingspan data")
trait_range <- read_csv("output/trait_range.csv", show_col_types = FALSE)
trait_sp <- trait_range %>%
  mutate(species = str_replace_all(species, " ", "_")) %>%
  group_by(species) %>%
  summarise(
    Family = dplyr::first(Family[!is.na(Family)]),
    WS_L   = dplyr::first(WS_L[!is.na(WS_L)]),
    WS_U   = dplyr::first(WS_U[!is.na(WS_U)]),
    .groups = "drop"
  )
sp_ws <- trait_sp %>% filter((is.finite(WS_L) & WS_L > 0) | (is.finite(WS_U) & WS_U > 0)) %>% pull(species) %>% sort()
log("sp_ws: %d species (upper %d, lower %d)", length(sp_ws),
    sum(is.finite(trait_sp$WS_U) & trait_sp$WS_U > 0), sum(is.finite(trait_sp$WS_L) & trait_sp$WS_L > 0))

# Phylogeny-matched species
log("PHASE 1c: phylogeny-matched species")
matching   <- read_csv("updatedata/phylogeny_matching.csv", show_col_types = FALSE)
tree_final <- ape::read.tree("updatedata/phylogeny_matched.tre")
sp_247 <- sort(intersect(unique(tree_final$tip.label), sp_426))
sp_ws_phy <- sort(intersect(sp_ws, sp_247))
log("sp_247: %d; with wingspan: %d", length(sp_247), length(sp_ws_phy))

saveRDS(list(sp_426 = sp_426, sp_ws = sp_ws, sp_247 = sp_247, sp_ws_phy = sp_ws_phy), "output/subsets.rds")
writeLines(c(
  "Species subsets (script/06_attrition_tables.R)",
  sprintf("568 migratory species -> %d with seasonal suitability maps", length(sp_426)),
  sprintf("-> %d with wingspan data (upper or lower)", length(sp_ws)),
  sprintf("-> %d with wingspan data and matched to the phylogeny (%d exact, %d congeneric tips overall)",
          length(sp_ws_phy), sum(matching$match_type == "exact"), sum(matching$match_type == "genus_proxy"))
), "output/subsets_methodology_notes.txt")

# Table S1: attrition tests
log("PHASE 2: building Table S1 (stepwise attrition)")

lat_df <- metrics %>% dplyr::select(species, season, mean_abs_lat)
log("lat_df: %d rows", nrow(lat_df))

# Use unfiltered 426 species
master_cov <- out %>%
  group_by(species) %>%
  summarise(
    prop_mean = { pt <- prop_tropics; if (all(is.na(pt))) NA_real_ else mean(pt, na.rm = TRUE) },
    range_km2 = mean(range_km2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(lat_df %>% group_by(species) %>% summarise(abs_lat = mean(mean_abs_lat, na.rm = TRUE), .groups = "drop"), by = "species") %>%
  left_join(trait_sp %>% dplyr::select(species, Family), by = "species")

log("master_cov: %d species, %d with finite prop_mean, %d with Family assigned",
    nrow(master_cov), sum(is.finite(master_cov$prop_mean)), sum(!is.na(master_cov$Family)))

fmt_cell <- function(retained_vals, excluded_vals) {
  retained_vals <- retained_vals[is.finite(retained_vals)]
  excluded_vals <- excluded_vals[is.finite(excluded_vals)]
  if (length(retained_vals) < 2 || length(excluded_vals) < 2) {
    return(sprintf("%.3g vs %.3g (insufficient n for test)", median(retained_vals), median(excluded_vals)))
  }
  wt <- suppressWarnings(wilcox.test(retained_vals, excluded_vals))
  cd <- effsize::cliff.delta(retained_vals, excluded_vals)
  sprintf("%.3g vs %.3g (δ = %+.2f, P = %.3f)", median(retained_vals), median(excluded_vals), cd$estimate, wt$p.value)
}

family_chisq_p <- function(retained_sp, excluded_sp) {
  fam <- master_cov %>% filter(species %in% c(retained_sp, excluded_sp), !is.na(Family)) %>%
    mutate(group = ifelse(species %in% retained_sp, "Retained", "Excluded"))
  tab <- table(fam$group, fam$Family)
  if (any(dim(tab) < 2) || nrow(fam) < 2) return(NA_character_)
  ct <- suppressWarnings(chisq.test(tab))
  if (any(ct$expected < 5)) {
    ct_sim <- suppressWarnings(chisq.test(tab, simulate.p.value = TRUE, B = 5000))
    return(sprintf("%.3f (simulated)", ct_sim$p.value))
  }
  sprintf("%.3f", ct$p.value)
}

build_step_row <- function(step_label, retained_sp, excluded_sp) {
  if (length(excluded_sp) == 0) {
    return(tibble(
      `Filtering step` = step_label, Retained = length(retained_sp), Excluded = 0,
      abs_lat = "no exclusions", prop_mean = "no exclusions", `log10(range_km2)` = "no exclusions",
      `Family composition P` = NA_character_
    ))
  }
  ret_cov <- master_cov %>% filter(species %in% retained_sp)
  exc_cov <- master_cov %>% filter(species %in% excluded_sp)
  tibble(
    `Filtering step` = step_label,
    Retained = length(retained_sp),
    Excluded = length(excluded_sp),
    abs_lat = fmt_cell(ret_cov$abs_lat, exc_cov$abs_lat),
    prop_mean = fmt_cell(ret_cov$prop_mean, exc_cov$prop_mean),
    `log10(range_km2)` = fmt_cell(log10(ret_cov$range_km2), log10(exc_cov$range_km2)),
    `Family composition P` = family_chisq_p(retained_sp, excluded_sp)
  )
}

row1 <- tibble(`Filtering step` = "568 (migratory species list) -> 426 (with suitability maps)",
               Retained = 426L, Excluded = 142L,
               abs_lat = "not assessable", prop_mean = "not assessable", `log10(range_km2)` = "not assessable",
               `Family composition P` = NA_character_)
row2 <- build_step_row("426 (suitability maps) -> wingspan data", sp_ws, setdiff(sp_426, sp_ws))
row3 <- build_step_row("wingspan data -> phylogenetically matched", sp_ws_phy, setdiff(sp_ws, sp_ws_phy))

table_s1 <- bind_rows(row1, row2, row3)
write_csv(table_s1, "output/TableS1_attrition.csv")
log("saved output/TableS1_attrition.csv (3 rows)")

# Table S2: nested re-estimation
log("PHASE 3: building Table S2 (nested-subset re-estimation)")

df_lat_sp <- master_cov %>% dplyr::select(species, abs_lat) %>%
  left_join(trait_sp %>% dplyr::select(species, WS_U), by = "species")

final_df_lat_full <- df_wb %>% left_join(df_lat_sp, by = "species")

fit_reversal_model <- function(subset_sp, data_obs) {
  d <- data_obs %>%
    filter(species %in% subset_sp, is.finite(WS_U), WS_U > 0, is.finite(abs_lat), is.finite(range_km2)) %>%
    mutate(species = factor(species), log_WS_U = log10(WS_U))
  m <- lmerTest::lmer(log10(range_km2) ~ abs_lat * log_WS_U + prop_within + season + (1 | species), data = d)
  b <- lme4::fixef(m)
  V <- as.matrix(vcov(m))
  b_ws <- b["log_WS_U"]; se_ws <- sqrt(V["log_WS_U", "log_WS_U"])
  b_int <- b["abs_lat:log_WS_U"]; se_int <- sqrt(V["abs_lat:log_WS_U", "abs_lat:log_WS_U"])
  cov_ws_int <- V["log_WS_U", "abs_lat:log_WS_U"]

  # Reversal latitude, delta method
  r <- -b_ws / b_int
  var_r <- (1 / b_int)^2 * se_ws^2 + (b_ws / b_int^2)^2 * se_int^2 -
    2 * (1 / b_int) * (b_ws / b_int^2) * cov_ws_int
  se_r <- sqrt(var_r)

  list(
    n_species = length(unique(d$species)), n_obs = nrow(d),
    ws_est = b_ws, ws_lwr = b_ws - 1.96 * se_ws, ws_upr = b_ws + 1.96 * se_ws,
    int_est = b_int, int_lwr = b_int - 1.96 * se_int, int_upr = b_int + 1.96 * se_int,
    rev_est = r, rev_lwr = r - 1.96 * se_r, rev_upr = r + 1.96 * se_r
  )
}

fmt_est <- function(est, lwr, upr) sprintf("%.3f [%.3f, %.3f]", est, lwr, upr)
fmt_rev <- function(est, lwr, upr) sprintf("%.1f° [%.1f°, %.1f°]", est, lwr, upr)

subset_specs <- list(
  list(label = "Wingspan data", sp = sp_ws),
  list(label = "Wingspan data, phylogeny-matched", sp = sp_ws_phy)
)

s2_rows <- purrr::map(subset_specs, function(spec) {
  res <- fit_reversal_model(spec$sp, final_df_lat_full)
  tibble(
    `Species subset` = spec$label,
    `n species (n obs)` = sprintf("%d (%d)", res$n_species, res$n_obs),
    `Wingspan at equator` = fmt_est(res$ws_est, res$ws_lwr, res$ws_upr),
    `Wingspan x |latitude|` = fmt_est(res$int_est, res$int_lwr, res$int_upr),
    `Reversal latitude (deg)` = fmt_rev(res$rev_est, res$rev_lwr, res$rev_upr)
  )
})
table_s2 <- bind_rows(s2_rows)
write_csv(table_s2, "output/TableS2_nested.csv")
log("saved output/TableS2_nested.csv (2 rows)")

log("ALL PHASES DONE at +%.1f min total", as.numeric(Sys.time() - t_start, units = "mins"))

cat("\n\n=== TABLE S1 (markdown) ===\n")
print(knitr::kable(table_s1, format = "markdown"))
cat("\n\n=== TABLE S2 (markdown) ===\n")
print(knitr::kable(table_s2, format = "markdown"))
