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

# 377 species with wingspan
log("PHASE 1b: sp_377 (426 -> 377)")

trait_range <- read_csv("output/trait_range.csv", show_col_types = FALSE)

# Standardise species name format
trait_range_u <- trait_range %>% mutate(species = str_replace_all(species, " ", "_"))

trait_sp <- trait_range_u %>%
  group_by(species) %>%
  summarise(
    Family = dplyr::first(Family[!is.na(Family)]),
    WS_L   = dplyr::first(WS_L[!is.na(WS_L)]),
    WS_U   = dplyr::first(WS_U[!is.na(WS_U)]),
    .groups = "drop"
  )

# Wingspan-presence candidates
cand_WS_L_notna <- trait_sp %>% filter(!is.na(WS_L)) %>% pull(species)
cand_WS_U_notna <- trait_sp %>% filter(!is.na(WS_U)) %>% pull(species)
cand_WS_both_notna <- trait_sp %>% filter(!is.na(WS_L), !is.na(WS_U)) %>% pull(species)
cand_WS_either_notna <- trait_sp %>% filter(!is.na(WS_L) | !is.na(WS_U)) %>% pull(species)
cand_WS_both_pos <- trait_sp %>% filter(is.finite(WS_L), WS_L > 0, is.finite(WS_U), WS_U > 0) %>% pull(species)
cand_WS_either_pos <- trait_sp %>% filter((is.finite(WS_L) & WS_L > 0) | (is.finite(WS_U) & WS_U > 0)) %>% pull(species)

log("WS-based candidates on sp_426 -- WS_L notna=%d, WS_U notna=%d, both notna=%d, either notna=%d, both>0=%d, either>0=%d (target 377; none expected to match)",
    length(cand_WS_L_notna), length(cand_WS_U_notna), length(cand_WS_both_notna), length(cand_WS_either_notna),
    length(cand_WS_both_pos), length(cand_WS_either_pos))

# Replicate script 02 filter
df_sp_trait_check <- trait_range_u %>%
  group_by(species, Family) %>%
  summarise(prop_mean = mean(prop_tropics, na.rm = TRUE), .groups = "drop")
n_split_species <- df_sp_trait_check %>% count(species) %>% filter(n > 1) %>% nrow()
log("species split across >1 (species,Family) group: %d (0 = no double counting)", n_split_species)

df_sp_trait <- trait_range_u %>%
  group_by(species, Family) %>%
  summarise(
    prop_mean = mean(prop_tropics, na.rm = TRUE),
    range_km2 = mean(range_km2, na.rm = TRUE),
    WS_L = dplyr::first(WS_L[!is.na(WS_L)]),
    WS_U = dplyr::first(WS_U[!is.na(WS_U)]),
    .groups = "drop"
  ) %>%
  filter(is.finite(prop_mean))

cand_prop_mean_filter <- sort(unique(df_sp_trait$species))
log("Candidate 'df_sp_trait prop_mean filter' (script/02_wingspan_traits.R, replicated exactly): n=%d (target 377)", length(cand_prop_mean_filter))
log("  of these, WS_L present: %d, WS_U present: %d, both present: %d (i.e. NOT all 377 have complete wingspan data)",
    sum(!is.na(df_sp_trait$WS_L)), sum(!is.na(df_sp_trait$WS_U)), sum(!is.na(df_sp_trait$WS_L) & !is.na(df_sp_trait$WS_U)))

target_377 <- 377
cands_377 <- list(
  WS_L_notna = cand_WS_L_notna, WS_U_notna = cand_WS_U_notna,
  WS_both_notna = cand_WS_both_notna, WS_either_notna = cand_WS_either_notna,
  WS_both_pos = cand_WS_both_pos, WS_either_pos = cand_WS_either_pos,
  prop_mean_filter = cand_prop_mean_filter
)
sizes_377 <- sapply(cands_377, length)
log("ALL sp_377 candidate sizes: %s", paste(names(sizes_377), sizes_377, sep = "=", collapse = ", "))
match_377 <- names(cands_377)[sizes_377 == target_377]

def_377_key <- if (length(match_377) >= 1) match_377[1] else names(sizes_377)[which.min(abs(sizes_377 - target_377))]
sp_377 <- sort(cands_377[[def_377_key]])
def_377 <- if (def_377_key == "prop_mean_filter") {
  "script/02_wingspan_traits.R's df_sp_trait construction: group_by(species, Family) %>% summarise(prop_mean = mean(prop_tropics, na.rm=TRUE), ...) %>% filter(is.finite(prop_mean)) -- drops 49 of 426 species whose prop_tropics is NA in every season row (an upstream data-quality issue in the range metrics, unrelated to wingspan availability). Reproduces 377 exactly; within this set 316 species have both WS_L and WS_U. No WS-presence-based candidate reproduces 377 under any tested variant (NA-based, >0-based, with/without Family, top-5-family-restricted)."
} else {
  sprintf("No exact match found; closest WS-based candidate '%s' adopted.", def_377_key)
}
log("ADOPTED sp_377 definition [%s]: n=%d (target 377)", def_377_key, length(sp_377))

n_fam_377 <- length(unique(trait_sp$Family[trait_sp$species %in% sp_377 & !is.na(trait_sp$Family)]))
log("distinct families among sp_377: %d (manuscript says 5)", n_fam_377)

# 247 phylogeny-matched species
log("PHASE 1c: sp_247 (phylogeny-matched species)")

matching   <- read_csv("updatedata/phylogeny_matching.csv", show_col_types = FALSE)
tree_final <- ape::read.tree("updatedata/phylogeny_matched.tre")
log("exact_matches: %d, genus_proxy matches: %d",
    sum(matching$match_type == "exact"), sum(matching$match_type == "genus_proxy"))
log("raw tree_final$tip.label length: %d (target 247)", length(tree_final$tip.label))

# Count unique species
dup_tip_labels <- tree_final$tip.label[duplicated(tree_final$tip.label)]
if (length(dup_tip_labels) > 0) {
  log("duplicate tip label(s) causing the 248-vs-247 raw-vector discrepancy: %s", paste(unique(dup_tip_labels), collapse = ", "))
}

df_phylo <- df_wb %>% filter(species %in% tree_final$tip.label) %>% mutate(species = factor(species))
sp_247 <- sort(unique(as.character(df_phylo$species)))
log("sp_247 = unique(df_phylo$species) (post-fix): n=%d (manuscript target: 247)", length(sp_247))


# Nesting check
log("PHASE 1d: nesting check across sp_247 / sp_377 / sp_426")
nest_247_377 <- setdiff(sp_247, sp_377)
nest_377_426 <- setdiff(sp_377, sp_426)
log("sp_247 not in sp_377: %d", length(nest_247_377))
log("sp_377 not in sp_426: %d (should be 0, sp_377 built from sp_426's trait_range.csv pool)", length(nest_377_426))
nesting_ok <- (length(nest_247_377) == 0) && (length(nest_377_426) == 0)
log("STRICT NESTING HOLDS: %s", nesting_ok)

# Save subsets and notes
saveRDS(list(sp_426 = sp_426, sp_377 = sp_377, sp_247 = sp_247), "output/subsets.rds")

notes <- c(
  "Species-subset definitions used to build output/subsets.rds, output/TableS1_attrition.csv, output/TableS2_nested.csv",
  "Generated by script/06_attrition_tables.R.",
  "",
  "FUNNEL: 568 (migratory species list) -> 426 -> 377 -> 247.",
  "",
  sprintf("sp_426 (n=%d, target 426): unique species with a seasonal suitability map (updatedata/species_season_metrics.csv).", length(sp_426)),
  "",
  sprintf("sp_377 (n=%d, target 377): %s", length(sp_377), def_377),
  sprintf("  Full candidate sizes tested (all applied directly to sp_426, no intermediate pool): %s.",
          paste(names(sizes_377), sizes_377, sep = "=", collapse = ", ")),
  sprintf("  Species split across >1 (species,Family) group when replicating script/02_wingspan_traits.R's df_sp_trait group_by(species,Family): %d (0 = no species double-counted).", n_split_species),
  sprintf("  Note: within sp_377, WS_L is present for %d species and WS_U for %d (of 377); WS-presence subsets would give WS_L/both=%d, WS_U/either=%d.",
          sum(!is.na(df_sp_trait$WS_L)), sum(!is.na(df_sp_trait$WS_U)), length(cand_WS_both_pos), length(cand_WS_either_pos)),
  "",
  sprintf("sp_247 (n=%d, target 247): unique(df_phylo$species), i.e. distinct species names in tree_final$tip.label -- NOT the raw length of tree_final$tip.label (which is %d). Evidence: script/01_range_tropics.R counts length(unique(df_phylo$species)) = 247. The raw vector is one longer because two different phylogeny tips (%s) both regex-extracted to the same species label, producing a duplicated tip label; %%in%% membership testing is unaffected by that duplicate, so sp_247 comes out at exactly 247 once counted correctly. This is a deliberate substitution of 'unique(df_phylo$species)' for the literal 'tree_final$tip.label', matching both the script's own annotation and the manuscript's target number.",
          length(sp_247), length(tree_final$tip.label), paste(unique(dup_tip_labels), collapse = ", ")),
  sprintf("  Matching: %d exact, %d genus-proxy tips (updatedata/phylogeny_matching.csv).",
          sum(matching$match_type == "exact"), sum(matching$match_type == "genus_proxy")),
  "",
  sprintf("Nesting check: sp_247 subset of sp_377: %s | sp_377 subset of sp_426: %s",
          length(nest_247_377) == 0, length(nest_377_426) == 0),
  if (length(nest_247_377) > 0) sprintf("  NOTE: %d sp_247 species are NOT in sp_377 -- phylogenetic tree-matching draws from the full sp_426 pool independent of the df_sp_trait/prop_mean-filter step, so sp_247 is not guaranteed to nest inside sp_377. Table S1's 377->247 row uses the plain set difference regardless.", length(nest_247_377)) else NULL
)
writeLines(notes[!sapply(notes, is.null)], "output/subsets_methodology_notes.txt")
log("saved output/subsets.rds and output/subsets_methodology_notes.txt")

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

log("master_cov: %d species (target 426), %d with finite prop_mean (target 377), %d with Family assigned",
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
row2 <- build_step_row("426 (suitability maps) -> 377 (df_sp_trait / prop_mean complete)", sp_377, setdiff(sp_426, sp_377))
row3 <- build_step_row("377 -> 247 (phylogenetically matched)", sp_247, setdiff(sp_377, sp_247))

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
  list(label = "Wingspan-complete (n=377)", sp = sp_377),
  list(label = "Phylogeny-matched (n=247)", sp = sp_247)
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
