# Sample attrition tables.

suppressMessages({
  library(terra)
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

t_start <- Sys.time()
dir.create("output/checkpoints", showWarnings = FALSE, recursive = TRUE)
ckpt_path <- function(name) file.path("output/checkpoints", paste0(name, ".rds"))
save_ckpt <- function(obj, name) { saveRDS(obj, ckpt_path(name)); invisible(obj) }
has_ckpt  <- function(name) file.exists(ckpt_path(name))
load_ckpt <- function(name) readRDS(ckpt_path(name))

log <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), sprintf(...), "\n", sep = "")

# 426 species with maps
log("PHASE 1a: sp_426 from cached output/range_tropics.csv")
out <- read_csv("output/range_tropics.csv", show_col_types = FALSE)
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

path_ras <- "data/SuitabilityMaps_MigratorySpecies"
files <- list.files(path_ras, pattern = "^Binary_S[1-4].*", full.names = TRUE, recursive = TRUE)
meta <- tibble(file = files) %>%
  mutate(stem = tools::file_path_sans_ext(basename(file))) %>%
  mutate(
    season  = str_match(stem, "^Binary_(S[1-4])")[, 2],
    species = str_match(stem, "^Binary_(S[1-4])(.*)$")[, 3]
  ) %>%
  dplyr::select(file, species, season)

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

# Replicate script 03 filter
df_sp_trait_check <- trait_range_u %>%
  group_by(species, Family) %>%
  summarise(prop_mean = mean(prop_tropics, na.rm = TRUE), .groups = "drop")
n_split_species <- df_sp_trait_check %>% count(species) %>% filter(n > 1) %>% nrow()
log("species split across >1 (species,Family) group: %d (0 = coordinator's group-by hypothesis is NOT the cause)", n_split_species)

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
log("Candidate 'df_sp_trait prop_mean filter' (script/03_wingspan_traits.R lines 64-73, replicated exactly): n=%d (target 377)", length(cand_prop_mean_filter))
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
  "script/03_wingspan_traits.R's df_sp_trait construction: group_by(species, Family) %>% summarise(prop_mean = mean(prop_tropics, na.rm=TRUE), ...) %>% filter(is.finite(prop_mean)) -- drops 49 of 426 species whose prop_tropics is NA in every season row (an upstream data-quality issue in script/02_range_tropics.R's calc_metrics(), unrelated to wingspan availability). Reproduces 377 EXACTLY, but note this does NOT match the user's stated 'both WS_L and WS_U present' rule -- within this 377-species set only 316 have both WS present. No WS-presence-based candidate reproduces 377 under any tested variant (NA-based, >0-based, with/without Family, top-5-family-restricted)."
} else {
  sprintf("No exact match found; closest WS-based candidate '%s' adopted.", def_377_key)
}
log("ADOPTED sp_377 definition [%s]: n=%d (target 377)", def_377_key, length(sp_377))

n_fam_377 <- length(unique(trait_sp$Family[trait_sp$species %in% sp_377 & !is.na(trait_sp$Family)]))
log("distinct families among sp_377: %d (manuscript says 5)", n_fam_377)

# 247 phylogeny-matched species
log("PHASE 1c: sp_247 (phylogenetic matching, deterministic fix applied)")

tree_path <- "data/phylogenic/ntDegen359_fossils_smith_brown_strategyA.tre"
tree <- ape::read.nexus(tree_path)
tip_mapping <- tibble(original_label = tree$tip.label) %>%
  mutate(extracted_name = str_extract(original_label, "[A-Z][a-z]+_[a-z]+(?=(_|$))"))

final_df_genus <- df_wb %>%
  distinct(species) %>%
  mutate(genus = str_extract(species, "^[A-Z][a-z]+"))

tip_mapping_genus <- tip_mapping %>%
  mutate(genus = str_extract(extracted_name, "^[A-Z][a-z]+")) %>%
  filter(!is.na(genus))

match_genus_proxies <- function(unmatched_species, available_tips_for_proxy, deterministic = TRUE) {
  if (deterministic) {
    unmatched_species_indexed <- unmatched_species %>%
      arrange(genus, species) %>%
      group_by(genus) %>% mutate(spec_rank = row_number()) %>% ungroup()
    available_tips_indexed <- available_tips_for_proxy %>%
      arrange(genus, original_label) %>%
      group_by(genus) %>% mutate(tip_rank = row_number()) %>% ungroup()
  } else {
    # Old order-dependent matching
    unmatched_species_indexed <- unmatched_species %>%
      group_by(genus) %>% mutate(spec_rank = row_number()) %>% ungroup()
    available_tips_indexed <- available_tips_for_proxy %>%
      group_by(genus) %>% mutate(tip_rank = row_number()) %>% ungroup()
  }
  unmatched_species_indexed %>%
    inner_join(available_tips_indexed %>% dplyr::select(genus, original_label, tip_rank),
               by = c("genus" = "genus", "spec_rank" = "tip_rank")) %>%
    mutate(match_type = "genus_proxy") %>%
    dplyr::select(species, original_label, match_type)
}

run_tree_matching <- function(deterministic) {
  exact_matches <- final_df_genus %>%
    distinct(species, genus) %>%
    inner_join(tip_mapping_genus, by = c("species" = "extracted_name", "genus" = "genus")) %>%
    mutate(match_type = "exact")

  unmatched_species <- final_df_genus %>%
    distinct(species, genus) %>%
    filter(!species %in% exact_matches$species)

  used_exact_tips <- exact_matches$original_label
  available_tips_for_proxy <- tip_mapping_genus %>% filter(!original_label %in% used_exact_tips)

  genus_proxies <- match_genus_proxies(unmatched_species, available_tips_for_proxy, deterministic)

  all_matches <- bind_rows(
    exact_matches %>% dplyr::select(species, original_label, match_type),
    genus_proxies %>% dplyr::select(species, original_label, match_type)
  )
  all_matches_unique <- all_matches %>% group_by(original_label) %>% slice(1) %>% ungroup()
  list(exact_n = nrow(exact_matches), proxy_n = nrow(genus_proxies), all_matches_unique = all_matches_unique)
}

res_fixed  <- run_tree_matching(deterministic = TRUE)
res_buggy  <- run_tree_matching(deterministic = FALSE)  # Comparison only

log("exact_matches: %d (fixed) / %d (pre-fix, should be identical -- exact matching doesn't use row_number)",
    res_fixed$exact_n, res_buggy$exact_n)
log("genus_proxy matches: %d (fixed) / %d (pre-fix)", res_fixed$proxy_n, res_buggy$proxy_n)
log("total matched species: %d (fixed) / %d (pre-fix)",
    nrow(res_fixed$all_matches_unique), nrow(res_buggy$all_matches_unique))

changed_identity <- length(union(
  setdiff(res_fixed$all_matches_unique$species, res_buggy$all_matches_unique$species),
  setdiff(res_buggy$all_matches_unique$species, res_fixed$all_matches_unique$species)
))
log("species whose match identity changed between pre-fix and fixed logic: %d", changed_identity)

tree_final <- keep.tip(tree, res_fixed$all_matches_unique$original_label)
tree_final$tip.label <- res_fixed$all_matches_unique$species[match(tree_final$tip.label, res_fixed$all_matches_unique$original_label)]
log("raw tree_final$tip.label length: %d (target 247)", length(tree_final$tip.label))

# Count unique species
dup_tip_labels <- tree_final$tip.label[duplicated(tree_final$tip.label)]
if (length(dup_tip_labels) > 0) {
  log("duplicate tip label(s) causing the 248-vs-247 raw-vector discrepancy: %s", paste(unique(dup_tip_labels), collapse = ", "))
}

df_phylo <- df_wb %>% filter(species %in% tree_final$tip.label) %>% mutate(species = factor(species))
sp_247 <- sort(unique(as.character(df_phylo$species)))
log("sp_247 = unique(df_phylo$species) (post-fix): n=%d (manuscript target: 247)", length(sp_247))

save_ckpt(tree_final, "tree_final_fixed")
save_ckpt(list(res_fixed = res_fixed, res_buggy = res_buggy, changed_identity = changed_identity,
                dup_tip_labels = dup_tip_labels), "phase1d_matching")

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
  "Generated by script/07_attrition_tables.R.",
  "",
  "FUNNEL (simplified per user direction): 568 (not locally computable) -> 426 -> 377 -> 247.",
  "An earlier version of this script/analysis used a 6-number funnel (568/426/403/377/247/207) with",
  "an intermediate 'seasonal-shift' (403) and 'climate-baseline BPMM' (207) step. Both steps have been",
  "REMOVED per explicit user direction (the user knows their own pipeline; those two numbers are no",
  "longer part of the funnel). This note documents only the current, simplified 3-step funnel.",
  "",
  sprintf("sp_426 (n=%d, target 426): unique species with a Binary_S[1-4] suitability raster in data/SuitabilityMaps_MigratorySpecies (from output/range_tropics.csv).", length(sp_426)),
  "",
  sprintf("sp_377 (n=%d, target 377): %s", length(sp_377), def_377),
  sprintf("  Full candidate sizes tested (all applied directly to sp_426, no intermediate pool): %s.",
          paste(names(sizes_377), sizes_377, sep = "=", collapse = ", ")),
  sprintf("  Species split across >1 (species,Family) group when replicating script/03_wingspan_traits.R's df_sp_trait group_by(species,Family): %d -- the coordinator's hypothesis that this group-by silently double-counts a species was tested explicitly and is NOT the cause (0 splits found).", n_split_species),
  sprintf("  IMPORTANT CAVEAT: the adopted sp_377 does not satisfy 'both WS_L and WS_U present' -- within it, WS_L is present for %d species and WS_U for %d (of 377). If the user's WS-based rule is confirmed to be correct after all, the true sp_377 for that rule is not reproducible from local code (closest candidates: WS_L/both=%d, WS_U/either=%d).",
          sum(!is.na(df_sp_trait$WS_L)), sum(!is.na(df_sp_trait$WS_U)), length(cand_WS_both_pos), length(cand_WS_either_pos)),
  "",
  sprintf("sp_247 (n=%d, target 247): unique(df_phylo$species), i.e. distinct species names in tree_final$tip.label -- NOT the raw length of tree_final$tip.label (which is %d after the bugfix). Evidence: script/02_range_tropics.R line 337 has the comment \"#247\" directly under print(paste(\"final species\", length(unique(df_phylo$species)))). The raw vector is one longer because two different phylogeny tips (%s) both regex-extracted to the same species label, producing a duplicated tip label; %%in%% membership testing is unaffected by that duplicate, so sp_247 comes out at exactly 247 once counted correctly. This is a deliberate substitution of 'unique(df_phylo$species)' for the literal 'tree_final$tip.label', matching both the script's own annotation and the manuscript's target number.",
          length(sp_247), length(tree_final$tip.label), paste(unique(dup_tip_labels), collapse = ", ")),
  sprintf("  Fix applied (script/02_range_tropics.R lines ~285-300): row_number() genus-proxy-matching replaced with a deterministic tie-break (arrange by species name / tip label alphabetically within genus before ranking). Exact matches=%d, genus-proxy matches=%d.", res_fixed$exact_n, res_fixed$proxy_n),
  sprintf("  Pre-fix (row_number on incoming data-frame order) gave the same total (exact=%d, proxy=%d, total=%d), but %d species' specific match identity differed between the two runs -- confirms the bug changes WHICH species get matched, not HOW MANY (per-genus match count is bounded by min(unmatched, available tips) regardless of tie-break order). sp_247 is unaffected by the change to sp_377's definition in this revision (tree matching is built from df_wb/sp_426 directly, independent of sp_377), so it is unchanged from the prior (6-number) version of this analysis.", res_buggy$exact_n, res_buggy$proxy_n, nrow(res_buggy$all_matches_unique), changed_identity),
  "",
  sprintf("Nesting check: sp_247 subset of sp_377: %s | sp_377 subset of sp_426: %s",
          length(nest_247_377) == 0, length(nest_377_426) == 0),
  if (length(nest_247_377) > 0) sprintf("  NOTE: %d sp_247 species are NOT in sp_377 -- phylogenetic tree-matching draws from the full sp_426 pool independent of the df_sp_trait/prop_mean-filter step, so sp_247 is not guaranteed to nest inside sp_377. Table S1's 377->247 row uses the plain set difference regardless.", length(nest_247_377)) else NULL
)
writeLines(notes[!sapply(notes, is.null)], "log/subsets_methodology_notes.txt")
log("saved output/subsets.rds and log/subsets_methodology_notes.txt")

# Table S1: attrition tests
log("PHASE 2: building Table S1 (stepwise attrition)")

if (has_ckpt("lat_df")) {
  lat_df <- load_ckpt("lat_df")
  log("lat_df loaded from checkpoint: %d rows", nrow(lat_df))
} else {
  log("computing lat_df (area-weighted mean |latitude| per species-season, ~4 min)...")
  t_l0 <- Sys.time()
  lat_metrics_one <- function(f) {
    r <- rast(f)
    if (!is.lonlat(r)) r <- project(r, "EPSG:4326", method = "near")
    occ <- (r == 1)
    if (global(occ, "sum", na.rm = TRUE)[1, 1] == 0) return(tibble(mean_abs_lat = NA_real_))
    lat <- init(r, "y")
    a <- cellSize(r, unit = "km")
    area_total <- global(mask(a, occ, maskvalues = 0), "sum", na.rm = TRUE)[1, 1]
    mean_abs_lat <- global(mask(abs(lat) * a, occ, maskvalues = 0), "sum", na.rm = TRUE)[1, 1] / area_total
    tibble(mean_abs_lat = mean_abs_lat)
  }
  lat_df <- meta %>%
    mutate(metrics = purrr::map(file, purrr::possibly(lat_metrics_one, otherwise = tibble(mean_abs_lat = NA_real_)))) %>%
    unnest(metrics) %>%
    dplyr::select(species, season, mean_abs_lat)
  log("lat_df done in %.1f min", as.numeric(Sys.time() - t_l0, units = "mins"))
  save_ckpt(lat_df, "lat_df")
}

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
