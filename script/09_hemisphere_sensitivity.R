# Hemisphere sensitivity test.

setwd("/Users/hlii0385/Desktop/Phd_Haiyu_LI/Rapoport-s-rule-and-Bergmann-s-rule-of-Migratory-butterflies")

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(stringr)
  library(lme4)
  library(lmerTest)
})

cat("================================================================\n")
cat("  Hemispheric Sensitivity Analysis\n")
cat("  Reference: Fig 2b (Bergmann) & Fig 2c (Wingspan×Lat×Range)\n")
cat("================================================================\n\n")

df_ws <- tryCatch(
  read.csv("data/df_lat_sp.csv") %>%
    rename(family = Family) %>%
    mutate(species = str_replace_all(species, " ", "_")),
  error = function(e) { cat("ERROR loading df_lat_sp.csv:", conditionMessage(e), "\n"); NULL }
)

tr <- tryCatch(
  read.csv("output/trait_range.csv") %>%
    mutate(species = str_replace_all(species, " ", "_")),
  error = function(e) { cat("ERROR loading trait_range.csv:", conditionMessage(e), "\n"); NULL }
)

if (is.null(df_ws) || is.null(tr)) stop("Cannot proceed without base data files.")

df_range <- tr %>%
  group_by(species) %>%
  summarise(range_km2_mean = mean(range_km2, na.rm = TRUE), .groups = "drop")

df <- df_ws %>%
  inner_join(df_range, by = "species") %>%
  filter(
    is.finite(WS_L), WS_L > 0,
    is.finite(WS_U), WS_U > 0,
    is.finite(range_km2_mean), range_km2_mean > 0
  ) %>%
  mutate(
    wingspan_upper  = log10(WS_U),
    wingspan_lower  = log10(WS_L),
    range_size_log10 = log10(range_km2_mean),
    abs_latitude    = abs_lat
  )

cat(sprintf("N after merging wingspan + range: %d species\n\n", nrow(df)))

signed_lat_cache <- "output/north/signed_lat_per_species.csv"

if (file.exists(signed_lat_cache)) {
  cat(sprintf("Loading signed centroid latitude from %s ...\n", signed_lat_cache))
  signed_lat_df <- read.csv(signed_lat_cache) %>%
    mutate(species = str_replace_all(species, " ", "_"))
  cat(sprintf("  Loaded %d species\n", nrow(signed_lat_df)))

} else {
  cat("Cache not found — computing signed centroid latitude from rasters...\n")
  suppressPackageStartupMessages(library(terra))

  ras_path  <- "data/SuitabilityMaps_MigratorySpecies"
  ras_files <- list.files(ras_path, pattern = "^Binary_S[1-4].*\\.tif$",
                          full.names = TRUE, recursive = FALSE)

  ras_meta <- tibble(file = ras_files) %>%
    mutate(
      stem    = tools::file_path_sans_ext(basename(file)),
      season  = str_extract(stem, "S[1-4]"),
      species = str_replace(stem, "^Binary_S[1-4]", "")
    )

  one_per_sp <- ras_meta %>%
    arrange(species, season) %>%
    group_by(species) %>% slice(1) %>% ungroup() %>%
    filter(species %in% df$species)

  cat(sprintf("  Will process %d species\n", nrow(one_per_sp)))

  get_signed_lat <- function(f) {
    tryCatch({
      r <- rast(f)
      if (!is.lonlat(r)) r <- project(r, "EPSG:4326", method = "near")
      lat    <- init(r, "y")
      masked <- mask(lat, r, maskvalues = 0)
      as.numeric(global(masked, "mean", na.rm = TRUE)[1, 1])
    }, error = function(e) NA_real_)
  }

  signed_lats <- numeric(nrow(one_per_sp))
  for (i in seq_len(nrow(one_per_sp))) {
    if (i %% 50 == 0 || i == 1)
      cat(sprintf("    [%d / %d]\n", i, nrow(one_per_sp)))
    signed_lats[i] <- get_signed_lat(one_per_sp$file[i])
  }

  signed_lat_df <- one_per_sp %>%
    mutate(centroid_latitude = signed_lats) %>%
    dplyr::select(species, centroid_latitude)

  write.csv(signed_lat_df, signed_lat_cache, row.names = FALSE)
  cat(sprintf("  Saved to %s\n", signed_lat_cache))
}

df <- df %>%
  left_join(signed_lat_df, by = "species")

n_missing <- sum(is.na(df$centroid_latitude))
cat(sprintf("  Missing signed lat: %d species (dropped)\n", n_missing))
df <- df %>% filter(!is.na(centroid_latitude))
cat(sprintf("  Final analysis N = %d species\n\n", nrow(df)))

cat("================================================================\n")
cat("  STEP 1 — Sample Size Check\n")
cat("================================================================\n\n")

df <- df %>%
  mutate(
    hemisphere = ifelse(centroid_latitude >= 0, "North", "South"),
    lat_zone   = ifelse(abs_latitude < 18, "tropics (<18°)", "extratropics (≥18°)")
  )

cat("── Total N per hemisphere ──\n")
hemi_tab <- table(df$hemisphere)
print(hemi_tab)

cat("\n── N per hemisphere × latitude zone ──\n")
cross_tab <- table(df$hemisphere, df$lat_zone)
print(cross_tab)

cat("\n── Summary (range of abs_latitude) ──\n")
df %>%
  group_by(hemisphere, lat_zone) %>%
  summarise(
    n          = n(),
    lat_min    = round(min(abs_latitude), 1),
    lat_max    = round(max(abs_latitude), 1),
    .groups    = "drop"
  ) %>%
  as.data.frame() %>%
  print()

cat("\n================================================================\n")
cat("  STEP 2 — Decision Fork\n")
cat("================================================================\n\n")

n_south        <- sum(df$hemisphere == "South")
n_south_trop   <- sum(df$hemisphere == "South" & df$lat_zone == "tropics (<18°)")
n_south_extrat <- sum(df$hemisphere == "South" & df$lat_zone == "extratropics (≥18°)")

cat(sprintf("Southern Hemisphere total N = %d\n", n_south))
cat(sprintf("  tropics (<18°)    : %d\n", n_south_trop))
cat(sprintf("  extratropics (≥18°): %d\n", n_south_extrat))

use_3A <- (n_south >= 30) &&
          (n_south_trop  >= 10) &&
          (n_south_extrat >= 10)

if (use_3A) {
  cat("\nDecision: STEP 3A (formal three-way interaction)\n")
  cat(sprintf("  Reason: South N=%d with adequate spread across the 18° threshold\n", n_south))
} else {
  cat("\nDecision: STEP 3B (separate-hemisphere models)\n")
  cat(sprintf("  Reason: South N=%d", n_south))
  if (n_south_trop < 10 || n_south_extrat < 10)
    cat(sprintf("; imbalanced across 18° threshold (trop=%d, extrat=%d)",
                n_south_trop, n_south_extrat))
  cat("\n")
}

if (use_3A) {
  cat("\n================================================================\n")
  cat("  STEP 3A — Formal hemisphere × wingspan × latitude interaction\n")
  cat("================================================================\n\n")

  df$hemisphere <- factor(df$hemisphere, levels = c("North", "South"))

  for (ws_type in c("upper", "lower")) {
    ws_col <- if (ws_type == "upper") "wingspan_upper" else "wingspan_lower"
    cat(sprintf("\n── wingspan_%s ──\n", ws_type))

    df_mod <- df[is.finite(df[[ws_col]]), ]

    f_pooled <- as.formula(sprintf(
      "range_size_log10 ~ %s * abs_latitude + (1 | family)", ws_col))
    f_hemi   <- as.formula(sprintf(
      "range_size_log10 ~ %s * abs_latitude * hemisphere + (1 | family)", ws_col))

    m_pooled <- tryCatch(
      lmer(f_pooled, data = df_mod, REML = FALSE),
      error   = function(e) { cat("  lmer (pooled) failed:", conditionMessage(e), "\n"); NULL },
      warning = function(w) {
        cat("  lmer (pooled) warning:", conditionMessage(w), "\n")
        suppressWarnings(lmer(f_pooled, data = df_mod, REML = FALSE))
      }
    )

    m_hemi <- tryCatch(
      lmer(f_hemi, data = df_mod, REML = FALSE),
      error   = function(e) { cat("  lmer (hemisphere) failed:", conditionMessage(e), "\n"); NULL },
      warning = function(w) {
        cat("  lmer (hemisphere) warning:", conditionMessage(w), "\n")
        suppressWarnings(lmer(f_hemi, data = df_mod, REML = FALSE))
      }
    )

    if (!is.null(m_pooled) && !is.null(m_hemi)) {
      cat("\n  LRT: pooled vs hemisphere-interacted model\n")
      print(anova(m_pooled, m_hemi))

      cat("\n  Fixed-effect coefficients containing 'hemisphere':\n")
      s  <- summary(m_hemi)$coefficients
      hemi_rows <- grep("hemisphere", rownames(s), value = TRUE)
      if (length(hemi_rows) > 0) {
        print(round(s[hemi_rows, , drop = FALSE], 4))
      } else {
        cat("  (no hemisphere terms found)\n")
      }

      cat("\n  Singular fit (pooled / hemi):", isSingular(m_pooled), "/", isSingular(m_hemi), "\n")
    }
  }

} else {
  cat("\n================================================================\n")
  cat("  STEP 3B — Separate-hemisphere models\n")
  cat("================================================================\n\n")

  fit_hemi_model <- function(df_sub, hemi_label, ws_col) {
    n <- nrow(df_sub)
    cat(sprintf("\n  [%s | %s]  N = %d\n", hemi_label, ws_col, n))

    if (n < 10) {
      cat("  SKIP: fewer than 10 observations\n")
      return(invisible(NULL))
    }

    fml <- as.formula(sprintf(
      "range_size_log10 ~ %s * abs_latitude + (1 | family)", ws_col))

    m <- tryCatch(
      suppressWarnings(lmer(fml, data = df_sub, REML = TRUE)),
      error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL }
    )

    if (is.null(m)) return(invisible(NULL))

    sing <- isSingular(m)
    s    <- summary(m)$coefficients
    int_row <- grep(":", rownames(s))

    if (length(int_row) > 0) {
      est  <- s[int_row, "Estimate"]
      se   <- s[int_row, "Std. Error"]
      ci_l <- est - 1.96 * se
      ci_u <- est + 1.96 * se
      wide <- (ci_u - ci_l) > 1
      cat(sprintf("  Singular fit: %s\n", sing))
      cat(sprintf("  Interaction (%s × abs_lat): β = %.4f  95%%CI [%.4f, %.4f]%s\n",
                  ws_col, est, ci_l, ci_u,
                  ifelse(wide, "  ← CI uninformatively wide", "")))
    } else {
      cat("  Interaction term not found in summary.\n")
    }

    invisible(m)
  }

  df_north <- df %>% filter(hemisphere == "North")
  df_south <- df %>% filter(hemisphere == "South")

  for (ws_type in c("upper", "lower")) {
    ws_col <- if (ws_type == "upper") "wingspan_upper" else "wingspan_lower"
    cat(sprintf("\n════ wingspan_%s ════\n", ws_type))

    fit_hemi_model(df_north[is.finite(df_north[[ws_col]]), ], "North", ws_col)
    fit_hemi_model(df_south[is.finite(df_south[[ws_col]]), ], "South", ws_col)
  }

  cat("\n\n── Zero-crossing assessment ──\n")
  cat("The ~18° zero-crossing is observable in:\n")

  for (ws_type in c("upper", "lower")) {
    ws_col <- if (ws_type == "upper") "wingspan_upper" else "wingspan_lower"
    cat(sprintf("\n  wingspan_%s:\n", ws_type))

    for (hemi in c("North", "South")) {
      df_sub <- df %>%
        filter(hemisphere == hemi, is.finite(.data[[ws_col]]))
      n_below <- sum(df_sub$abs_latitude < 18)
      n_above <- sum(df_sub$abs_latitude >= 18)
      cat(sprintf("    %s: N=%d  (<18deg: n=%d, >=18deg: n=%d)\n",
                  hemi, nrow(df_sub), n_below, n_above))
    }
  }

  cat("\n  → Zero-crossing near 18° is ")
  if (n_south_trop < 5 || n_south_extrat < 5) {
    cat("NOT reliably estimable in the Southern Hemisphere\n")
    cat("    (insufficient coverage on both sides of 18°).\n")
    cat("    The pattern is effectively only estimable from Northern Hemisphere data.\n")
  } else {
    cat("potentially visible in both hemispheres, but Southern CI is wide.\n")
  }
}

cat("\n================================================================\n")
cat("  Analysis complete.\n")
cat("================================================================\n")
