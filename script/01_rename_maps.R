# Standardise species map names.

library(dplyr)
library(stringr)
library(purrr)
library(tibble)

tax_map <- tribble(
  ~species_prev, ~species_up,
  "Allancastria_cerisyi",      "Zerynthia_cerisyi",
  "Appias_drusilla",           "Glutophrissa_drusilla",
  "Colias_sareptensis",        "Colias_alfacariensis",
  "Eurema_lisa",               "Pyrisitia_lisa",
  "Eurema_nicippe",            "Abaeis_nicippe",
  "Euthalia_nais",             "Symphaedra_nais",
  "Junonia_lavinia",           "Junonia_genoveva",
  "Libytheana_carinenta",      "Prolibythea_carinenta",
  "Papilio_anchisiades",       "Heraclides_anchisiades",
  "Papilio_cresphontes",       "Heraclides_cresphontes",
  "Papilio_menatius",          "Pterourus_menatius",
  "Papilio_paeon",             "Heraclides_paeon",
  "Papilio_thoas",             "Heraclides_thoas",
  "Papilio_torquatus",         "Heraclides_torquatus",
  "Papilio_troilus",           "Pterourus_troilus",
  "Plebejus_acmon",            "Icaricia_acmon",
  "Polyommatus_semiargus",     "Cyaniris_semiargus",
  "Protographium_agesilaus",   "Eurytides_agesilaus",
  "Protographium_philolaus",   "Eurytides_philolaus",
  "Speyeria_callippe",         "Argynnis_callippe",
  "Hylephila_phylaeus",        "Hylephila_phyleus"
)

dir_in <- "data/SuitabilityMaps_MigratorySpecies"
files <- list.files(dir_in, pattern = "^Binary_S[1-4]", full.names = TRUE, recursive = TRUE)

plan <- tibble(path = files,
               dir  = dirname(path),
               file = basename(path)) %>%
  mutate(
    m = str_match(file, "^Binary_(S[1-4])(_?)([^.]+)(\\.[^.]+)$"),
    season = m[,2],
    sep    = m[,3],
    sp_prev = m[,4],
    ext    = m[,5]
  ) %>%
  dplyr::select(-m) %>%
  left_join(tax_map, by = c("sp_prev" = "species_prev")) %>%
  mutate(
    new_file = if_else(!is.na(species_up),
                       paste0("Binary_", season, sep, species_up, ext),
                       file),
    new_path = file.path(dir, new_file)
  )

to_change <- plan %>% filter(file != new_file)
to_change %>%  dplyr::select(file, new_file)

stopifnot(!any(duplicated(to_change$new_path)))

ok <- file.rename(from = to_change$path, to = to_change$new_path)
if (!all(ok)) warning("Some files could not be renamed. Check permissions / existing filenames.")
