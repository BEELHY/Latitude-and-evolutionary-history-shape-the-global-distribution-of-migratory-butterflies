# Hotspot metrics by zoogeographic realm.

library(dplyr)
library(tidyr)
library(ggplot2)

region2realm <- c(
  "North American"   = "Nearctic",
  "Mexican"          = "Nearctic",
  "Panamanian"       = "Panamanian",
  "Amazonian"        = "Neotropical",
  "South American"   = "Neotropical",
  "African"          = "Afrotropical",
  "Guineo-Congolian" = "Afrotropical",
  "Saharo-Arabian"   = "Saharo-Arabian",
  "Madagascan"       = "Madagascan",
  "Eurasian"         = "Palearctic",
  "Arctico-Siberian" = "Palearctic",
  "Chinese"          = "Sino-Japanese",
  "Tibetan"          = "Sino-Japanese",
  "Oriental"         = "Oriental",
  "Indo-Malayan"     = "Oriental",
  "Australian"       = "Australian",
  "Novozelandic"     = "Australian",
  "Papua-Melanesian" = "Papua-Melanesian"
)
write.csv(data.frame(Region = names(region2realm), Realm = unname(region2realm)),
          "output/region_comparison/Table_S_Region_to_Realm_Holt2013.csv", row.names = FALSE)

raw <- read.csv("output/region_comparison/1.计数Region_Richness_Distribution.csv",
                check.names = FALSE, fileEncoding = "UTF-8")

rich_cols <- setdiff(names(raw), c("Region", "Total_Migratory_Cells"))
S_vals    <- as.integer(rich_cols)
counts    <- as.matrix(raw[, rich_cols]); storage.mode(counts) <- "double"
rownames(counts) <- raw$Region

stopifnot(all(abs(rowSums(counts) - raw$Total_Migratory_Cells) < 1e-6))
stopifnot(setequal(rownames(counts), names(region2realm)))

global_counts <- colSums(counts)
global_total  <- sum(global_counts)
tail_share    <- rev(cumsum(rev(global_counts))) / global_total
S_star        <- S_vals[which(tail_share <= 0.03)[1]]
hotspot_idx   <- S_vals >= S_star
global_hot    <- sum(global_counts[hotspot_idx])

cat(sprintf("S* = %d | %s occupied cells | %s hotspot cells\n\n",
            S_star, format(global_total, big.mark = ","),
            format(global_hot, big.mark = ",")))

summarise_units <- function(mapping) {
  realms <- unique(mapping)
  lapply(realms, function(rl) {
    members <- names(mapping)[mapping == rl]
    cnt <- colSums(counts[members, , drop = FALSE])
    n   <- sum(cnt); hs <- sum(cnt[hotspot_idx])
    data.frame(
      Realm              = rl,
      N_regions          = length(members),
      Regions            = paste(members, collapse = " + "),
      Occupied_cells     = n,
      Pct_world_occupied = 100 * n / global_total,
      Hotspot_cells      = hs,
      Hotspot_share      = 100 * hs / global_hot,
      Max_richness       = max(S_vals[cnt > 0]),
      stringsAsFactors   = FALSE
    )
  }) %>% bind_rows() %>% arrange(desc(Hotspot_share))
}

realms <- summarise_units(region2realm)
realms$Hotspot_over_rep <- realms$Hotspot_share / realms$Pct_world_occupied

realms_out <- realms %>%
  mutate(across(c(Pct_world_occupied, Hotspot_share, Hotspot_over_rep),
                ~ round(.x, 2)))

write.csv(realms_out, "output/region_comparison/Table_S_Realm_Hotspot_Metrics.csv", row.names = FALSE)
print(realms_out %>% select(-Regions), row.names = FALSE)

COL_OCC <- "#E69F00"
COL_HOT <- "#0072B2"

fmt_pct <- function(x) {
  ifelse(x == 0,   "0%",
  ifelse(x < 0.01, "<0.01%",
  ifelse(x < 1,    paste0(formatC(x, format = "f", digits = 2), "%"),
                   paste0(formatC(x, format = "f", digits = 1), "%"))))
}

land <- read.csv("output/region_comparison/Table_S_Realm_Land_Cells.csv", stringsAsFactors = FALSE)
stopifnot(setequal(land$Realm, realms$Realm))

occ <- realms %>%
  left_join(land, by = "Realm") %>%
  mutate(
    Global_share = 100 * Occupied_cells / global_total,
    Within_realm = 100 * Occupied_cells / Land_cells
  )
stopifnot(all(occ$Occupied_cells <= occ$Land_cells))

occ <- occ %>% arrange(desc(Global_share))

write.csv(
  occ %>% select(Realm, Regions, Land_cells, Occupied_cells,
                 Global_share, Within_realm) %>%
    mutate(across(c(Global_share, Within_realm), ~ round(.x, 2))),
  "output/region_comparison/Table_S_Realm_Occupancy.csv", row.names = FALSE)

shares <- occ %>%
  select(Realm, Within_realm, Global_share) %>%
  pivot_longer(-Realm, names_to = "metric", values_to = "pct") %>%
  mutate(
    Realm  = factor(Realm, levels = rev(occ$Realm)),
    metric = factor(metric,
                    levels = c("Within_realm", "Global_share"),
                    labels = c("% of that realm's land cells that are occupied",
                               "% of the world's occupied cells")),
    lab    = fmt_pct(pct)
  )

p <- ggplot(shares, aes(pct, Realm, fill = metric)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.70) +
  geom_text(aes(label = lab), position = position_dodge(width = 0.78),
            hjust = -0.15, size = 2.6, colour = "grey30") +
  scale_fill_manual(values = c("% of that realm's land cells that are occupied" = COL_OCC,
                               "% of the world's occupied cells"                 = COL_HOT),
                    name = NULL) +
  guides(fill = guide_legend(reverse = TRUE)) +
  scale_x_continuous(limits = c(0, 46), breaks = seq(0, 40, 10),
                     expand = expansion(mult = c(0, 0.01)),
                     labels = function(x) paste0(x, "%")) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 10) +
  theme(
    legend.position    = "top",
    legend.key.size    = unit(0.8, "lines"),
    axis.line.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.2, colour = "grey93")
  )

ggsave("output/region_comparison/Figure_S_Realm_Global_Shares.pdf", p, width = 6.6, height = 4.6, device = pdf)
ggsave("output/region_comparison/Figure_S_Realm_Global_Shares.png", p, width = 6.6, height = 4.6, dpi = 400)

cat("\nWrote: Table_S_Realm_Hotspot_Metrics.csv,",
    "Figure_S_Realm_Global_Shares.{pdf,png},",
    "Table_S_Realm_Occupancy.csv,",
    "Table_S_Region_to_Realm_Holt2013.csv\n")
