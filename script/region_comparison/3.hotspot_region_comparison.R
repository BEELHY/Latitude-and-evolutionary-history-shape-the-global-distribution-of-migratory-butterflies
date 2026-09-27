# Hotspot metrics by region.

library(dplyr)
library(tidyr)
library(ggplot2)

in_file <- "output/region_comparison/1.计数Region_Richness_Distribution.csv"

raw <- read.csv(in_file, check.names = FALSE, fileEncoding = "UTF-8")

rich_cols <- setdiff(names(raw), c("Region", "Total_Migratory_Cells"))
S_vals    <- as.integer(rich_cols)

counts <- as.matrix(raw[, rich_cols])
rownames(counts) <- raw$Region
storage.mode(counts) <- "double"

stopifnot(all(abs(rowSums(counts) - raw$Total_Migratory_Cells) < 1e-6))

global_counts <- colSums(counts)
global_total  <- sum(global_counts)
tail_share    <- rev(cumsum(rev(global_counts))) / global_total

S_star <- S_vals[which(tail_share <= 0.03)[1]]
cat(sprintf("Occupied cells worldwide: %s\n", format(global_total, big.mark = ",")))
cat(sprintf("Hotspot cutoff S* = %d  (top %.2f%% of occupied cells)\n\n",
            S_star, 100 * tail_share[which(S_vals == S_star)]))

quantile_from_counts <- function(cnt, q) {
  S_vals[which(cumsum(cnt) >= sum(cnt) * q)[1]]
}

hotspot_idx    <- S_vals >= S_star
global_hotspot <- sum(global_counts[hotspot_idx])

metrics <- lapply(seq_len(nrow(counts)), function(i) {
  cnt <- counts[i, ]
  n   <- sum(cnt)
  hs  <- sum(cnt[hotspot_idx])
  data.frame(
    Region            = rownames(counts)[i],
    Occupied_cells    = n,
    Pct_world_occupied= 100 * n / global_total,
    Mean_richness     = sum(S_vals * cnt) / n,
    Median_richness   = quantile_from_counts(cnt, 0.50),
    P95_richness      = quantile_from_counts(cnt, 0.95),
    Max_richness      = max(S_vals[cnt > 0]),
    Hotspot_cells     = hs,
    Hotspot_share     = 100 * hs / global_hotspot,
    stringsAsFactors  = FALSE
  )
}) %>% bind_rows() %>% arrange(desc(Hotspot_share))

metrics$Hotspot_over_rep <- metrics$Hotspot_share / metrics$Pct_world_occupied

MIN_CELLS <- 1000
metrics$Flag <- ifelse(metrics$Occupied_cells < MIN_CELLS, "small-sample*", "")

metrics_out <- metrics %>%
  mutate(across(c(Pct_world_occupied, Mean_richness,
                  Hotspot_share, Hotspot_over_rep), ~ round(.x, 2)))

write.csv(metrics_out, "output/region_comparison/Table_S_Region_Hotspot_Metrics.csv", row.names = FALSE)
print(metrics_out, row.names = FALSE)

shares <- metrics %>%
  select(Region, Pct_world_occupied, Hotspot_share) %>%
  pivot_longer(-Region, names_to = "metric", values_to = "pct") %>%
  mutate(
    Region = factor(Region, levels = rev(metrics$Region)),
    metric = factor(metric,
                    levels = c("Pct_world_occupied", "Hotspot_share"),
                    labels = c("% of world's occupied cells",
                               "% of world's hotspot cells"))
  )

p_share <- ggplot(shares, aes(pct, Region, fill = metric)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.68) +
  scale_fill_manual(values = c("% of world's occupied cells" = "grey72",
                               "% of world's hotspot cells"  = "#0072B2"),
                    name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.06)),
                     labels = function(x) paste0(x, "%")) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 10) +
  theme(
    legend.position    = "top",
    legend.key.size    = unit(0.8, "lines"),
    axis.line.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.2, colour = "grey92")
  )

ggsave("output/region_comparison/Figure_S_Region_Global_Shares.pdf", p_share,
       width = 6.4, height = 5.0, device = pdf)
ggsave("output/region_comparison/Figure_S_Region_Global_Shares.png", p_share,
       width = 6.4, height = 5.0, dpi = 400)

exceed <- lapply(seq_len(nrow(counts)), function(i) {
  cnt <- counts[i, ]
  data.frame(
    Region = rownames(counts)[i],
    S      = S_vals,
    P      = rev(cumsum(rev(cnt))) / sum(cnt),
    stringsAsFactors = FALSE
  )
}) %>% bind_rows() %>% filter(P > 0)

top_regions <- c(head(metrics$Region[metrics$Flag == ""], 6), "North American")
top_regions <- unique(top_regions)

exceed <- exceed %>%
  mutate(
    grp       = ifelse(Region %in% top_regions, Region, "Other regions"),
    highlight = Region %in% top_regions
  )

pal <- c("#E69F00", "#56B4E9", "#009E73", "#D55E00",
         "#0072B2", "#CC79A7", "#000000")[seq_along(top_regions)]
names(pal) <- top_regions
pal <- c(pal, "Other regions" = "grey78")

p <- ggplot() +
  geom_line(data = filter(exceed, !highlight),
            aes(S, P, group = Region, colour = grp),
            linewidth = 0.35, alpha = 0.8) +
  geom_line(data = filter(exceed, highlight),
            aes(S, P, group = Region, colour = grp),
            linewidth = 0.9) +
  geom_vline(xintercept = S_star, linetype = "dashed",
             colour = "grey35", linewidth = 0.4) +
  annotate("text", x = S_star + 0.8, y = 0.62,
           label = sprintf("hotspot threshold\nS* = %d (top 3%% worldwide)", S_star),
           hjust = 0, size = 2.9, colour = "grey25", lineheight = 0.95) +
  scale_y_log10(breaks = c(1, 0.1, 0.01, 1e-3, 1e-4, 1e-5),
                labels = c("100%", "10%", "1%", "0.1%", "0.01%", "0.001%")) +
  scale_x_continuous(breaks = c(1, seq(8, 48, 8)),
                     expand = expansion(mult = c(0.01, 0.02))) +
  scale_colour_manual(values = pal, breaks = c(top_regions, "Other regions"),
                      name = NULL) +
  labs(
    x = expression("Migratory species richness threshold, " * italic(S)),
    y = expression(paste("Occupied cells with richness ", phantom() >= phantom(),
                         italic(S), " (% of region's occupied cells)"))
  ) +
  theme_classic(base_size = 10) +
  theme(
    legend.position      = c(0.99, 0.99),
    legend.justification = c(1, 1),
    legend.key.height    = unit(0.85, "lines"),
    legend.background    = element_rect(fill = alpha("white", 0.85), colour = NA),
    axis.line            = element_line(linewidth = 0.35),
    panel.grid.major.y   = element_line(linewidth = 0.2, colour = "grey92")
  )

ggsave("output/region_comparison/Figure_S_Region_Exceedance.pdf", p, width = 6.8, height = 4.6, device = pdf)
ggsave("output/region_comparison/Figure_S_Region_Exceedance.png", p, width = 6.8, height = 4.6, dpi = 400)

cat("\nWrote: Table_S_Region_Hotspot_Metrics.csv,",
    "Figure_S_Region_Global_Shares.{pdf,png},",
    "Figure_S_Region_Exceedance.{pdf,png}\n")
