# Figure 3 tree and inset.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(ggplot2)
  library(ape)
  library(phytools)
  library(ggtree)
  library(ggtreeExtra)
  library(ggnewscale)
  library(patchwork)
  library(viridis)
  library(png)
  library(grid)
  library(ggimage)
})

out_dir <- "output/Manuscript/reproducibility_code"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Load data
trait_range <- read_csv("output/trait_range.csv", show_col_types = FALSE)
df_lat_sp   <- read_csv("output/df_lat_sp.csv", show_col_types = FALSE)
tree_final  <- read.tree("updatedata/phylogeny_matched.tre")
dup_idx <- which(duplicated(tree_final$tip.label))
if (length(dup_idx) > 0) {
  cat(sprintf("Dropping %d duplicated tip label(s): %s\n", length(dup_idx),
              paste(tree_final$tip.label[dup_idx], collapse = ", ")))
  tree_final <- drop.tip(tree_final, dup_idx)
}

final_df <- trait_range %>%
  mutate(species = str_replace_all(species, " ", "_")) %>%
  group_by(species) %>%
  mutate(prop_mean = mean(prop_tropics, na.rm = TRUE)) %>%
  ungroup()

df_sp_level <- final_df %>%
  inner_join(df_lat_sp %>% dplyr::select(species, abs_lat), by = "species") %>%
  group_by(species) %>%
  summarise(mean_range = mean(range_km2), Family = first(Family), .groups = "drop") %>%
  # Fill missing families from genus
  mutate(Family = case_when(
    is.na(Family) & str_detect(species, "^Acraea_|^Cirrochroa_|^Heliconius_") ~ "Nymphalidae",
    is.na(Family) & str_detect(species, "^Pontia_") ~ "Pieridae",
    TRUE ~ Family
  ))

df_phylo_final <- df_sp_level %>%
  mutate(log_Range = log10(mean_range)) %>%
  dplyr::filter(species %in% tree_final$tip.label) %>%
  as.data.frame()

cat(sprintf("Figure 3: %d species matched to the %d-tip phylogeny\n",
            nrow(df_phylo_final), length(tree_final$tip.label)))

species_to_keep <- intersect(tree_final$tip.label, df_phylo_final$species)
tree_pruned <- keep.tip(tree_final, species_to_keep)

range_vec <- setNames(df_phylo_final$log_Range, df_phylo_final$species)[tree_pruned$tip.label]
anc_res <- fastAnc(tree_pruned, range_vec)

df_for_tree_branches <- data.frame(label = names(range_vec), Range_Color_Branch = as.numeric(range_vec))
df_for_tree_nodes <- data.frame(node = as.integer(names(anc_res)), Range_Color_Node = as.numeric(anc_res))

df_plot_final <- df_phylo_final %>% rename(label = species) %>% mutate(Range_Val = log_Range)

df_for_fruit_ring <- df_plot_final %>%
  dplyr::filter(label %in% tree_pruned$tip.label) %>%
  dplyr::select(label, Family) %>%
  rename(Family_Identity_Ring = Family) %>% as.data.frame()

# Fan tree
BASE_SIZE <- 12
DPI <- 600
range_lim <- range(df_plot_final$Range_Val)

p_iter <- ggtree(tree_pruned, layout = "fan", open.angle = 15, linewidth = 0.5)
p_iter$data <- p_iter$data %>%
  left_join(df_for_tree_branches, by = "label") %>%
  left_join(df_for_tree_nodes, by = "node") %>%
  mutate(Final_Evolutionary_Value = coalesce(Range_Color_Branch, Range_Color_Node))

p_iter <- p_iter +
  aes(color = Final_Evolutionary_Value) +
  geom_tree(linewidth = 0.8) +
  scale_color_viridis_c(option = "viridis", limits = range_lim,
                        name = expression(atop("Range size", "(" * log[10] * ", " * km^2 * ")")))

p_iter <- p_iter +
  new_scale_fill() +
  geom_fruit(data = df_plot_final, geom = geom_tile,
             mapping = aes(y = label, fill = Range_Val),
             width = 10, offset = 0.05, linewidth = 0) +
  scale_fill_viridis_c(option = "viridis", limits = range_lim, guide = "none")

p_iter <- p_iter +
  new_scale_fill() +
  geom_fruit(data = df_for_fruit_ring, geom = geom_tile,
             mapping = aes(y = label, fill = Family_Identity_Ring),
             width = 5, offset = 0.08, linewidth = 0) +
  scale_fill_brewer(palette = "Set3", guide = "none")

# Root edge opens centre
p_iter <- p_iter +
  scale_x_continuous(expand = expansion(mult = c(0.2, 0.1))) +
  geom_rootedge(rootedge = 75, colour = NA)  # invisible: only opens the centre

# No labels or photos here
p_iter <- p_iter +
  theme(
    legend.position = "right",
    legend.title = element_text(size = BASE_SIZE),
    legend.text = element_text(size = BASE_SIZE),
    legend.key.height = unit(0.55, "cm"),
    legend.key.width = unit(0.35, "cm"),
    plot.title = element_blank(),
    plot.margin = margin(4, 4, 4, 4)
  )

# Export legend separately
g_full <- ggplotGrob(p_iter)
legend_grob <- g_full$grobs[[grep("guide-box", g_full$layout$name)[1]]]
png(file.path(out_dir, "Figure3_legend.png"), width = 4, height = 4,
    units = "in", res = DPI, bg = "transparent")
grid::grid.draw(legend_grob)
dev.off()
p_iter <- p_iter + theme(legend.position = "none")

# Histogram, all species
p_hist <- ggplot(df_plot_final, aes(x = Range_Val)) +
  geom_histogram(aes(fill = after_stat(x)), bins = 30, color = "white", linewidth = 0.2, show.legend = FALSE) +
  scale_fill_viridis_c(option = "viridis", limits = range_lim, oob = scales::squish) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = expression("Range size (" * log[10] * ", " * km^2 * ")"), y = "Species") +
  theme_classic(base_size = BASE_SIZE - 2) +
  theme(panel.background = element_blank(), plot.background = element_blank(),
        axis.line = element_line(linewidth = 0.4, colour = "black"),
        axis.ticks = element_line(linewidth = 0.4, colour = "black"),
        axis.text = element_text(colour = "black"))

# Export tree and histogram
IMG_PX <- 2200
ggsave(file.path(out_dir, "Figure3_base.png"), p_iter,
       width = IMG_PX / 300, height = IMG_PX / 300, dpi = DPI, bg = "white")
ggsave(file.path(out_dir, "Figure3_hist.png"), p_hist,
       width = 1.95, height = 1.5,  # Native size keeps text at 10 pt
       dpi = DPI, bg = "transparent")
cat(sprintf("Saved: %s and Figure3_hist.png\n", file.path(out_dir, "Figure3_base.png")))
cat("Next: run script/18_figure3_compose.py to overlay family labels + butterfly photos\n")
