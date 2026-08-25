# Export phylogeny + species-level data for Figure 2d (Wing Size vs. Range Size synthesis)
# so a colleague can fit GLM/OLS models with phylogenetic control (e.g. PGLS via caper/nlme)
# outside this pipeline.
#
# Run this after script/05.plots.R in the same R session — it needs tree_final and
# df_phylo_final, both created there.

library(ape)
library(dplyr)
library(readr)

dir.create("output/phylo_export", showWarnings = FALSE, recursive = TRUE)

species_to_export <- intersect(tree_final$tip.label, df_phylo_final$species)
tree_export <- keep.tip(tree_final, species_to_export)

df_export <- df_phylo_final %>%
  filter(species %in% species_to_export) %>%
  dplyr::select(species, Family, mean_range, log_Range, WS_U, WS_L, prop_mean, abs_lat)

stopifnot(setequal(tree_export$tip.label, df_export$species))

write.tree(tree_export, "output/phylo_export/figure2d_tree.nwk")
write_csv(df_export, "output/phylo_export/figure2d_species_data.csv")

cat(sprintf(
  "Exported %d species:\n  tree -> output/phylo_export/figure2d_tree.nwk\n  data -> output/phylo_export/figure2d_species_data.csv\n",
  length(species_to_export)
))
