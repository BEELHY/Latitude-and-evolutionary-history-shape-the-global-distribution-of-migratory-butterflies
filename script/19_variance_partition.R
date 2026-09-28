# Integrated model, variance partitioning.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ape)
  library(brms)
  library(ggplot2)
})

# Species-season data
metrics <- read_csv("updatedata/species_season_metrics.csv", show_col_types = FALSE)
df_lat_sp <- read_csv("output/df_lat_sp.csv", show_col_types = FALSE)
tree <- read.tree("updatedata/phylogeny_matched.tre")
tree <- drop.tip(tree, which(duplicated(tree$tip.label)))
if (!is.ultrametric(tree)) tree <- phytools::force.ultrametric(tree)

d <- metrics %>%
  filter(is.finite(prop_tropics), range_km2 > 0) %>%
  group_by(species) %>%
  mutate(prop_mean = mean(prop_tropics), prop_within = prop_tropics - prop_mean) %>%
  ungroup() %>%
  left_join(df_lat_sp %>% select(species, WS_U), by = "species") %>%
  filter(species %in% tree$tip.label, is.finite(WS_U), WS_U > 0,
         is.finite(mean_bio4), is.finite(mean_bio15), is.finite(mean_elev)) %>%
  mutate(across(c(mean_bio4, mean_bio15, mean_elev, prop_mean, prop_within), ~ as.numeric(scale(.x)), .names = "{.col}_z"),
         log_WS_U_z = as.numeric(scale(log10(WS_U))),
         season = factor(season), species_phylo = species)
A <- vcv.phylo(keep.tip(tree, unique(d$species)))
A <- A / max(diag(A))  # Correlation scale
cat(sprintf("Data: %d observations, %d species\n", nrow(d), length(unique(d$species))))

# Integrated phylogenetic model
m_int <- brm(
  log10(range_km2) ~ mean_bio4_z + mean_bio15_z + mean_elev_z + prop_mean_z + prop_within_z +
    log_WS_U_z + season + (1 | gr(species_phylo, cov = A)) + (1 | species),
  data = d, data2 = list(A = A), family = gaussian(),
  prior = c(prior(normal(0, 1), class = "b"),
            prior(exponential(1), class = "sd"),
            prior(exponential(1), class = "sigma")),
  chains = 4, iter = 6000, warmup = 2000, cores = 4, seed = 1,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)
print(summary(m_int))

# Variance components per draw
dr <- as_draws_df(m_int)
X <- model.matrix(~ mean_bio4_z + mean_bio15_z + mean_elev_z + prop_mean_z + prop_within_z +
                    log_WS_U_z + season, data = d)
env_terms <- c("mean_bio4_z", "mean_bio15_z", "mean_elev_z", "prop_mean_z", "prop_within_z")
sea_terms <- grep("^season", colnames(X), value = TRUE)
B <- as.matrix(dr[, paste0("b_", colnames(X)[-1])])
colnames(B) <- colnames(X)[-1]
f_env <- X[, env_terms] %*% t(B[, env_terms])
f_ws  <- X[, "log_WS_U_z", drop = FALSE] %*% t(B[, "log_WS_U_z", drop = FALSE])
f_sea <- X[, sea_terms] %*% t(B[, sea_terms])
v_env <- apply(f_env, 2, var); v_ws <- apply(f_ws, 2, var); v_sea <- apply(f_sea, 2, var)
v_fix <- apply(f_env + f_ws + f_sea, 2, var)
v_phy <- dr$sd_species_phylo__Intercept^2; v_sp <- dr$sd_species__Intercept^2; v_res <- dr$sigma^2
tot <- v_fix + v_phy + v_sp + v_res
parts <- data.frame(Environment = v_env / tot, `Wingspan` = v_ws / tot, Season = v_sea / tot,
                    `Shared (fixed effects)` = (v_fix - v_env - v_ws - v_sea) / tot,
                    Phylogeny = v_phy / tot, `Species (non-phylogenetic)` = v_sp / tot,
                    Residual = v_res / tot, check.names = FALSE)
vp <- parts %>%
  pivot_longer(everything(), names_to = "component", values_to = "prop") %>%
  group_by(component) %>%
  summarise(mean = mean(prop), lo = quantile(prop, 0.025), hi = quantile(prop, 0.975), .groups = "drop")
print(vp)

# Environment vs phylogeny
dd <- parts$Environment - parts$Phylogeny
env_vs_phy <- data.frame(diff_mean = mean(dd), lo = quantile(dd, 0.025), hi = quantile(dd, 0.975),
                         p_phy_ge_env = mean(parts$Phylogeny >= parts$Environment),
                         ratio_median = median(parts$Phylogeny / parts$Environment))
print(env_vs_phy)

dir.create("output/phylo_export", showWarnings = FALSE, recursive = TRUE)
write_csv(vp, "output/phylo_export/variance_partition.csv")
write_csv(env_vs_phy, "output/phylo_export/environment_vs_phylogeny.csv")
fe <- as.data.frame(fixef(m_int)); fe$term <- rownames(fe)
write_csv(fe, "output/phylo_export/integrated_model_fixef.csv")
cat(sprintf("Max Rhat: %.3f\n", max(brms::rhat(m_int), na.rm = TRUE)))

# Figure: variance components
lv <- c("Environment", "Wingspan", "Season", "Shared (fixed effects)", "Phylogeny", "Species (non-phylogenetic)", "Residual")
p <- ggplot(vp %>% mutate(component = factor(component, levels = rev(lv))), aes(mean, component)) +
  geom_col(aes(fill = component), width = 0.65, show.legend = FALSE) +
  scale_fill_manual(values = c(Environment = "#1B7837", Wingspan = "#762A83", Season = "grey70",
                               `Shared (fixed effects)` = "grey85", Phylogeny = "#E08214",
                               `Species (non-phylogenetic)` = "#FDB863", Residual = "grey50")) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.2, orientation = "y") +
  scale_x_continuous(labels = scales::percent) +
  labs(x = "Proportion of variance in log10 range size (95% CrI)", y = NULL) +
  theme_classic(base_size = 12)
dir.create("output/SI", showWarnings = FALSE)
ggsave("output/SI/FigS_variance_partition.png", p, width = 6.5, height = 3.8, dpi = 400, bg = "white")
