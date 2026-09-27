# Exploratory plots for range size, wingspan and phylogeny.

library(ggplot2)
library(phytools)
library(dplyr)
library(ggeffects)
library(tidyr)
library(viridis)
library(ggtree)
library(tidybayes)
library(broom.mixed)
library(ggnewscale)
library(ggtreeExtra)
library(sjPlot)
library(ggeffects)
library(ggplot2)
library(stringr)
library(car)

df_plot2_1 <- df_lat_sp %>%
  filter(is.finite(WS_L), WS_L > 0, is.finite(WS_U), WS_U > 0) %>%
  pivot_longer(
    cols = c(WS_L, WS_U),
    names_to = "WS_type",
    values_to = "WS_value"
  ) %>%
  mutate(WS_type = factor(WS_type, levels = c("WS_U", "WS_L"), labels = c("Upper (WS_U)", "Lower (WS_L)")))

fix_L <- fixef(m_lat_L)
fix_U <- fixef(m_lat_U)

pred_L <- predict_response(m_lat_L, terms = "abs_lat",back_transform = FALSE) %>%
  as_tibble() %>%
  mutate(WS_type = "Lower (WS_L)")

pred_U <- predict_response(m_lat_U, terms = "abs_lat",back_transform = FALSE) %>%
  as_tibble() %>%
  mutate(WS_type = "Upper (WS_U)")

df_preds2_1 <- bind_rows(pred_L, pred_U) %>%
  rename(abs_lat = x, log10_WS = predicted)

ggplot() +
  geom_ribbon(data = df_preds2_1, aes(x = abs_lat, ymin = conf.low, ymax = conf.high, fill = WS_type),
              alpha = 0.2) +
  geom_point(data = df_plot2_1, aes(x = abs_lat, y = log10(WS_value), color = WS_type, shape = Family),alpha=0.5,size=1.3
           ) +
  geom_line(data = df_preds2_1, aes(x = abs_lat, y = log10_WS, color = WS_type),
            size = 1.2) +
  scale_color_manual(values = c("Lower (WS_L)" = "#56B4E9", "Upper (WS_U)" = "#D55E00")) +
  scale_fill_manual(values = c("Lower (WS_L)" = "#56B4E9", "Upper (WS_U)" = "#D55E00")) +
  scale_shape_manual(values = c(16, 17, 15, 18, 25)) +
  labs(
    title = "Bergmann's Rule in migratry Butterflies",
    subtitle = "Points show raw data; Ribbons show 95% CI from LMM (Fixed Effects)",
    x = "Absolute Latitude",
    y = expression(log[10] * " Wing Span (mm)"),
    color = "Wing Span Type",
    fill = "Wing Span Type",
    shape = "Family"
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 12)
  )

pred_pattern <- predict_response(m_pattern, terms = "prop_mean", back_transform = FALSE) %>%
  as_tibble() %>%
  rename(prop_mean = x, log10_range = predicted)

pred_mechanism <- predict_response(m_mechanism, terms = "mean_bio4", back_transform = FALSE) %>%
  as_tibble() %>%
  rename(mean_bio4 = x, log10_range = predicted)

ggplot() +
  geom_ribbon(data = pred_pattern, aes(x = prop_mean, ymin = conf.low, ymax = conf.high),
              fill = "#D55E00", alpha = 0.2) +
  geom_point(data = final_df, aes(x = prop_mean, y = log10(range_km2), color = mean_bio4),
             alpha = 0.3, size = 1) +
  scale_color_viridis_c(option = "viridis") +

  geom_line(data = pred_pattern, aes(x = prop_mean, y = log10_range),
            color = "#D55E00", size = 1.2) +

  labs(
    title = "Figure 2.2: Rapoport Pattern",
    subtitle = "Relationship between range size and tropical occupancy",
    x = "Proportion of range within tropics",
    y = expression(log[10] * " Range Size (" * km^2 * ")"),
    color = "Temp Seasonality\n(Bio4)"
  ) +
  theme_classic()

ggplot() +
  geom_ribbon(data = pred_mechanism, aes(x = mean_bio4, ymin = conf.low, ymax = conf.high),
              fill = "#D55E00", alpha = 0.2) +
  geom_point(data = final_df, aes(x = mean_bio4, y = log10(range_km2), color = prop_mean),
             alpha = 0.3, size = 1) +
  scale_color_viridis_c(option = "viridis")+
  geom_line(data = pred_mechanism, aes(x = mean_bio4, y = log10_range),
            color = "#D55E00", size = 1.2) +
  labs(
    title = "Figure 2.3: Rapoport Mechanism",
    subtitle = "Climatic variability (Bio4) driving range expansion",
    x = "Temperature Seasonality (Bio4)",
    y = expression(log[10] * " Range Size (" * km^2 * ")")
  ) +
  theme_classic()

mech_compare <- bind_rows(
  broom.mixed::tidy(m_mechanism_full_z, effects = "fixed", conf.int = TRUE) %>%
    mutate(model = "Full data (no phylogeny)"),
  broom.mixed::tidy(m_mechanism_subset_z, effects = "fixed", conf.int = TRUE) %>%
    mutate(model = "Phylo subset (no phylogeny)"),
  broom.mixed::tidy(m_mechanism_phylo, conf.int = TRUE) %>%
    mutate(model = "Phylo subset (phylogenetic GLMM)")
) %>%
  filter(term == "mean_bio4_z") %>%
  mutate(model = factor(model, levels = c(
    "Full data (no phylogeny)",
    "Phylo subset (no phylogeny)",
    "Phylo subset (phylogenetic GLMM)"
  )))

ggplot(mech_compare, aes(x = estimate, y = model, color = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.15, linewidth = 0.9) +
  geom_point(size = 3) +
  scale_color_manual(values = c(
    "Full data (no phylogeny)"          = "#16A085",
    "Phylo subset (no phylogeny)"       = "#2E86C1",
    "Phylo subset (phylogenetic GLMM)"  = "#E67E22"
  )) +
  labs(
    title = "Figure 2c robustness check",
    subtitle = "Standardized effect of temperature seasonality (Bio4) on log10 Range Size,\nwith vs. without phylogenetic control",
    x = "Standardized coefficient (95% CI/CrI)",
    y = NULL
  ) +
  theme_classic() +
  theme(legend.position = "none")

df_synthesis <- final_df %>%
  inner_join(df_lat_sp %>% dplyr::select(species, WS_L, WS_U, abs_lat), by = "species") %>%
  filter(!is.na(WS_U), !is.na(range_km2))

df_sp_level <- df_synthesis %>%
  group_by(species) %>%
  summarise(
    mean_range = mean(range_km2),
    prop_mean=first(prop_mean),
    WS_U = first(WS_U),
    WS_L = first(WS_L),
    abs_lat = first(abs_lat),
    .groups = "drop"
  )

ggplot(df_sp_level, aes(x = log10(WS_U), y = log10(mean_range))) +
  geom_point(aes(color = prop_mean), size = 2.5, alpha = 0.7) +
  geom_smooth(method = "lm", color = "black", linetype = "dashed", se = TRUE, fill = "grey80") +
  scale_color_viridis_c(option = "plasma", name = "Tropical propotion") +
  scale_shape_manual(values = c(16, 17, 15, 18, 25))+
labs(
    title = "Figure 2.4: Synthesis of Size and Range Dynamics",
    subtitle = "Trade-off between Wing Span and Geographic Distribution",
    x = expression(log[10] * " Wing Span (Upper, mm)"),
    y = expression(log[10] * " Mean Range Size (" * km^2 * ")")
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

final_df_lat_clean<-final_df_lat
final_df_lat_clean$log_WS<-log10(final_df_lat$WS_L)
final_df_lat_clean$log_range<-log10(final_df_lat$range_km2)

m_clean <- lmer(log_range ~ abs_lat*log_WS + prop_within + season + (1 | species),
                     data = final_df_lat_clean)
m_clean_vif <- lmer(log_range ~ abs_lat+log_WS + prop_within + season + (1 | species),
                data = final_df_lat_clean)

predictions_clean <- ggpredict(m_clean, terms = c("log_WS", "abs_lat [0, 18.2, 40, 60]"))

performance::check_model(m_clean)
performance::check_model(m_clean_vif)

vif(m_clean_vif)

clean_plot <- plot(predictions_clean) +
  labs(
    title = "Effect of Wing Size on Range Size Across Latitudes",
    x = "Wing Size (log10)",
    y = "Predicted Range Size (log10 km²)",
    color = "Absolute Latitude"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 12),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

print(clean_plot)

clean_plot <- plot(predictions_clean) +
  labs(
    title = "Effect of Wing Size on Range Size Across Latitudes",
    x = "Wing Size (log10)",
    y = "Predicted Range Size (log10 km²)",
    color = "Absolute Latitude",
    fill = "Absolute Latitude"
  ) +
  scale_color_viridis_d(option = "viridis") +
  scale_fill_viridis_d(option = "viridis") +
  theme_classic()

print(clean_plot)

interact_compare <- bind_rows(
  broom.mixed::tidy(m_interact_full_z, effects = "fixed", conf.int = TRUE) %>%
    mutate(model = "Full data (no phylogeny)"),
  broom.mixed::tidy(m_interact_subset_z, effects = "fixed", conf.int = TRUE) %>%
    mutate(model = "Phylo subset (no phylogeny)"),
  broom.mixed::tidy(m_interact_phylo, conf.int = TRUE) %>%
    mutate(model = "Phylo subset (phylogenetic GLMM)")
) %>%
  filter(term == "abs_lat_z:log_WS_U_z") %>%
  mutate(model = factor(model, levels = c(
    "Full data (no phylogeny)",
    "Phylo subset (no phylogeny)",
    "Phylo subset (phylogenetic GLMM)"
  )))

p_interact_compare <- ggplot(interact_compare, aes(x = estimate, y = model, color = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.15, linewidth = 0.9) +
  geom_point(size = 3) +
  scale_color_manual(values = c(
    "Full data (no phylogeny)"          = "#16A085",
    "Phylo subset (no phylogeny)"       = "#2E86C1",
    "Phylo subset (phylogenetic GLMM)"  = "#E67E22"
  )) +
  labs(
    title = "Wing size x Latitude interaction: robustness check",
    subtitle = "Standardized abs_lat x log(WS_U) interaction on log10 Range Size,\nwith vs. without phylogenetic control",
    x = "Standardized interaction coefficient (95% CI/CrI)",
    y = NULL
  ) +
  theme_classic() +
  theme(legend.position = "none")

print(p_interact_compare)

ws_mean  <- mean(log10(df_model_WS_U_lat$WS_U))
ws_sd    <- sd(log10(df_model_WS_U_lat$WS_U))
lat_mean <- mean(df_model_WS_U_lat$abs_lat)
lat_sd   <- sd(df_model_WS_U_lat$abs_lat)
lat_vals <- c(0, 18.2, 40, 60)
lat_z    <- (lat_vals - lat_mean) / lat_sd

n_grid <- 50
pred_grid <- expand.grid(
  log_WS_U_z    = seq(min(df_model_WS_U_lat$log_WS_U_z), max(df_model_WS_U_lat$log_WS_U_z), length.out = n_grid),
  abs_lat_z     = lat_z,
  prop_within_z = 0,
  season        = factor("S1", levels = levels(df_model_WS_U_lat$season))
)
lat_labels <- paste0(lat_vals, "°")
pred_grid$lat_group <- rep(factor(lat_labels, levels = lat_labels), each = n_grid)

pred_fitted <- fitted(m_interact_phylo, newdata = pred_grid, re_formula = NA, summary = TRUE)
pred_grid <- pred_grid %>%
  mutate(
    fit      = pred_fitted[, "Estimate"],
    lwr      = pred_fitted[, "Q2.5"],
    upr      = pred_fitted[, "Q97.5"],
    WS_U_raw = 10^(log_WS_U_z * ws_sd + ws_mean)
  )

p_interact <- ggplot(pred_grid, aes(x = WS_U_raw, y = fit, color = lat_group, fill = lat_group)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.1) +
  scale_x_log10() +
  scale_color_viridis_d(option = "viridis") +
  scale_fill_viridis_d(option = "viridis") +
  labs(
    title = "Wing Size x Latitude Interaction on Range Size (phylogeny-controlled)",
    subtitle = "Predicted from m_interact_phylo; other covariates held at reference values",
    x = "Wing Span, Upper (mm, log scale)",
    y = expression("Predicted " * log[10] * " Range Size (" * km^2 * ")"),
    color = "Absolute Latitude", fill = "Absolute Latitude"
  ) +
  theme_classic()

print(p_interact)

dir.create("output/phylo_export", showWarnings = FALSE, recursive = TRUE)
write_csv(interact_compare, "output/phylo_export/figure2d_interaction_phylo_comparison.csv")
ggsave("output/phylo_export/figure2d_interaction_phylo_comparison.png", p_interact_compare, width = 8, height = 4, dpi = 200)
ggsave("output/phylo_export/figure2d_wingsize_latitude_interaction.png", p_interact, width = 7, height = 5, dpi = 200)

family_data <- df_lat_sp %>%
  mutate(species = str_replace_all(species, " ", "_")) %>%
  distinct(species, Family)

df_phylo_final <- df_sp_level %>%
  mutate(species = str_replace_all(species, " ", "_")) %>%
  left_join(family_data, by = "species") %>%
  mutate(log_Range = log10(mean_range)) %>%
  filter(species %in% tree_final$tip.label) %>%
  as.data.frame()

df_tree_side <- df_phylo_final %>%
  rename(label = species) %>%
  dplyr::select(label, Family) %>%
  as.data.frame()

df_fruit_side <- df_phylo_final %>%
  rename(label = species,
         Family_Bar = Family,
         Range_Val = log_Range) %>%
  dplyr::select(label, Family_Bar, Range_Val) %>%
  as.data.frame()

p_final <- ggtree(tree_final, layout = "fan", open.angle = 15, linewidth = 0.3) %<+% df_tree_side

p_final <- p_final + geom_tree(aes(color = Family), linewidth = 0.6)

p_final <- p_final +
  geom_fruit(
    data = df_fruit_side,
    geom = geom_bar,
    mapping = aes(y = label, x = Range_Val, fill = Family_Bar),
    stat = "identity",
    orientation = "y",
    pwidth = 0.3,
    offset = 0.1,
    axis.params = list(
      axis = "x",
      text.size = 2,
      title = "log10 Range",
      title.size = 3
    )
  )

p_final <- p_final +
  scale_color_brewer(palette = "Set1", name = "Butterfly Family") +
  scale_fill_brewer(palette = "Set1", name = "Butterfly Family") +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  ) +
  labs(title = "Phylogenetic Distribution of Geographic Range")

print(p_final)

model_fixef <- m_rapoport_with_size_U %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  mutate(.variable = str_remove(.variable, "b_")) %>%
  filter(.variable != "Intercept")

p_3_2 <- ggplot(model_fixef, aes(y = .variable, x = .value)) +
  stat_pointinterval(.width = c(.66, .95), color = "steelblue") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Figure 3.2: Standardized Coefficients for Range Size",
    subtitle = "Points show posterior mean; bars show 66% and 95% CIs",
    x = "Standardized Effect Size (Estimate)",
    y = "Predictors"
  ) +
  theme_classic()

p_3_2

hyp <- as_tibble(VarCorr(m_rapoport_with_size_U)$species_phylo$sd[,1]^2)
res <- as_tibble(VarCorr(m_rapoport_with_size_U)$residual$sd[,1]^2)

df_var <- data.frame(
  Source = c("Phylogeny", "Residuals"),
  Variance = c(mean(hyp$value), mean(res$value))
) %>%
  mutate(Percentage = Variance / sum(Variance) * 100)

p_3_3 <- ggplot(df_var, aes(x = "", y = Percentage, fill = Source)) +
  geom_bar(stat = "identity", width = 0.5) +
  coord_polar("y", start = 0) +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")),
            position = position_stack(vjust = 0.5), color = "white") +
  scale_fill_manual(values = c("Phylogeny" = "#2c7bb6", "Residuals" = "#d7191c")) +
  theme_void() +
  labs(title = "Figure 3.3: Variance Partitioning of Range Size",
       subtitle = "Phylogenetic Signal vs. Residual Variation")

p_3_3

df_plot_final <- df_phylo_final %>%
  rename(label = species) %>%
  mutate(Range_Val = log10(mean_range)) %>%
  as.data.frame()
species_to_keep <- intersect(tree_final$tip.label, df_plot_final$label)
tree_pruned <- keep.tip(tree_final, species_to_keep)

range_vec <- df_plot_final$Range_Val
names(range_vec) <- df_plot_final$label
range_vec <- range_vec[tree_pruned$tip.label]
anc_res <- fastAnc(tree_pruned, range_vec)

df_for_tree_branches <- data.frame(label = names(range_vec), Range_Color_Branch = as.numeric(range_vec))
df_for_tree_nodes <- data.frame(node = as.integer(names(anc_res)), Range_Color_Node = as.numeric(anc_res))

df_for_fruit_ring <- df_plot_final %>%
  filter(label %in% tree_pruned$tip.label) %>%
  dplyr::select(label, Family) %>%
  rename(Family_Identity_Ring = Family) %>%
  as.data.frame()

df_for_fruit_text <- df_for_fruit_ring %>%
  group_by(Family_Identity_Ring) %>%
  summarise(label = label[ceiling(n()/2)], .groups = "drop") %>%
  rename(Family_Name_Text = Family_Identity_Ring)

p_iter <- ggtree(tree_pruned, layout = "fan", open.angle = 15, linewidth = 0.5)

p_iter$data <- p_iter$data %>%
  left_join(df_for_tree_branches, by = "label") %>%
  left_join(df_for_tree_nodes, by = "node") %>%
  mutate(Final_Evolutionary_Value = coalesce(Range_Color_Branch, Range_Color_Node))
p_iter <- p_iter +
  aes(color = Final_Evolutionary_Value) +
  geom_tree(linewidth = 0.8) +
  scale_color_viridis_c(option = "viridis", name = "log10 Range Size")

p_iter <- p_iter +
  new_scale_fill() +
  geom_fruit(
    data = df_plot_final,
    geom = geom_tile,
    mapping = aes(y = label, fill = Range_Val),
    width = 10,
    offset = 0.05,
    linewidth = 0
  ) +
  scale_fill_viridis_c(option = "viridis", guide = "none")

p_iter <- p_iter +
  new_scale_fill() +
  geom_fruit(
    data = df_for_fruit_ring,
    geom = geom_tile,
    mapping = aes(y = label, fill = Family_Identity_Ring),
    width = 5,
    offset = 0.08,
    linewidth = 0
  ) +
  scale_fill_brewer(palette = "Set3",  guide = "none")

p_iter <- p_iter +
  geom_fruit(
    data = df_for_fruit_text,
    geom = geom_text,
    mapping = aes(
      y = label,
      label = Family_Name_Text,
      angle = ifelse(angle > 180, angle + 90, angle - 90)
    ),
    offset = 0.65,
    size = 4.5,
    fontface = "bold",
    hjust = 0.5,
    check_overlap = TRUE
  )

p_iter <- p_iter +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15)
  ) +
  labs(title = "Phylogenetic Signal in Butterfly Geographic Range")

print(p_iter)

library(patchwork)
p_iter <- p_iter +
  scale_x_continuous(expand = expansion(mult = c(0.2, 0.1))) +
  geom_rootedge(rootedge = 75)

p_hist <- ggplot(df_plot_final, aes(x = Range_Val)) +
  geom_histogram(aes(fill = after_stat(x)),
                 bins = 30,
                 color = "white",
                 show.legend = FALSE) +
  scale_fill_viridis_c(option = "plasma") +

  labs(x = "log(Distribution Range)", y = "Density") +

  theme_minimal(base_size = 10) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 8, face = "bold"),
    axis.text = element_text(size = 8, face = "bold")
  )

final_plot <- p_iter +
  inset_element(
    p_hist,
    left = 0.32,
    bottom = 0.35,
    right = 0.54,
    top = 0.58,
    align_to = 'full'
  )

print(final_plot)
