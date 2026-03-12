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
library(stringr)


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
  scale_color_manual(values = c("Lower (WS_L)" = "#377eb8", "Upper (WS_U)" = "#e41a1c")) +
  scale_fill_manual(values = c("Lower (WS_L)" = "#377eb8", "Upper (WS_U)" = "#e41a1c")) +
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


#Rapoport

pred_pattern <- predict_response(m_pattern, terms = "prop_mean", back_transform = FALSE) %>%
  as_tibble() %>%
  rename(prop_mean = x, log10_range = predicted)

pred_mechanism <- predict_response(m_mechanism, terms = "mean_bio4", back_transform = FALSE) %>%
  as_tibble() %>%
  rename(mean_bio4 = x, log10_range = predicted)

ggplot() +
  geom_ribbon(data = pred_pattern, aes(x = prop_mean, ymin = conf.low, ymax = conf.high), 
              fill = "#377eb8", alpha = 0.2) +
  geom_point(data = final_df, aes(x = prop_mean, y = log10(range_km2)), 
             color = "grey30", alpha = 0.2, size = 1) +
  
  geom_line(data = pred_pattern, aes(x = prop_mean, y = log10_range), 
            color = "#377eb8", size = 1.2) +
  labs(
    title = "Figure 2.2: Rapoport Pattern",
    subtitle = "Relationship between range size and tropical occupancy",
    x = "Proportion of range within tropics",
    y = expression(log[10] * " Range Size (" * km^2 * ")")
  ) +
  theme_classic()


ggplot() +
  geom_ribbon(data = pred_mechanism, aes(x = mean_bio4, ymin = conf.low, ymax = conf.high), 
              fill = "#e41a1c", alpha = 0.2) +
  geom_point(data = final_df, aes(x = mean_bio4, y = log10(range_km2), color = prop_mean), 
             alpha = 0.3, size = 1) +
  scale_color_viridis_c(option = "plasma")+
  geom_line(data = pred_mechanism, aes(x = mean_bio4, y = log10_range), 
            color = "#e41a1c", size = 1.2) +
  labs(
    title = "Figure 2.3: Rapoport Mechanism",
    subtitle = "Climatic variability (Bio4) driving range expansion",
    x = "Temperature Seasonality (Bio4)",
    y = expression(log[10] * " Range Size (" * km^2 * ")")
  ) +
  theme_classic()

#size vs Range

df_synthesis <- final_df %>%
  inner_join(df_lat_sp %>% dplyr::select(species, WS_L, WS_U, abs_lat), by = "species") %>%
  filter(!is.na(WS_U), !is.na(range_km2))

df_sp_level <- df_synthesis %>%
  group_by(species) %>%
  summarise(
    mean_range = mean(range_km2),
    prop_mean=first(prop_mean),
    WS_U = first(WS_U),
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


#phylo

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
    mapping = aes(y = label, x = Range_Val, fill = Family_Bar), # 使用新列名
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



#p3-1 new
df_plot_final <- df_phylo_final %>%
  rename(label = species) %>%
  mutate(Range_Val = log10(mean_range)) %>%
  as.data.frame()
species_to_keep <- intersect(tree_final$tip.label, df_plot_final$label)
tree_pruned <- keep.tip(tree_final, species_to_keep)

# 2. 祖先重建 (ASR)
range_vec <- df_plot_final$Range_Val
names(range_vec) <- df_plot_final$label
range_vec <- range_vec[tree_pruned$tip.label]
anc_res <- fastAnc(tree_pruned, range_vec)

# 3. 【数据隔离核心】：给每一层要用的变量取不同的名字
# A. 树枝颜色专用 (注入树内部)
df_for_tree_branches <- data.frame(label = names(range_vec), Range_Color_Branch = as.numeric(range_vec))
df_for_tree_nodes <- data.frame(node = as.integer(names(anc_res)), Range_Color_Node = as.numeric(anc_res))

# B. 外圈色环专用 (外部传入 geom_fruit)
df_for_fruit_ring <- df_plot_final %>%
  filter(label %in% tree_pruned$tip.label) %>%
  dplyr::select(label, Family) %>%
  rename(Family_Identity_Ring = Family) %>% # 唯一列名
  as.data.frame()

# C. 文字标签专用 (外部传入 geom_fruit)
df_for_fruit_text <- df_for_fruit_ring %>%
  group_by(Family_Identity_Ring) %>%
  summarise(label = label[ceiling(n()/2)], .groups = "drop") %>%
  rename(Family_Name_Text = Family_Identity_Ring) # 唯一列名


# A. 初始化树
p_iter <- ggtree(tree_pruned, layout = "fan", open.angle = 15, linewidth = 0.5)

# B. 手动注入树枝颜色数据 (只注入颜色需要的数值)
p_iter$data <- p_iter$data %>%
  left_join(df_for_tree_branches, by = "label") %>%
  left_join(df_for_tree_nodes, by = "node") %>%
  mutate(Final_Evolutionary_Value = coalesce(Range_Color_Branch, Range_Color_Node))

# C. 绘制彩色全树 (使用 Final_Evolutionary_Value)
p_iter <- p_iter + 
  aes(color = Final_Evolutionary_Value) + 
  geom_tree(linewidth = 0.8) +
  scale_color_viridis_c(option = "viridis", name = "log10 Range\n(Evolutionary)")

# D. 添加 Family 厚色环 (使用 Family_Identity_Ring)
p_iter <- p_iter +
  new_scale_fill() + 
  geom_fruit(
    data = df_for_fruit_ring,
    geom = geom_tile,
    mapping = aes(y = label, fill = Family_Identity_Ring),
    width = 1.2,      # 大宽度使其连成环
    offset = 0.1,
    linewidth = 0
  ) +
  scale_fill_brewer(palette = "Set1",  guide = "none") 

# E. 在圆弧上方标注 Family 名字 (使用 Family_Name_Text)
p_iter <- p_iter +
  geom_fruit(
    data = df_for_fruit_text,
    geom = geom_text,
    mapping = aes(y = label, label = Family_Name_Text),
    offset = 0.5,      # 放在色环外侧
    size = 4,
    fontface = "bold",
    hjust = 0.5,
    check_overlap = TRUE
  )

# F. 最终修饰
p_iter <- p_iter +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15)
  ) +
  labs(title = "Phylogenetic Signal in Butterfly Geographic Range")

# 打印
print(p_iter)

