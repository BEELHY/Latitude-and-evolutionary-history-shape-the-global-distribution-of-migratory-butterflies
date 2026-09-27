# Figure for phylogenetic sensitivity GAMs.

suppressMessages({
  library(mgcv)
  library(gratia)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(ragg)
  library(spdep)
  library(ggrepel)
})

PUB_FONT <- "Helvetica"

bundle <- readRDS("output/phylo_export/figure4d_phylo_sensitivity_models.rds")
model_df <- bundle$data

coords <- as.matrix(model_df[, c("lon1", "lat1")])
nb <- dnearneigh(coords, d1 = 0, d2 = 1.5)
listw <- nb2listw(nb, style = "W", zero.policy = TRUE)

moran_label <- function(mod) {
  mt <- moran.test(residuals(mod, type = "deviance"), listw, zero.policy = TRUE)
  i_val <- unname(mt$estimate["Moran I statistic"])
  p_val <- mt$p.value
  p_str <- if (p_val < 0.001) "< 0.001" else sprintf("= %.3f", p_val)
  sprintf("Moran's I = %.3f (p %s)", i_val, p_str)
}

moran_full   <- moran_label(bundle$m_full)
moran_subset <- moran_label(bundle$m_subset)
moran_phylo  <- moran_label(bundle$m_phylo)

get_smooth_curve <- function(mod, term_name, var_name) {
  dat <- smooth_estimates(mod, select = term_name)
  x_val <- if (".value" %in% names(dat)) dat$.value else dat[[var_name]]
  y_val <- if (".estimate" %in% names(dat)) dat$.estimate else dat$partial
  x_scaled <- (x_val - min(x_val)) / (max(x_val) - min(x_val))
  data.frame(Variable = var_name, X_Relative = x_scaled, Effect = y_val)
}

get_linear_curve <- function(mod, data, var_name) {
  var_seq <- seq(min(data[[var_name]], na.rm = TRUE), max(data[[var_name]], na.rm = TRUE), length.out = 100)
  pred_df <- data[1:100, ]
  pred_df[[var_name]] <- var_seq
  terms_out <- predict(mod, newdata = pred_df, type = "terms")
  data.frame(
    Variable = var_name,
    X_Relative = seq(0, 1, length.out = 100),
    Effect = terms_out[, var_name]
  )
}

build_spotlight_data <- function(mod, data, target_group = "built", linear_vars = character(0)) {
  curve_bio4  <- get_smooth_curve(mod, "s(Bio_4)", "Bio_4")
  curve_bio15 <- get_smooth_curve(mod, "s(Bio_15)", "Bio_15")
  curve_elev  <- get_smooth_curve(mod, "s(Elevation)", "Elevation")
  curve_hii   <- get_linear_curve(mod, data, "HII")

  extra_curves <- lapply(linear_vars, function(v) get_linear_curve(mod, data, v))

  all_curves <- bind_rows(curve_bio4, curve_bio15, curve_elev, curve_hii, extra_curves)

  lu_levels <- levels(data$Landuse)
  pred_df_lu <- data[1:length(lu_levels), ]
  pred_df_lu$Landuse <- factor(lu_levels, levels = lu_levels)
  lu_terms <- predict(mod, newdata = pred_df_lu, type = "terms")
  landuse_eff <- data.frame(Landuse = lu_levels, Shift = lu_terms[, "Landuse"])
  ref_shift <- landuse_eff$Shift[landuse_eff$Landuse == "cropland"]
  landuse_eff$Shift <- landuse_eff$Shift - ref_shift

  plot_data <- tidyr::crossing(all_curves, landuse_eff) %>%
    mutate(
      Shifted_Effect = Effect + Shift,
      Is_Target = (Landuse == target_group)
    )
  plot_data
}

theme_pub <- function() {
  theme_minimal(base_size = 10, base_family = PUB_FONT) +
    theme(
      plot.title       = element_text(size = 10.5, face = "plain", family = PUB_FONT, margin = margin(b = 2), lineheight = 1.05),
      plot.subtitle    = element_text(size = 9, family = PUB_FONT, colour = "grey30", margin = margin(b = 4)),
      axis.title       = element_text(size = 10, family = PUB_FONT),
      axis.text        = element_text(size = 8.5, family = PUB_FONT, colour = "black"),
      legend.text      = element_text(size = 9, family = PUB_FONT),
      legend.title     = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, colour = "grey88"),
      axis.line        = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks       = element_line(linewidth = 0.3, colour = "black"),
      plot.margin      = margin(6, 8, 6, 6)
    )
}

VAR_COLORS <- c("Bio_4" = "#d73027", "Bio_15" = "#4575b4", "Elevation" = "#1b7837",
                 "HII" = "#762a83", "NRI" = "#e08214", "NTI" = "#35978f")

make_panel <- function(plot_data, title, subtitle, y_limits, end_label_vars = character(0)) {
  target_data <- plot_data |> filter(Is_Target)
  target_main <- target_data |> filter(!Variable %in% end_label_vars)
  target_end  <- target_data |> filter(Variable %in% end_label_vars)

  p <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_line(
      data = plot_data |> filter(!Is_Target),
      aes(x = X_Relative, y = Shifted_Effect, group = interaction(Variable, Landuse)),
      color = "grey60", alpha = 0.3, linewidth = 0.45
    ) +
    geom_line(
      data = target_main,
      aes(x = X_Relative, y = Shifted_Effect, color = Variable),
      linewidth = 1.0
    )

  p <- p + scale_color_manual(values = VAR_COLORS)

  if (nrow(target_end) > 0) {
    for (v in unique(target_end$Variable)) {
      p <- p + geom_line(
        data = target_end |> filter(Variable == v),
        aes(x = X_Relative, y = Shifted_Effect),
        color = VAR_COLORS[[v]], linewidth = 1.0
      )
    }

    end_points <- target_end |>
      group_by(Variable) |>
      filter(X_Relative == max(X_Relative)) |>
      ungroup()
    p <- p +
      geom_text_repel(
        data = end_points,
        aes(x = X_Relative, y = Shifted_Effect, label = Variable),
        color = VAR_COLORS[end_points$Variable],
        hjust = 0, nudge_x = 0.03, direction = "y", segment.size = 0.3,
        min.segment.length = 0, size = 3, fontface = "bold", family = PUB_FONT,
        seed = 1
      ) +
      scale_x_continuous(expand = expansion(mult = c(0.02, 0.14)))
  } else {
    p <- p + scale_x_continuous(expand = expansion(mult = 0.02))
  }

  p +
    coord_cartesian(ylim = y_limits, clip = "off") +
    theme_pub() +
    theme(legend.position = "none") +
    labs(title = title, subtitle = subtitle,
         x = "Relative gradient (0 to 1)", y = "Partial effect (intercept shifted)")
}

pd_full   <- build_spotlight_data(bundle$m_full, model_df)
pd_subset <- build_spotlight_data(bundle$m_subset, model_df)
pd_phylo  <- build_spotlight_data(bundle$m_phylo, model_df)

fmt_coef <- function(mod, term) {
  pt <- summary(mod)$p.table
  est <- pt[term, "Estimate"]
  p_val <- pt[term, "Pr(>|t|)"]
  p_str <- if (p_val < 0.001) "< 0.001" else sprintf("= %.3f", p_val)
  sprintf("%s: b = %.3f, p %s", term, est, p_str)
}
nri_nti_note <- paste(fmt_coef(bundle$m_phylo, "NRI"), fmt_coef(bundle$m_phylo, "NTI"), sep = "\n")

y_ab <- range(c(pd_full$Shifted_Effect, pd_subset$Shifted_Effect))
y_c  <- range(pd_phylo$Shifted_Effect)

p_full   <- make_panel(pd_full,   "Full pool (~400 sp.)\nno phylogeny", moran_full,   y_ab)
p_subset <- make_panel(pd_subset, "247-sp. tree subset\nno phylogeny", moran_subset, y_ab)
p_phylo  <- make_panel(pd_phylo,  "247-sp. subset + NRI/NTI", moran_phylo, y_c) +
  theme(legend.position = "bottom") +
  geom_label(
    data = data.frame(x = Inf, y = Inf, label = nri_nti_note),
    aes(x = x, y = y, label = label), inherit.aes = FALSE,
    hjust = 1, vjust = 1, size = 3, family = PUB_FONT, lineheight = 1.2,
    linewidth = 0, fill = "white", alpha = 0.85
  )

final_figure <- ((p_full | p_subset) / p_phylo) +
  plot_layout(guides = "collect", heights = c(1, 1.3)) &
  theme(legend.position = "bottom")

final_figure <- final_figure +
  plot_annotation(
    title = "Relative drivers of global butterfly richness:\nphylogenetic sensitivity test",
    subtitle = "Highlighted group: built vs. background land-use categories (grey)\n1° grid, n = 4,792 cells | NRI/NTI via picante::ses.mpd/ses.mntd (null = 'richness', 199 runs)",
    tag_levels = "A"
  ) &
  theme(
    text = element_text(family = PUB_FONT),
    plot.title = element_text(size = 13, face = "bold", family = PUB_FONT, lineheight = 1.1),
    plot.subtitle = element_text(size = 8.5, family = PUB_FONT, colour = "grey30", lineheight = 1.15),
    plot.tag = element_text(size = 12, face = "bold", family = PUB_FONT)
  )

dir.create("output/BAM/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("output/phylo_export", showWarnings = FALSE, recursive = TRUE)

FIG_WIDTH_MM  <- 170
FIG_HEIGHT_MM <- 235

ggsave(
  "output/BAM/figures/Support_Figure_1__Main_Spotlight.png", final_figure,
  width = FIG_WIDTH_MM, height = FIG_HEIGHT_MM, units = "mm",
  dpi = 600, bg = "white", device = ragg::agg_png
)

ggsave(
  "output/phylo_export/figure4d_phylo_sensitivity_spotlight_comparison.pdf", final_figure,
  width = FIG_WIDTH_MM, height = FIG_HEIGHT_MM, units = "mm",
  device = "pdf"
)

cat("Saved: output/BAM/figures/Support_Figure_1__Main_Spotlight.png (600 dpi, Helvetica, replaced)\n")
cat("Saved: output/phylo_export/figure4d_phylo_sensitivity_spotlight_comparison.pdf (vector)\n")
