# Figure 3: assemble tree, family labels, butterfly photos, histogram inset and colour bar as vector
# artwork at the Nature double-column width (183 x 170 mm; text 7 pt). Photos are embedded rasters.

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggtree)
  library(ggtreeExtra)
  library(png)
  library(grid)
})

out_dir <- "output/Manuscript/reproducibility_code"
img_dir <- "updatedata/butterfly_family_images"
parts <- readRDS(file.path(out_dir, "Figure3_parts.rds"))
svglite::svglite(tempfile(fileext = ".svg"))  # measure grobs on a device that knows Arial

FIG_W <- 183; FIG_H <- 170        # mm
BASE_SIZE <- 7
R_MM <- 64                         # radius of the family ring on the page
TEXT_R <- 1.06; IMG_R <- 1.23; IMG_SIZE <- 0.23   # label/photo radius and photo size (~15 mm), relative to R

# Vector PDF with embedded Arial: SVG (svglite) converted by librsvg, since R's Cairo PDF device is unavailable here
save_vector_pdf <- function(plot, file, width, height) {
  svg <- tempfile(fileext = ".svg")
  ggsave(svg, plot, width = width, height = height, units = "mm", device = svglite::svglite)
  rsvg::rsvg_pdf(svg, file)
}

# Set3 colours used for the family ring
fam_col <- c(Hesperiidae = "#8DD3C7", Lycaenidae = "#FFFFB3", Nymphalidae = "#BEBADA",
             Papilionidae = "#FB8072", Pieridae = "#80B1D3")

# Locate the ring centre, radius and family mid-angles on a test render of the tree
probe_in <- 6
probe <- tempfile(fileext = ".png")
ggsave(probe, parts$tree, width = probe_in, height = probe_in, dpi = 300, bg = "white")
a <- readPNG(probe)[, , 1:3]
ink <- which(a[, , 1] < 0.98 | a[, , 2] < 0.98 | a[, , 3] < 0.98, arr.ind = TRUE)
x0 <- min(ink[, 2]); x1 <- max(ink[, 2]); y0 <- min(ink[, 1]); y1 <- max(ink[, 1])
cx <- (x0 + x1) / 2; cy <- (y0 + y1) / 2
r_px <- ((x1 - x0) / 2 + (y1 - y0) / 2) / 2
cols <- t(col2rgb(fam_col)) / 255
theta <- seq(0, 2 * pi, length.out = 2001)[-2001]
px <- round(cx + 0.985 * r_px * cos(theta)); py <- round(cy + 0.985 * r_px * sin(theta))
fam_at <- sapply(seq_along(theta), function(i) {
  rgb <- a[py[i], px[i], ]
  d <- colSums((t(cols) - rgb)^2)
  if (min(d) < 3 * (40 / 255)^2) names(fam_col)[which.min(d)] else NA
})
fam_angle <- sapply(names(fam_col), function(f) {
  t <- theta[fam_at %in% f]
  -atan2(mean(sin(t)), mean(cos(t)))   # image rows run downwards; page y runs upwards
})

# Scale the tree so that the ring radius is R_MM, and place the ring centre on the page
tree_mm <- R_MM * (probe_in * 300) / r_px
ctr <- c(FIG_W / 2 - 4, FIG_H / 2)
tx0 <- ctr[1] - tree_mm * cx / (probe_in * 300)
ty1 <- ctr[2] + tree_mm * cy / (probe_in * 300)

fig <- ggplot() +
  coord_fixed(xlim = c(0, FIG_W), ylim = c(0, FIG_H), expand = FALSE) +
  theme_void() +
  annotation_custom(ggplotGrob(parts$tree), xmin = tx0, xmax = tx0 + tree_mm,
                    ymin = ty1 - tree_mm, ymax = ty1)

for (f in names(fam_angle)) {
  th <- fam_angle[[f]]
  # photo
  img <- readPNG(file.path(img_dir, paste0(f, ".png")))
  long <- R_MM * IMG_SIZE
  w <- if (ncol(img) >= nrow(img)) long else long * ncol(img) / nrow(img)
  h <- if (ncol(img) >= nrow(img)) long * nrow(img) / ncol(img) else long
  ix <- ctr[1] + R_MM * IMG_R * cos(th); iy <- ctr[2] + R_MM * IMG_R * sin(th)
  fig <- fig + annotation_raster(img, xmin = ix - w / 2, xmax = ix + w / 2, ymin = iy - h / 2, ymax = iy + h / 2)
  # label, tangent to the ring and kept upright
  rot <- (th * 180 / pi + 90) %% 360
  if (rot > 90 && rot <= 270) rot <- rot - 180
  fig <- fig + annotate("text", x = ctr[1] + R_MM * TEXT_R * cos(th), y = ctr[2] + R_MM * TEXT_R * sin(th),
                        label = f, angle = rot, size = BASE_SIZE / .pt, family = "Arial")
}

# Histogram inset in the open centre and colour bar at the top right, next to the ring
hw <- 44; hh <- 30   # mm; fills the open centre of the tree
fig <- fig + annotation_custom(ggplotGrob(parts$hist), xmin = ctr[1] + 0.05 * R_MM - hw / 2,
                               xmax = ctr[1] + 0.05 * R_MM + hw / 2, ymin = ctr[2] - hh / 2, ymax = ctr[2] + hh / 2)
lw <- 1.52 * 25.4; lh <- 0.27 * 25.4
LEG_Y <- ctr[2] + 0.83 * R_MM   # top right, level with the upper part of the ring
fig <- fig + annotation_custom(ggplotGrob(parts$legend), xmin = FIG_W - lw - 1, xmax = FIG_W - 1,
                               ymin = LEG_Y, ymax = LEG_Y + lh)

invisible(dev.off())
save_vector_pdf(fig, file.path(out_dir, "Figure3_final.pdf"), width = FIG_W, height = FIG_H)
ggsave(file.path(out_dir, "Figure3_final.png"), fig, width = FIG_W, height = FIG_H, units = "mm", dpi = 600,
       device = ragg::agg_png, bg = "white")
cat("Saved Figure 3 (PDF and 600 dpi PNG) to", out_dir, "\n")
