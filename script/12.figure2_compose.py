#!/usr/bin/env python3
# Compose Figure 2 panels into final figure.

from pathlib import Path
from PIL import Image

BASE = Path("output/Manuscript/reproducibility_code")
PANEL_DIR = BASE / "panels"
GUTTER = 24

names = ["panel_a", "panel_b", "panel_c", "panel_d"]
imgs = [Image.open(PANEL_DIR / f"{n}.png").convert("RGB") for n in names]

w, h = imgs[0].size
assert all(im.size == (w, h) for im in imgs), "Panels must all be the same size " \
    "(set consistent width/height in ggsave() in script/11.figure2_final.R)"

grid_w = w * 2 + GUTTER
grid_h = h * 2 + GUTTER

canvas = Image.new("RGB", (grid_w, grid_h), "white")
canvas.paste(imgs[0], (0, 0))
canvas.paste(imgs[1], (w + GUTTER, 0))
canvas.paste(imgs[2], (0, h + GUTTER))
canvas.paste(imgs[3], (w + GUTTER, h + GUTTER))

out_png = BASE / "Figure2_final.png"
out_pdf = BASE / "Figure2_final.pdf"

canvas.save(out_png, dpi=(400, 400))
canvas.save(out_pdf, resolution=400.0)

print(f"Saved: {out_png} ({canvas.size[0]}x{canvas.size[1]} px)")
print(f"Saved: {out_pdf}")
