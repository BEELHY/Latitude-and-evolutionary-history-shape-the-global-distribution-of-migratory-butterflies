#!/usr/bin/env python3
# Compose Figure 3.

import os
from pathlib import Path
import math
from PIL import Image, ImageDraw, ImageFont
import numpy as np

BASE = Path("output/Manuscript/reproducibility_code")
BASE_PNG = BASE / "Figure3_base.png"
OUT_PNG = BASE / "Figure3_final.png"
IMG_DIR = Path("updatedata/butterfly_family_images")
# Arial Bold, first found
FONT_CANDIDATES = [
    os.environ.get("FIG_FONT", ""),
    "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
    "/usr/share/fonts/truetype/msttcorefonts/Arial_Bold.ttf",
    "C:/Windows/Fonts/arialbd.ttf",
]
FONT_PATH = next((f for f in FONT_CANDIDATES if f and os.path.exists(f)), None)

# Set3 colours by family
FAMILY_COLORS = {
    "Hesperiidae": (0x8D, 0xD3, 0xC7),
    "Lycaenidae":  (0xFF, 0xFF, 0xB3),
    "Nymphalidae": (0xBE, 0xBA, 0xDA),
    "Papilionidae": (0xFB, 0x80, 0x72),
    "Pieridae":    (0x80, 0xB1, 0xD3),
}

# Find tree centre and radius
img = Image.open(BASE_PNG).convert("RGB")
arr = np.array(img)
W, H = img.size

mask = np.ones((H, W), dtype=bool)

non_white = np.any(arr < 250, axis=2) & mask
ys, xs = np.nonzero(non_white)
new_x0, new_x1 = xs.min(), xs.max()
new_y0, new_y1 = ys.min(), ys.max()
center = ((new_x0 + new_x1) / 2, (new_y0 + new_y1) / 2)
R = ((new_x1 - new_x0) / 2 + (new_y1 - new_y0) / 2) / 2
print(f"Tree bbox: x[{new_x0},{new_x1}] y[{new_y0},{new_y1}] center={center} R={R:.1f}")

# Sample family ring
SAMPLE_RADIUS = R * 0.985  # Inside outer edge
N_SAMPLES = 2000

def nearest_family(rgb):
    best, best_d = None, 1e9
    for fam, c in FAMILY_COLORS.items():
        d = sum((int(rgb[i]) - c[i]) ** 2 for i in range(3))
        if d < best_d:
            best, best_d = fam, d
    return best if best_d < 40 ** 2 * 3 else None  # Not a ring colour

angle_lists = {fam: [] for fam in FAMILY_COLORS}
for i in range(N_SAMPLES):
    theta = 2 * math.pi * i / N_SAMPLES
    px = int(round(center[0] + SAMPLE_RADIUS * math.cos(theta)))
    py = int(round(center[1] + SAMPLE_RADIUS * math.sin(theta)))
    if 0 <= px < W and 0 <= py < H:
        fam = nearest_family(arr[py, px])
        if fam is not None:
            angle_lists[fam].append(theta)

fam_angle = {}
for fam, angles in angle_lists.items():
    if not angles:
        print(f"WARNING: no ring pixels matched to {fam}")
        continue
    s = sum(math.sin(a) for a in angles) / len(angles)
    c = sum(math.cos(a) for a in angles) / len(angles)
    fam_angle[fam] = math.atan2(s, c)
    print(f"{fam}: {len(angles)} samples, mid-angle = {math.degrees(fam_angle[fam]):.1f} deg")

# Placement
TEXT_R_FRAC = 1.06  # Labels outside ring
IMG_R_FRAC = 1.13  # Photos beyond labels
IMG_SIZE_FRAC = 0.075  # Photo size vs radius
DPI = 600
FONT_SIZE_PX = round(12 / 72 * DPI)  # 12 pt

PAD = int(R * 0.35)
canvas = Image.new("RGBA", (W + 2 * PAD, H + 2 * PAD), (255, 255, 255, 255))
canvas.paste(img, (PAD, PAD))
center = (center[0] + PAD, center[1] + PAD)

font = ImageFont.truetype(FONT_PATH, FONT_SIZE_PX) if FONT_PATH else ImageFont.load_default(FONT_SIZE_PX)

# Histogram in centre
hist = Image.open(BASE / "Figure3_hist.png").convert("RGBA")
hw, hh = hist.size
canvas.alpha_composite(hist, (int(center[0] - hw / 2), int(center[1] - hh / 2)))

for fam, a in fam_angle.items():
    deg = math.degrees(a)
    # Tangent label angle
    text_deg = (-deg + 90) % 360
    if text_deg > 90 and text_deg <= 270:
        text_deg -= 180  # Keep text upright

    # Photo
    icx = center[0] + R * IMG_R_FRAC * math.cos(a)
    icy = center[1] + R * IMG_R_FRAC * math.sin(a)
    photo = Image.open(IMG_DIR / f"{fam}.png").convert("RGBA")
    ow, oh = photo.size
    long_side = R * IMG_SIZE_FRAC
    if ow >= oh:
        iw, ih = long_side, long_side * oh / ow
    else:
        ih, iw = long_side, long_side * ow / oh
    photo = photo.resize((max(1, int(iw)), max(1, int(ih))), Image.LANCZOS)
    canvas.alpha_composite(photo, (int(icx - iw / 2), int(icy - ih / 2)))

    # Label
    tx = center[0] + R * TEXT_R_FRAC * math.cos(a)
    ty = center[1] + R * TEXT_R_FRAC * math.sin(a)
    txt_layer = Image.new("RGBA", (8 * FONT_SIZE_PX, 2 * FONT_SIZE_PX), (0, 0, 0, 0))
    d = ImageDraw.Draw(txt_layer)
    d.text((0, 0), fam, font=font, fill=(0, 0, 0, 255))
    txt_layer = txt_layer.crop(txt_layer.getbbox())
    rotated = txt_layer.rotate(text_deg, expand=True, resample=Image.BICUBIC)
    canvas.alpha_composite(rotated, (int(tx - rotated.width / 2), int(ty - rotated.height / 2)))

# Trim margins
final = canvas.convert("RGB")
_a = np.array(final); _nz = np.any(_a < 250, axis=2); _y, _x = np.nonzero(_nz)
final = final.crop((_x.min(), _y.min(), _x.max() + 1, _y.max() + 1))

# Legend in free corner
leg = Image.open(BASE / "Figure3_legend.png").convert("RGBA")
_la = np.array(leg); _ly, _lx = np.nonzero(_la[:, :, 3] > 0)
leg = leg.crop((_lx.min(), _ly.min(), _lx.max() + 1, _ly.max() + 1))
lw, lh = leg.size
M = 50  # Clearance
occ = np.any(np.array(final) < 250, axis=2)
occ = np.pad(occ, ((0, max(0, lh + M - occ.shape[0])), (0, lw + M)), constant_values=False)
x_best = final.width
for x in range(final.width, -1, -5):
    if occ[0:lh + M, max(0, x - M):x + lw].any():
        break
    x_best = x
new_w = max(final.width, x_best + lw)
fin = Image.new("RGBA", (new_w, final.height), (255, 255, 255, 255))
fin.paste(final, (0, 0))
fin.alpha_composite(leg, (x_best, 0))
final = fin.convert("RGB")

# White border
B = 80
bordered = Image.new("RGB", (final.width + 2 * B, final.height + 2 * B), "white")
bordered.paste(final, (B, B))
final = bordered
final_arr = np.array(final)
nz = np.any(final_arr < 250, axis=2)
ys2, xs2 = np.nonzero(nz)
b = 60
crop_box = (max(xs2.min() - b, 0), max(ys2.min() - b, 0),
            min(xs2.max() + b, final.width), min(ys2.max() + b, final.height))
final = final.crop(crop_box)

final.save(OUT_PNG, dpi=(DPI, DPI))
final.save(BASE / "Figure3_final.pdf", resolution=DPI)
print(f"Saved: {OUT_PNG} ({final.width}x{final.height})")
