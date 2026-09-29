# Figure 1: seasonal-switching richness (Robinson maps, latitude profile, seasonal net change).

import os
import sys
import glob
import re
import warnings
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.warp import reproject, Resampling
from rasterio.transform import from_bounds
import matplotlib.pyplot as plt
import matplotlib as mpl
import cartopy.crs as ccrs

warnings.filterwarnings('ignore')

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial'],
    'font.size': 7,
    'axes.labelsize': 7,
    'axes.titlesize': 7,
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'axes.linewidth': 0.5,
    'xtick.major.width': 0.5,
    'ytick.major.width': 0.5,
    'xtick.major.size': 2.5,
    'ytick.major.size': 2.5,
    'pdf.fonttype': 42,
    'savefig.facecolor': 'white',
    'figure.facecolor': 'white'
})

# Paths (raw maps not shared)
shp_path = "data/world-administrative-boundaries/world-administrative-boundaries.shp"
csv_path = "output/df_lat_sp.csv"
data_dir = "data/SuitabilityMaps_MigratorySpecies"
output_dir = "output/Manuscript/figure1_panels"
os.makedirs(output_dir, exist_ok=True)
cache = sys.argv[1] if len(sys.argv) > 1 else None  # optional .npz of precomputed grids

if cache and os.path.exists(cache):
    z = np.load(cache)
    migratory_richness_map = z['rich'].astype('float32')
    S = {s: z[f'S{s}'].astype('float32') for s in [1, 2, 3, 4]}
else:
    df_ref = pd.read_csv(csv_path)
    keep = set(df_ref.species)
    all_files = glob.glob(os.path.join(data_dir, "Binary_S*.tif"))
    species_dict = {}
    for f in all_files:
        m = re.match(r"Binary_S(\d)(.*)\.tif$", os.path.basename(f))
        if m and m.group(2) in keep:
            species_dict.setdefault(m.group(2), {})[int(m.group(1))] = f

    # Native-resolution master grid
    with rasterio.open(all_files[0]) as src:
        height, width = int(180 / src.res[1]), int(360 / src.res[0])
    master_transform = from_bounds(-180, -90, 180, 90, width, height)

    migratory_richness_map = np.zeros((height, width), dtype='float32')
    S = {s: np.zeros((height, width), dtype='float32') for s in [1, 2, 3, 4]}
    for sp, seasons in species_dict.items():
        t_max = len(seasons)
        temp = {s: np.zeros((height, width), dtype='float32') for s in [1, 2, 3, 4]}
        for s, path in seasons.items():
            with rasterio.open(path) as src:
                reproject(src.read(1), temp[s], src_transform=src.transform, src_crs=src.crs,
                          dst_transform=master_transform, dst_crs=src.crs, resampling=Resampling.nearest)
        p_sum = sum(temp.values())
        mig_mask = np.where((p_sum > 0) & (p_sum < t_max), 1, 0).astype('float32')
        if np.sum(mig_mask) == 0:
            continue
        migratory_richness_map += mig_mask
        for s in [1, 2, 3, 4]:
            S[s] += (temp[s] * mig_mask)

height = migratory_richness_map.shape[0]
max_richness = int(np.max(migratory_richness_map))
world = gpd.read_file(shp_path).dissolve()

# Seasonal net change: continuous viridis on a symmetric log scale so that small changes stay visible
labels = [("S2–S1", S[2] - S[1]), ("S3–S2", S[3] - S[2]), ("S4–S3", S[4] - S[3]), ("S1–S4", S[1] - S[4])]
c_min = int(min(d.min() for _, d in labels)); c_max = int(max(d.max() for _, d in labels))
print(f"Maximum richness: {max_richness}; net change range {c_min} to {c_max}")
c_norm = mpl.colors.SymLogNorm(linthresh=1, linscale=0.5, vmin=c_min, vmax=c_max, base=10)

robinson = ccrs.Robinson(central_longitude=0)
pc = ccrs.PlateCarree()
LAND = '#e8e8e8'
P_STEP = 2  # every 2nd cell (5 arc-min) for rendering
FS = 7      # Nature: 5-7 pt text at the final width of 183 mm
FS_TAG = 8  # panel labels 8 pt bold


def base_map(ax, equator=True):
    ax.set_extent([-180, 180, -60, 90], crs=pc)
    ax.spines['geo'].set_visible(False)
    ax.add_geometries(world.geometry, crs=pc, facecolor=LAND, edgecolor='none', zorder=1)
    if equator:
        ax.plot([-180, 180], [0, 0], transform=pc, color='0.35', lw=0.4, ls=(0, (3, 3)), zorder=3)


def raster(ax, arr, **kw):
    return ax.imshow(arr[::P_STEP, ::P_STEP], origin='upper', extent=[-180, 180, -90, 90],
                     transform=pc, interpolation='nearest', regrid_shape=2400, zorder=2, **kw)


def title(x, y, letter, text):
    fig.text(x, y, letter, fontsize=FS_TAG, fontweight='bold', va='top', ha='left')
    fig.text(x + 0.13 / FW, y, text, fontsize=FS, va='top', ha='left', linespacing=1.15)


def end_bar(cax, mappable, lo, hi, label=None, mid=None):
    """Horizontal colour bar without ticks; minimum and maximum written at its two ends."""
    cb = fig.colorbar(mappable, cax=cax, orientation='horizontal', ticks=[])
    cb.ax.minorticks_off()
    cb.ax.xaxis.set_ticks([], minor=True)
    cb.outline.set_linewidth(0.4)
    cax.text(-0.03, 0.5, lo, transform=cax.transAxes, ha='right', va='center', fontsize=FS)
    cax.text(1.03, 0.5, hi, transform=cax.transAxes, ha='left', va='center', fontsize=FS)
    if label:
        cax.text(0.5, 1.5, label, transform=cax.transAxes, ha='center', va='bottom', fontsize=FS)
    if mid is not None:
        x = cb.norm(mid[0]) if not callable(getattr(cb.norm, 'inverse', None)) else cb.norm(mid[0])
        cax.plot([x, x], [0, 1], transform=cax.transAxes, color='white', lw=0.5)
        cax.text(x, -0.4, mid[1], transform=cax.transAxes, ha='center', va='top', fontsize=FS)
    return cb


# Layout in inches from the top-left corner, following slide 1 of Figure1_final.pptx, at the Nature
# double-column width (183 mm; maximum height 170 mm).
# Maps span 60 S - 90 N; panel b shares the vertical extent of map a, so latitudes line up.
FW, FH = 183 / 25.4, 6.62  # 168 mm tall
fig = plt.figure(figsize=(FW, FH))
ROB_W = 2 * robinson.transform_point(180, 0, pc)[0]
rob_y = lambda lat: robinson.transform_point(0, lat, pc)[1]
MAP_ASPECT = ROB_W / (rob_y(90) - rob_y(-60))


def box(x, y, w, h, **kw):
    return fig.add_axes([x / FW, 1 - (y + h) / FH, w / FW, h / FH], **kw)


def at(x, y):
    return x / FW, 1 - y / FH


A_TOP, A_W = 0.34, 5.0
A_H = A_W / MAP_ASPECT
B_X, B_W = 5.5, 1.62
title(*at(0.02, 0.04), 'a', 'Seasonal-switching richness')
ax_a = box(0.0, A_TOP, A_W, A_H, projection=robinson)
cax_a = box(3.3, 0.16, 1.0, 0.07)
title(*at(B_X - 0.55, 0.04), 'b', 'Mean seasonal-switching\nrichness (n species)')
ax_b = box(B_X, A_TOP, B_W, A_H)
C_TOP = A_TOP + A_H + 0.3
title(*at(0.02, C_TOP), 'c', 'Net species flux (gain − loss)')
cax_c = box(3.3, C_TOP + 0.12, 1.0, 0.07)  # titled like the colour bar of a
MAP_W = 3.55
MAP_H = MAP_W / MAP_ASPECT
ROW1 = C_TOP + 0.42
ROW2 = ROW1 + MAP_H + 0.2
axes_c = [box(x, y, MAP_W, MAP_H, projection=robinson) for y in (ROW1, ROW2) for x in (0.0, FW - MAP_W)]

# a
base_map(ax_a)
im_a = raster(ax_a, np.ma.masked_equal(migratory_richness_map, 0), cmap='viridis', vmin=1, vmax=max_richness)
end_bar(cax_a, im_a, '1', str(max_richness), label='Species (n)')

# b: mean richness across occupied cells within each latitude band, plotted on the Robinson
# vertical coordinate so that each latitude sits level with map a
lats = np.linspace(90, -90, height)
lat_mean = np.array([r[r > 0].mean() if (r > 0).any() else 0.0 for r in migratory_richness_map])
ys = np.array([rob_y(l) for l in lats])
ax_b.fill_betweenx(ys, 0, lat_mean, color=mpl.cm.viridis(0.55), alpha=0.35, lw=0)
ax_b.plot(lat_mean, ys, color=mpl.cm.viridis(0.3), lw=0.7)
ax_b.axhline(0, color='0.35', lw=0.4, ls=(0, (3, 3)))
ax_b.set_ylim(rob_y(-60), rob_y(90))
ax_b.set_yticks([rob_y(90), 0, rob_y(-60)])
ax_b.set_yticklabels(['90°', '0°', '−60°'])
ax_b.get_yticklabels()[1].set_va('top')  # sits just below the equator line, as in the slide
ax_b.set_xlim(0, np.ceil(lat_mean.max() / 5) * 5)
ax_b.set_xticks([0, 15])
for side in ('top', 'right'):
    ax_b.spines[side].set_visible(False)
ax_b.tick_params(labelsize=FS)
# equator carried across the gap between a and b
y_eq = 1 - (A_TOP + A_H * rob_y(90) / (rob_y(90) - rob_y(-60))) / FH
fig.add_artist(mpl.lines.Line2D([A_W / FW, B_X / FW], [y_eq, y_eq], color='0.35', lw=0.4, ls=(0, (3, 3))))

# c
for ax, (lab, d) in zip(axes_c, labels):
    base_map(ax, equator=False)
    raster(ax, np.ma.masked_equal(d, 0), cmap='viridis', norm=c_norm)
    ax.set_title(lab, pad=2, fontsize=FS)
sm = mpl.cm.ScalarMappable(norm=c_norm, cmap='viridis')
end_bar(cax_c, sm, f'−{abs(c_min)}', str(c_max), label='Species (n)', mid=(0, '0'))

for ext in ('png', 'pdf'):
    fig.savefig(os.path.join(output_dir, f"Figure1.{ext}"), dpi=600)
plt.close(fig)
print(f"Saved {output_dir}/Figure1.png/.pdf")
