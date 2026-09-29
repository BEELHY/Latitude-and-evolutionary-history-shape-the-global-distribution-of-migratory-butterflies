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
    'xtick.labelsize': 6.5,
    'ytick.labelsize': 6.5,
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

# Seasonal net change, shown in classes so that small and large changes are both visible
labels = [("S2–S1", S[2] - S[1]), ("S3–S2", S[3] - S[2]), ("S4–S3", S[4] - S[3]), ("S1–S4", S[1] - S[4])]
c_min = int(min(d.min() for _, d in labels)); c_max = int(max(d.max() for _, d in labels))
print(f"Maximum richness: {max_richness}; net change range {c_min} to {c_max}")
c_bounds = [c_min - 0.5, -10.5, -5.5, -2.5, -0.5, 0.5, 2.5, 5.5, 10.5, c_max + 0.5]
c_cols = [mpl.cm.PuOr_r(x) for x in (0.0, 0.13, 0.26, 0.38, 0.5, 0.62, 0.74, 0.87, 1.0)]
c_cols[4] = '#e8e8e8'  # no change: shown as land
c_cmap = mpl.colors.ListedColormap(c_cols)
c_norm = mpl.colors.BoundaryNorm(c_bounds, c_cmap.N)

robinson = ccrs.Robinson(central_longitude=0)
pc = ccrs.PlateCarree()
LAND = '#e8e8e8'
P_STEP = 2  # every 2nd cell (5 arc-min) for rendering


def base_map(ax):
    ax.set_global()
    ax.spines['geo'].set_visible(False)
    ax.add_geometries(world.geometry, crs=pc, facecolor=LAND, edgecolor='none', zorder=1)
    ax.plot([-180, 180], [0, 0], transform=pc, color='0.35', lw=0.4, ls=(0, (3, 3)), zorder=3)


def raster(ax, arr, **kw):
    return ax.imshow(arr[::P_STEP, ::P_STEP], origin='upper', extent=[-180, 180, -90, 90],
                     transform=pc, interpolation='nearest', regrid_shape=2400, zorder=2, **kw)


def letter(ax, s, x=-0.02, y=1.0):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=9, fontweight='bold', va='bottom', ha='left')


MM = 1 / 25.4
fig = plt.figure(figsize=(180 * MM, 170 * MM))
ax_a = fig.add_axes([0.00, 0.575, 0.745, 0.42], projection=robinson)
cax_a = fig.add_axes([0.035, 0.65, 0.20, 0.011])
ax_b = fig.add_axes([0.815, 0.625, 0.17, 0.35])
pos_c = [[0.00, 0.325, 0.495, 0.235], [0.505, 0.325, 0.495, 0.235],
         [0.00, 0.075, 0.495, 0.235], [0.505, 0.075, 0.495, 0.235]]
axes_c = [fig.add_axes(p, projection=robinson) for p in pos_c]
cax_c = fig.add_axes([0.17, 0.055, 0.66, 0.011])

# a
base_map(ax_a)
im_a = raster(ax_a, np.ma.masked_equal(migratory_richness_map, 0), cmap='viridis', vmin=1, vmax=max_richness)
letter(ax_a, 'a', 0.0, 0.93)
cb_a = fig.colorbar(im_a, cax=cax_a, orientation='horizontal', ticks=[1, 20, 40, max_richness])
cb_a.set_label('Seasonal-switching richness\n(number of species)', labelpad=2)
cb_a.ax.xaxis.set_label_position('top')
cb_a.outline.set_linewidth(0.4)

# b: mean richness across occupied cells within each latitude band
lats = np.linspace(90, -90, height)
lat_mean = np.array([r[r > 0].mean() if (r > 0).any() else 0.0 for r in migratory_richness_map])
ax_b.fill_betweenx(lats, 0, lat_mean, color=mpl.cm.viridis(0.55), alpha=0.35, lw=0)
ax_b.plot(lat_mean, lats, color=mpl.cm.viridis(0.3), lw=0.7)
ax_b.axhline(0, color='0.35', lw=0.4, ls=(0, (3, 3)))
ax_b.set_ylim(-60, 90)
ax_b.set_yticks([-60, -30, 0, 30, 60, 90])
ax_b.set_yticklabels(['60° S', '30° S', '0°', '30° N', '60° N', '90° N'])
ax_b.set_xlim(0, np.ceil(lat_mean.max() / 5) * 5)
ax_b.set_xlabel('Mean richness per\noccupied cell (species)')
for side in ('top', 'right'):
    ax_b.spines[side].set_visible(False)
letter(ax_b, 'b', -0.5, 1.0)

# c
for ax, (lab, d) in zip(axes_c, labels):
    base_map(ax)
    raster(ax, np.ma.masked_equal(d, 0), cmap=c_cmap, norm=c_norm)
    ax.set_title(lab, pad=1)
letter(axes_c[0], 'c', 0.0, 1.0)
cb_c = fig.colorbar(mpl.cm.ScalarMappable(norm=c_norm, cmap=c_cmap), cax=cax_c, orientation='horizontal',
                    spacing='uniform', boundaries=c_bounds)
cb_c.set_ticks([(lo + hi) / 2 for lo, hi in zip(c_bounds[:-1], c_bounds[1:])])
cb_c.set_ticklabels([f'−{abs(c_min)} to −11', '−10 to −6', '−5 to −3', '−2 to −1', '0',
                     '1 to 2', '3 to 5', '6 to 10', f'11 to {c_max}'])
cb_c.ax.tick_params(length=0, labelsize=5.5)
cb_c.set_label('Net change in richness (species gained − species lost)', labelpad=3)
cb_c.outline.set_linewidth(0.4)

for ext in ('png', 'pdf'):
    fig.savefig(os.path.join(output_dir, f"Figure1.{ext}"), dpi=600)
plt.close(fig)
print(f"Saved {output_dir}/Figure1.png/.pdf")
