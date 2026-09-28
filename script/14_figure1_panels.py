# Figure 1 panels (Robinson maps).

import os
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
from matplotlib import cm
import seaborn as sns
import cartopy.crs as ccrs
from tqdm import tqdm

warnings.filterwarnings('ignore')

# Journal style
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial'],
    'font.size': 18,
    'axes.labelsize': 20,
    'axes.titlesize': 22,
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
    'savefig.facecolor': 'white',
    'figure.facecolor': 'white'
})

# Paths (raw maps not shared)
shp_path = "data/world-administrative-boundaries/world-administrative-boundaries.shp"
csv_path = "output/df_lat_sp.csv"
data_dir = "data/SuitabilityMaps_MigratorySpecies"
output_dir = "output/Manuscript/figure1_panels"
os.makedirs(output_dir, exist_ok=True)

# Load boundaries and rasters
world = gpd.read_file(shp_path).dissolve()

df_ref = pd.read_csv(csv_path)
sp_to_fam = pd.Series(df_ref.Family.values, index=df_ref.species).to_dict()
all_files = glob.glob(os.path.join(data_dir, "Binary_S*.tif"))

species_dict = {}
for f in all_files:
    m = re.match(r"Binary_S(\d)(.*)\.tif$", os.path.basename(f))
    if m:
        s, sp = int(m.group(1)), m.group(2)
        if sp in sp_to_fam:
            if sp not in species_dict:
                species_dict[sp] = {}
            species_dict[sp][s] = f

# Native-resolution master grid
with rasterio.open(all_files[0]) as src:
    height, width = int(180 / src.res[1]), int(360 / src.res[0])
    master_transform = from_bounds(-180, -90, 180, 90, width, height)

migratory_richness_map = np.zeros((height, width), dtype='float32')
S = {s: np.zeros((height, width), dtype='float32') for s in [1, 2, 3, 4]}

for sp, seasons in tqdm(species_dict.items(), desc="Species"):
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

max_richness = int(np.max(migratory_richness_map))
print(f"Maximum richness: {max_richness} species")

P_STEP = 2  # Render every 2nd cell
proj_robinson = ccrs.Robinson(central_longitude=0)
plate_carree = ccrs.PlateCarree()

# Panel a: Robinson map
fig_a, ax_a = plt.subplots(figsize=(16, 9), subplot_kw={'projection': proj_robinson})
ax_a.set_global()
ax_a.spines['geo'].set_linewidth(1.2)
ax_a.spines['geo'].set_edgecolor('black')
ax_a.add_geometries(world.geometry, crs=plate_carree,
                    facecolor='white', edgecolor='black', linewidth=0.3, alpha=0.3, zorder=1)
masked_a = np.ma.masked_equal(migratory_richness_map, 0)
ax_a.imshow(masked_a[::P_STEP, ::P_STEP], origin='upper',
            extent=[-180, 180, -90, 90], transform=plate_carree,
            cmap='viridis', zorder=2)
gl_a = ax_a.gridlines(crs=plate_carree, draw_labels=False, linewidth=1, color='black', linestyle='--', alpha=0.6)
gl_a.ylocs = [0]  # Equator only
gl_a.xlocs = []
plt.savefig(os.path.join(output_dir, "Panel_A_Map_Robinson_Clean.png"), dpi=300, bbox_inches='tight')
plt.close()

# Panel a colour bar
fig_cb_a, ax_cb_a = plt.subplots(figsize=(8, 1.8))
norm_a = mpl.colors.Normalize(vmin=1, vmax=max_richness)
cb_a = fig_cb_a.colorbar(mpl.cm.ScalarMappable(norm=norm_a, cmap='viridis'), cax=ax_cb_a,
                         orientation='horizontal', ticks=[1, 10, 20, 30, 40, 50, max_richness])
cb_a.ax.xaxis.set_ticks_position('top')
cb_a.set_label('Seasonal-switching richness', fontsize=18, labelpad=12)
cb_a.ax.xaxis.set_label_position('top')
cb_a.ax.tick_params(labelsize=15, length=6, width=1.2)
plt.savefig(os.path.join(output_dir, "Fig1_Panel_A_Colorbar_Standalone.png"), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(output_dir, "Fig1_Panel_A_Colorbar_Standalone.pdf"), bbox_inches='tight', transparent=True)
plt.close()

# Panel b: latitude profile
lats = np.linspace(90, -90, height)
lat_mean = np.zeros(height)
for i in range(height):
    row_vals = migratory_richness_map[i, :]
    occupied = row_vals[row_vals > 0]
    lat_mean[i] = np.mean(occupied) if len(occupied) > 0 else 0.0

fig_b, ax_b = plt.subplots(figsize=(8, 10))
ax_b.plot(lat_mean, lats, color=cm.viridis(0.3), linewidth=3.5)
ax_b.fill_betweenx(lats, 0, lat_mean, color=cm.viridis(0.6), alpha=0.25)
ax_b.set_ylim(-60, 90)
ax_b.set_yticks([-60, -30, 0, 30, 60, 90])
ax_b.set_yticklabels(['60°S', '30°S', '0°', '30°N', '60°N', '90°N'], fontsize=18)
ax_b.set_xticks(np.arange(0, int(np.max(lat_mean)) + 5, 5))
ax_b.tick_params(axis='x', labelsize=16)
ax_b.set_xlabel("Mean richness (occupied cells)", fontsize=18, labelpad=10)
ax_b.set_ylabel("Latitude", fontsize=18, labelpad=10)
ax_b.axhline(0, color='black', linestyle='--', linewidth=1.2, alpha=0.6)
sns.despine(ax=ax_b)
plt.savefig(os.path.join(output_dir, "Fig1_Panel_B_Profile.png"), dpi=300, bbox_inches='tight')
plt.close()

# Panel c: seasonal net change
fig_c, axes_c = plt.subplots(2, 2, figsize=(18, 12),
                             subplot_kw={'projection': proj_robinson},
                             gridspec_kw={'hspace': 0.12, 'wspace': 0.05})
labels = [("S2–S1", S[2] - S[1]), ("S3–S2", S[3] - S[2]), ("S4–S3", S[4] - S[3]), ("S1–S4", S[1] - S[4])]
v_lim = max(np.abs(d).max() for _, d in labels)  # Full range, no clipping
for i, (label, data) in enumerate(labels):
    ax = axes_c.flatten()[i]
    ax.set_global()
    ax.spines['geo'].set_linewidth(1.0)
    ax.spines['geo'].set_edgecolor('black')
    ax.add_geometries(world.geometry, crs=plate_carree,
                      facecolor='white', edgecolor='black', linewidth=0.2, alpha=0.3, zorder=1)
    data_m = np.ma.masked_where(np.abs(data) < 0.1, data)
    ax.imshow(data_m[::P_STEP, ::P_STEP], origin='upper',
              extent=[-180, 180, -90, 90], transform=plate_carree,
              cmap='viridis', vmin=-v_lim, vmax=v_lim, zorder=2)
    gl_c = ax.gridlines(crs=plate_carree, draw_labels=False, linewidth=0.8, color='black', linestyle='--', alpha=0.5)
    gl_c.ylocs = [0]  # Equator only
    gl_c.xlocs = []
    ax.set_title(label, fontsize=20, pad=8)
plt.savefig(os.path.join(output_dir, "Panel_C_Seasons_Robinson_Clean.png"), dpi=300, bbox_inches='tight')
plt.close()
print(f"Panel c colour limit: +/-{v_lim:.1f}")


# Panel c colour bars
def export_c_colorbar(v_range, ticks, suffix):
    fig_cb_c, ax_cb_c = plt.subplots(figsize=(10, 1.8))
    norm_c = mpl.colors.Normalize(vmin=-v_range, vmax=v_range)
    cb_c = fig_cb_c.colorbar(mpl.cm.ScalarMappable(norm=norm_c, cmap='viridis'), cax=ax_cb_c,
                             orientation='horizontal', ticks=ticks)
    cb_c.set_label('Net species flux (loss ← → gain)', fontsize=20, labelpad=12)
    cb_c.ax.tick_params(labelsize=16, length=6, width=1.2)
    plt.savefig(os.path.join(output_dir, f"Fig1_Panel_C_Colorbar_{suffix}.png"), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, f"Fig1_Panel_C_Colorbar_{suffix}.pdf"), bbox_inches='tight', transparent=True)
    plt.close()


v_int = int(v_lim)
export_c_colorbar(v_int, [-v_int, -20, -10, 0, 10, 20, v_int], f"pm{v_int}")
print(f"Panels saved to {output_dir}; assembled manually.")
