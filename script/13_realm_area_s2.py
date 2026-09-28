# Equal-area realm metrics (Fig. S2).

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
from rasterio.mask import mask
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm

warnings.filterwarnings('ignore')

# Paths (raw maps not shared)
DATA_DIR = "data/SuitabilityMaps_MigratorySpecies"
CSV_PATH = "output/df_lat_sp.csv"
ZOO_SHP = "data/CMEC_regions/Regions.shp"
OUTPUT_DIR = "output/realm_area"

os.makedirs(OUTPUT_DIR, exist_ok=True)
GRID_RES = 0.5

# Regions to realms (Holt 2013)
REGION_TO_REALM = {
    'North American': 'Nearctic',
    'Mexican': 'Nearctic',
    'Panamanian': 'Panamanian',
    'South American': 'Neotropical',
    'Amazonian': 'Neotropical',
    'Eurasian': 'Palaearctic',
    'Arctico-Siberian': 'Palaearctic',
    'Tibetan': 'Palaearctic',
    'Chinese': 'Sino-Japanese',
    'Japanese': 'Sino-Japanese',
    'Oriental': 'Oriental',
    'Indo-Malayan': 'Oriental',
    'African': 'Afrotropical',
    'Guineo-Congolian': 'Afrotropical',
    'Madagascan': 'Madagascan',
    'Saharo-Arabian': 'Saharo-Arabian',
    'Australian': 'Australian',
    'Novozelandic': 'Australian',
    'Papua-Melanesian': 'Australian',
    'Polynesian': 'Oceanian'
}

# Richness grid (0.5 degree)
print("[1/5] Seasonal-switching richness...")
df_ref = pd.read_csv(CSV_PATH)
sp_to_fam = pd.Series(df_ref.Family.values, index=df_ref.species).to_dict()
all_files = glob.glob(os.path.join(DATA_DIR, "Binary_S*.tif"))

species_dict = {}
for f in all_files:
    m = re.match(r"Binary_S(\d)(.*)\.tif$", os.path.basename(f))
    if m:
        s, sp = int(m.group(1)), m.group(2)
        if sp in sp_to_fam:
            if sp not in species_dict:
                species_dict[sp] = {}
            species_dict[sp][s] = f

with rasterio.open(all_files[0]) as src:
    master_crs = src.crs
    width = int(360 / GRID_RES)
    height = int(180 / GRID_RES)
    master_transform = from_bounds(-180, -90, 180, 90, width, height)

migratory_richness_map = np.zeros((height, width), dtype='float32')

for sp, seasons in tqdm(species_dict.items(), desc="Species"):
    t_max = len(seasons)
    temp = {s: np.zeros((height, width), dtype='float32') for s in [1, 2, 3, 4]}
    for s, path in seasons.items():
        with rasterio.open(path) as src:
            reproject(
                rasterio.band(src, 1),
                temp[s],
                src_transform=src.transform,
                src_crs=src.crs,
                dst_transform=master_transform,
                dst_crs=master_crs,
                resampling=Resampling.nearest
            )
    p_sum = sum(temp.values())
    mig_mask = np.where((p_sum > 0) & (p_sum < t_max), 1, 0).astype('float32')
    if np.sum(mig_mask) == 0:
        continue
    migratory_richness_map += mig_mask

# Spherical cell area (km2)
print("\n[2/5] Cell surface area...")
R = 6371.0088
d_lon_rad = np.radians(GRID_RES)
row_areas = np.zeros(height, dtype='float32')
for i in range(height):
    lat_top = 90.0 - i * GRID_RES
    lat_bot = 90.0 - (i + 1) * GRID_RES
    row_areas[i] = (R ** 2) * d_lon_rad * (np.sin(np.radians(lat_top)) - np.sin(np.radians(lat_bot)))

cell_area_km2 = np.repeat(row_areas[:, np.newaxis], width, axis=1).astype('float32')

# Global percentile thresholds
valid_occupied = migratory_richness_map[migratory_richness_map > 0]
top10_threshold = float(np.percentile(valid_occupied, 90))
top5_threshold = float(np.percentile(valid_occupied, 95))

# Temporary rasters for masking
temp_rich_tif = os.path.join(OUTPUT_DIR, "_temp_richness.tif")
temp_area_tif = os.path.join(OUTPUT_DIR, "_temp_area.tif")
profile = {
    'driver': 'GTiff', 'height': height, 'width': width, 'count': 1,
    'dtype': 'float32', 'crs': master_crs, 'transform': master_transform
}
with rasterio.open(temp_rich_tif, 'w', **profile) as dst:
    dst.write(migratory_richness_map, 1)
with rasterio.open(temp_area_tif, 'w', **profile) as dst:
    dst.write(cell_area_km2, 1)

# Extract per region
print("\n[3/5] Area and richness per region...")
gdf_regions = gpd.read_file(ZOO_SHP).to_crs(master_crs)

records_per_pixel = []

with rasterio.open(temp_rich_tif) as src_rich, rasterio.open(temp_area_tif) as src_area:
    for _, row in tqdm(gdf_regions.iterrows(), total=len(gdf_regions), desc="Regions"):
        region_name = str(row['Regions'])
        realm_name = REGION_TO_REALM.get(region_name, region_name)
        geom = [row['geometry']]

        try:
            r_img, _ = mask(src_rich, geom, crop=True)
            a_img, _ = mask(src_area, geom, crop=True)

            r_data = r_img[0]
            a_data = a_img[0]

            # Occupied cells (richness >= 1)
            mig_mask = (r_data > 0)
            if not np.any(mig_mask):
                continue

            r_valid = r_data[mig_mask].astype(int)
            a_valid = a_data[mig_mask]

            for rich, area in zip(r_valid, a_valid):
                records_per_pixel.append({
                    'Region': region_name,
                    'Realm': realm_name,
                    'Richness': rich,
                    'Area_km2': area
                })
        except Exception:
            continue

# Remove temporary rasters
for tmp in [temp_rich_tif, temp_area_tif]:
    if os.path.exists(tmp):
        os.remove(tmp)

df_all = pd.DataFrame(records_per_pixel)

# Summary tables
print("\n[4/5] Summary tables...")


def generate_summary(df, group_col):
    summary_list = []
    for name, group in df.groupby(group_col):
        tot_area = group['Area_km2'].sum()
        rich_vals = group['Richness'].values
        area_vals = group['Area_km2'].values

        area_ge_5 = group[group['Richness'] >= 5]['Area_km2'].sum()
        area_ge_10 = group[group['Richness'] >= 10]['Area_km2'].sum()
        area_ge_15 = group[group['Richness'] >= 15]['Area_km2'].sum()
        area_ge_20 = group[group['Richness'] >= 20]['Area_km2'].sum()
        area_top10 = group[group['Richness'] >= top10_threshold]['Area_km2'].sum()
        area_top5 = group[group['Richness'] >= top5_threshold]['Area_km2'].sum()

        summary_list.append({
            group_col: name,
            'Migratory_Area_km2': tot_area,
            'Mean_Richness': np.average(rich_vals, weights=area_vals),
            'Max_Richness': int(np.max(rich_vals)),
            # Fixed richness thresholds
            'Area_ge_5_km2': area_ge_5,
            'Pct_ge_5': (area_ge_5 / tot_area) * 100.0,
            'Area_ge_10_km2': area_ge_10,
            'Pct_ge_10': (area_ge_10 / tot_area) * 100.0,
            'Area_ge_15_km2': area_ge_15,
            'Pct_ge_15': (area_ge_15 / tot_area) * 100.0,
            'Area_ge_20_km2': area_ge_20,
            'Pct_ge_20': (area_ge_20 / tot_area) * 100.0,
            # Percentile thresholds
            'Area_Top10Pct_km2': area_top10,
            'Pct_Top10Pct': (area_top10 / tot_area) * 100.0,
            'Area_Top5Pct_km2': area_top5,
            'Pct_Top5Pct': (area_top5 / tot_area) * 100.0
        })

    df_res = pd.DataFrame(summary_list)
    # Share of global area
    tot_ge10_global = df_res['Area_ge_10_km2'].sum()
    tot_top5_global = df_res['Area_Top5Pct_km2'].sum()
    df_res['Global_Share_ge_10_Pct'] = (df_res['Area_ge_10_km2'] / tot_ge10_global) * 100.0
    df_res['Global_Share_Top5_Pct'] = (df_res['Area_Top5Pct_km2'] / tot_top5_global) * 100.0
    return df_res.sort_values('Migratory_Area_km2', ascending=False)


# Realm table
df_realm_summary = generate_summary(df_all, 'Realm')
csv_realm_summary = os.path.join(OUTPUT_DIR, "11_Realms_Baseline_Summary.csv")
df_realm_summary.to_csv(csv_realm_summary, index=False)

# Region table
df_region_summary = generate_summary(df_all, 'Region')
csv_region_summary = os.path.join(OUTPUT_DIR, "20_Regions_Baseline_Summary.csv")
df_region_summary.to_csv(csv_region_summary, index=False)

# Area by richness, per realm
df_pivot_realm = df_all.pivot_table(index='Realm', columns='Richness', values='Area_km2', aggfunc='sum').fillna(0.0)
df_pivot_realm['Total_Area_km2'] = df_pivot_realm.sum(axis=1)
df_pivot_realm = df_pivot_realm.sort_values('Total_Area_km2', ascending=False)
df_pivot_realm.to_csv(os.path.join(OUTPUT_DIR, "Realm_Richness_Area_km2_Distribution.csv"))

# Fig. S2 heatmap
print("\n[5/5] Fig. S2 heatmap...")
df_norm = df_pivot_realm.drop(columns=['Total_Area_km2'])
df_norm = df_norm.div(df_norm.sum(axis=1), axis=0) * 100.0

plt.figure(figsize=(16, 9))
sns.set_style("white")
plot_data = df_norm.iloc[:, :35].reindex(df_pivot_realm.index)

sns.heatmap(
    plot_data,
    cmap="rocket_r",
    annot=False,
    cbar_kws={'label': "Percentage of Realm's Migratory Area (%)", 'shrink': 0.8},
    mask=(plot_data == 0)
)
plt.title("Relative Area Composition of Migratory Richness across Zoogeographic Realms", fontsize=18, pad=20, fontweight='bold')
plt.xlabel("Seasonal-Switching Richness (Number of Migratory Species)", fontsize=13)
plt.ylabel("Zoogeographic Realms (Ranked by Total Area)", fontsize=13)
plt.tight_layout()
out_fig_s2 = os.path.join(OUTPUT_DIR, "Fig_S2_Area_Normalized_Heatmap.png")
plt.savefig(out_fig_s2, dpi=300, bbox_inches='tight')
plt.close()

# Report realm metrics
print(df_realm_summary[['Realm', 'Migratory_Area_km2', 'Max_Richness', 'Pct_ge_10', 'Global_Share_ge_10_Pct']]
      .to_string(index=False))
mask_np = df_realm_summary['Realm'].str.contains('Nearctic|Palaearctic', case=False)
print(f"\nNearctic + Palaearctic share of area with >= 10 species: "
      f"{df_realm_summary[mask_np]['Global_Share_ge_10_Pct'].sum():.1f}%")
print(f"Saved: {csv_realm_summary}, {out_fig_s2}")
