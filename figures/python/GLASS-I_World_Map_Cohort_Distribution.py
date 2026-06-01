# GLASS-I: World Map Cohort Distribution
# Author: C.M.S. Tesileanu
# Date: 2026-05-27

import matplotlib.pyplot as plt
import geopandas as gpd
from shapely.geometry import Point
import pandas as pd
import numpy as np
from matplotlib.patches import Wedge, Patch, Circle

# Load world map and remove Antarctica
url = "https://naciscdn.org/naturalearth/110m/cultural/ne_110m_admin_0_countries.zip"
world = gpd.read_file(url)
world = world[world['NAME'] != 'Antarctica']

# Patient data
data = {
    'City': [
        'Cardiff', 'Chicago', 'Cleveland', 'Detroit', 'Durham',
        'Farmington', 'Gainesville', 'Heidelberg',  'Hong Kong','Houston', 
        'Luxembourg','Melbourne', 'New York', 'Paris', 'Phoenix', 
        'Rochester', 'San Francisco', 'Sao Paolo', 'Seoul', 'Tokyo'
    ],
    'Astrocytoma': [1, 0, 3, 7, 6,
                    1, 0, 0, 2, 37,
                    5, 0, 2, 2, 8, 
                    1, 36, 3, 6, 4],
    'Oligodendroglioma': [0, 2, 1, 5, 23, 
                          0, 1, 1, 1, 17, 
                          3, 6, 5, 0, 10, 
                          12, 22, 0, 1, 5],
    'Latitude': [
        51.4816, 41.8781, 41.4993, 42.3314, 35.9940, 
        41.7360, 29.6516, 49.3988, 22.3193,29.7604, 
        49.6117, -37.8136, 40.7128, 48.8566, 33.4484,
        44.0121, 37.7749, -23.5505, 37.5665, 35.6895
    ],
    'Longitude': [
        -3.1791, -87.6298, -81.6944, -83.0458, -78.8986, 
        -72.6851, -82.3248, 8.6724, 114.1694,-95.3698, 
        6.1319, 144.9631, -74.0060, 2.3522, -112.0740,
        -92.4802, -122.4194, -46.6333, 126.9780, 139.6917
    ]
}

# DataFrame and GeoDataFrame
df = pd.DataFrame(data)
geometry = [Point(xy) for xy in zip(df['Longitude'], df['Latitude'])]
gdf_patients = gpd.GeoDataFrame(df, geometry=geometry)

# Plot setup
fig, ax = plt.subplots(figsize=(6, 12))
world.plot(ax=ax, facecolor='white', edgecolor='gray', linewidth=2)
ax.set_aspect('equal')

# Colors and transparency
colors = ['#800074', '#298C8C']
alpha = 0.6

# Draw pie charts
for idx, row in df.iterrows():
    x, y = row['Longitude'], row['Latitude']
    raw_values = np.array([row['Astrocytoma'], row['Oligodendroglioma']])
    log_values = np.log1p(raw_values)
    total_log = log_values.sum()
    if total_log == 0:
        continue

    fracs = log_values / total_log
    radius = 0.2 + total_log * 2  # FIXED: correctly compute pie size from total_log
    theta1 = 0
    for i, frac in enumerate(fracs):
        theta2 = theta1 + frac * 360
        wedge = Wedge(center=(x, y), r=radius, theta1=theta1, theta2=theta2,
                      facecolor=colors[i], edgecolor='black', linewidth=1.5, alpha=alpha)
        ax.add_patch(wedge)
        theta1 = theta2

# Tumor type legend
tumor_legend = [
    Patch(facecolor=colors[0], edgecolor='black', linewidth=1.5, label='IDH-mutant Astrocytoma', alpha=alpha),
    Patch(facecolor=colors[1], edgecolor='black', linewidth=1.5, label='IDH-mutant Oligodendroglioma', alpha=alpha)
]
leg1 = ax.legend(handles=tumor_legend, loc='lower center', frameon=False)
ax.add_artist(leg1)

# Pie size scale legend
size_vals = [1, 5, 15, 50]
legend_x, legend_y = -150, -45  # Lower-left corner of map for better padding

for i, val in enumerate(size_vals):
    log_val = np.log1p(val)
    radius = 0.2 + log_val * 2
    circ = Circle((legend_x, legend_y + i * 16), radius=radius,
                  facecolor='white', edgecolor='black', linewidth=1.5)
    ax.add_patch(circ)
    ax.text(legend_x + 10, legend_y + i * 16, f'{val} patient{"s" if val > 1 else ""}', 
            va='center', fontsize=11)

# Final touches
ax.set_xlim(-169, 191)
ax.set_ylim(-90, 90)
ax.axis('off')

plt.tight_layout()
plt.savefig('~/idh_tumor_map.pdf')
plt.show()
