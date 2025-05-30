import numpy as np
import rasterio
from pyproj import Geod
import folium
from pystac_client import Client

# ── PARAMETERS ────────────────────────────────────────────────────────────────
lat, lon           = 44.12851, -73.75022               # center point
distance           = 10                                  # metres
dem_path           = "USGS_1_n45w074_20190416.tif"          # downloaded DEM tile
output_html        = "slope_map.html"
# ── COMPUTE OFFSET POINTS ────────────────────────────────────────────────────
geod = Geod(ellps="WGS84")
points = {}
for name, az in [("N",   0), ("E",  90), ("S", 180), ("W", 270), ("C", None)]:
    if name == "C":
        points["C"] = (lat, lon)
    else:
        lon2, lat2, _ = geod.fwd(lon, lat, az, distance)
        points[name] = (lat2, lon2)

# ── SAMPLE ELEVATIONS ───────────────────────────────────────────────────────
with rasterio.open(dem_path) as src:
    elevs = {}
    for key, (lat_i, lon_i) in points.items():
        row, col = src.index(lon_i, lat_i)
        elevs[key] = src.read(1)[row, col]
        print(f"{key}: {elevs[key]:.1f} m at ({lat_i:.6f}, {lon_i:.6f})")

# ── CALCULATE SLOPE ──────────────────────────────────────────────────────────
dz_dx = (elevs["E"] - elevs["W"]) / (2 * distance)
dz_dy = (elevs["N"] - elevs["S"]) / (2 * distance)
slope_deg = np.degrees(np.arctan(np.hypot(dz_dx, dz_dy)))
print(f"Slope ≈ {slope_deg:.2f}°")

# ── PLOT ON FOLIUM MAP ───────────────────────────────────────────────────────
m = folium.Map(location=[lat, lon], zoom_start=17, tiles="Esri.WorldImagery")

for key, (lat_i, lon_i) in points.items():
    folium.CircleMarker(
        location=[lat_i, lon_i],
        radius=4,
        color="yellow",
        fill=True,
        fill_opacity=0.8,
        popup=f"{key}: {elevs[key]:.1f} m"
    ).add_to(m)

m.save(output_html)
print(f"Map saved to {output_html}")
