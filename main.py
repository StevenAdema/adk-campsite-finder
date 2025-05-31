import sys
from typing import List, Tuple
from dataclasses import dataclass
import rasterio
import numpy as np
from pystac_client import Client
import planetary_computer as pc
from pyproj import Transformer, Geod
from rasterio.windows import from_bounds
import folium
from typing import Optional
import math
from folium.plugins import MarkerCluster
from html_utils import create_map, add_layers_to_map, add_opacity_slider, add_markers_to_map, add_layer_control, add_click_to_copy, add_bounding_box_to_map
import os
import geopandas as gpd
from shapely.geometry import Point

# ─── DATACLASS ──────────────────────────────────────────────────────────────
@dataclass
class Coord:
    latitude:  float
    longitude: float
    ndvi:      float
    elevation: float = None
    slope:     float = None
    d_from_water: float = None
    d_from_trail: float = None

# ─── CONFIG ──────────────────────────────────────────────────────────────────
STAC_URL     = "https://planetarycomputer.microsoft.com/api/stac/v1"
CENTER_LON   = -73.75097
CENTER_LAT   = 44.13082
RADIUS_KM    = 0.5       # half width of square, in km
SAMPLE_N     = 100  # approx total points in the square
ELEV_RAD     = 10      # distance in metres to sample elevation and slope around each Coord

# ─── THRESHOLDS FOR FILTERING ──────────────────────────────────────────
NDVI_THRESH         = 0.25  # vegetation threshold 0 to 1 where 0 is no vegetation and 1 is dense vegetation
MIN_ELEV_THRESH     = 800   # min elevation in meters
MAX_ELEV_THRESH     = 1220  # max elevation in meters
SLOPE_THRESH        = 8     # slope threshold in degrees
DISTANCE_FROM_WATER = 5    # distance in metres from water bodies
DISTANCE_FROM_TRAIL = 45    # distance in metres from trails

# establish connection to Planetary Computer STAC
catalog = Client.open(STAC_URL)

def load_trail_layer(trail_file: str) -> gpd.GeoDataFrame:
    """
    Load a trail layer from a shapefile.
    
    Args:
        trail_file (str): Path to the trail shapefile.
    
    Returns:
        gpd.GeoDataFrame: GeoDataFrame containing the trail data.
    """
    try:
        trail_gdf = gpd.read_file(trail_file)
        print(f"Loaded trail layer with {len(trail_gdf)} features.")
        return trail_gdf
    except Exception as e:
        print(f"Error loading trail layer: {e}")
        return gpd.GeoDataFrame()  # Return an empty GeoDataFrame if loading fails

def create_coords_grid(
    center_lon: float, center_lat: float,
    radius_km: float, sample_n: int
) -> List[Coord]:
    """Create a grid of Coord (lat/lon, ndvi=None) in the target area."""
    item = pc.sign(next(catalog.search(
        collections=["naip"],
        intersects={"type":"Point","coordinates":[center_lon, center_lat]},
        sortby=[{"field":"properties.datetime","direction":"desc"}],
        limit=1
    ).items()))
    href = item.assets["image"].href

    with rasterio.open(href) as src:
        transformer_to_src = Transformer.from_crs("EPSG:4326", src.crs, always_xy=True)
        cx, cy = transformer_to_src.transform(center_lon, center_lat)
        half_side = radius_km * 1000.0

        bounds = (cx - half_side, cy - half_side, cx + half_side, cy + half_side)
        window = from_bounds(*bounds, transform=src.transform)
        win_transform = src.window_transform(window)

        h = int(window.height)
        w = int(window.width)
        grid_size = int(np.ceil(np.sqrt(sample_n)))
        rows = np.linspace(0, h-1, grid_size, dtype=int)
        cols = np.linspace(0, w-1, grid_size, dtype=int)
        rows_flat, cols_flat = np.meshgrid(rows, cols)
        rows_flat = rows_flat.ravel()
        cols_flat = cols_flat.ravel()

        xs, ys = rasterio.transform.xy(win_transform, rows_flat, cols_flat, offset="center")
        xs = np.array(xs); ys = np.array(ys)

        transformer_to_wgs = Transformer.from_crs(src.crs, "EPSG:4326", always_xy=True)
        lons, lats = transformer_to_wgs.transform(xs, ys)

    coords_list = [
        Coord(latitude=lat, longitude=lon, ndvi=None)
        for lat, lon in zip(lats, lons)
    ]
    return coords_list


def set_ndvi(coords: List[Coord]) -> None:
    """Set NDVI for each Coord object in the list, loading NAIP images covering all coords."""
    if not coords:
        return

    # Compute bounding box for all coords
    min_lon = min(c.longitude for c in coords)
    max_lon = max(c.longitude for c in coords)
    min_lat = min(c.latitude for c in coords)
    max_lat = max(c.latitude for c in coords)

    try:
        # Get ALL NAIP images that intersect with the bounding box
        search_results = catalog.search(
            collections=["naip"],
            intersects={
                "type": "Polygon",
                "coordinates": [[
                    [min_lon, min_lat],
                    [min_lon, max_lat],
                    [max_lon, max_lat],
                    [max_lon, min_lat],
                    [min_lon, min_lat]
                ]]
            },
            sortby=[{"field": "properties.datetime", "direction": "desc"}]
        )
        
        # Get all items, not just the first one
        items = list(search_results.items())
        if not items:
            print(f"No NAIP imagery found for bounding box")
            for c in coords:
                c.ndvi = 1
            return
        
        print(f"Found {len(items)} NAIP images covering the area")
        
        # Process each image and its overlapping coordinates
        processed_coords = set()
        
        for item in items:
            try:
                signed_item = pc.sign(item)
                href = signed_item.assets["image"].href
                
                with rasterio.open(href) as src:
                    to_src = Transformer.from_crs("EPSG:4326", src.crs, always_xy=True)
                    
                    # Find coordinates that fall within this image's bounds
                    for c in coords:
                        if id(c) in processed_coords:
                            continue  # Already processed this coordinate
                            
                        x, y = to_src.transform(c.longitude, c.latitude)
                        
                        # Check if point is within this raster's bounds
                        if (src.bounds.left <= x <= src.bounds.right and 
                            src.bounds.bottom <= y <= src.bounds.top):
                            
                            try:
                                vals = next(src.sample([(x, y)]))
                                # NAIP: Band 1=Red, Band 4=NIR
                                red = float(vals[0])
                                nir = float(vals[3])
                                
                                if red == 0 or nir == 0:
                                    c.ndvi = 1
                                else:
                                    denom = nir + red
                                    if denom == 0:
                                        c.ndvi = 0
                                    else:
                                        c.ndvi = (nir - red) / denom
                                
                                processed_coords.add(id(c))
                                
                            except Exception as e:
                                print(f"Error sampling coordinate ({c.latitude}, {c.longitude}): {e}")
                                c.ndvi = 1
                                processed_coords.add(id(c))
                                
            except Exception as e:
                print(f"Error processing NAIP image: {e}")
                continue
        
        # Set remaining unprocessed coordinates to 1
        for c in coords:
            if id(c) not in processed_coords:
                print(f"No coverage found for coordinate ({c.latitude}, {c.longitude})")
                c.ndvi = 1
                
    except Exception as e:
        print(f"Error in NDVI processing: {e}")
        for c in coords:
            c.ndvi = 1


def set_ndvi_from_local_tifs(coords: List[Coord], tif_directory: str) -> None:
    """Set NDVI for each Coord object using pre-downloaded TIF files."""
    if not coords:
        return

    # Compute bounding box for all coords
    min_lon = min(c.longitude for c in coords)
    max_lon = max(c.longitude for c in coords)
    min_lat = min(c.latitude for c in coords)
    max_lat = max(c.latitude for c in coords)

    # Find TIF files that intersect with the bounding box
    tif_files = [f for f in os.listdir(tif_directory) if f.endswith(".tif")]
    if not tif_files:
        print(f"No TIF files found in directory: {tif_directory}")
        for c in coords:
            c.ndvi = 1
        return

    print(f"Found {len(tif_files)} TIF files in {tif_directory}")

    processed_coords = set()

    for tif_file in tif_files:
        try:
            tif_path = os.path.join(tif_directory, tif_file)
            with rasterio.open(tif_path) as src:
                to_src = Transformer.from_crs("EPSG:4326", src.crs, always_xy=True)

                # Find coordinates that fall within this image's bounds
                for c in coords:
                    if id(c) in processed_coords:
                        continue  # Already processed this coordinate

                    x, y = to_src.transform(c.longitude, c.latitude)

                    # Check if point is within this raster's bounds
                    if (src.bounds.left <= x <= src.bounds.right and
                        src.bounds.bottom <= y <= src.bounds.top):

                        try:
                            vals = next(src.sample([(x, y)]))
                            # NAIP: Band 1=Red, Band 4=NIR
                            red = float(vals[0])
                            nir = float(vals[3])

                            if red == 0 or nir == 0:
                                c.ndvi = 1
                            else:
                                denom = nir + red
                                if denom == 0:
                                    c.ndvi = 0
                                else:
                                    c.ndvi = (nir - red) / denom

                            processed_coords.add(id(c))

                        except Exception as e:
                            print(f"Error sampling coordinate ({c.latitude}, {c.longitude}): {e}")
                            c.ndvi = 1
                            processed_coords.add(id(c))

        except Exception as e:
            print(f"Error processing TIF file {tif_file}: {e}")
            continue

    # Set remaining unprocessed coordinates to 1
    for c in coords:
        if id(c) not in processed_coords:
            print(f"No coverage found for coordinate ({c.latitude}, {c.longitude})")
            c.ndvi = 1


def set_elevation_and_slope(coords: List[Coord], distance = ELEV_RAD) -> None:
    """
    Batch-populate elevation & slope on each Coord using a 3DEP seamless DEM.
    - Uses degree-offset approximations instead of per-point Geod calls.
    - Samples all points in one src.sample() call.
    """
    # fetch the latest 3DEP‐Seamless DEM covering our coords
    catalog = Client.open(STAC_URL)
    item = pc.sign(
        next(
            catalog.search(
                collections=["3dep-seamless"],
                intersects={
                    "type": "Point",
                    "coordinates": [coords[0].longitude, coords[0].latitude]
                },
                sortby=[{"field": "properties.datetime", "direction": "desc"}],
                limit=1
            ).items()
        )
    )
    href = item.assets["data"].href

    # open DEM once and prepare batch sampling
    with rasterio.open(href) as src:
        to_src = Transformer.from_crs("EPSG:4326", src.crs, always_xy=True)
        n = len(coords)

        # build flattened lists of [N, E, S, W, C] for every coord
        lats, lons = [], []
        for c in coords:
            # approximate metre→degree offsets
            dlat = distance / 111320.0
            dlon = distance / (111320.0 * math.cos(math.radians(c.latitude)))
            # append N, E, S, W, and center
            lats.extend([c.latitude + dlat,
                         c.latitude,
                         c.latitude - dlat,
                         c.latitude,
                         c.latitude])
            lons.extend([c.longitude,
                         c.longitude + dlon,
                         c.longitude,
                         c.longitude - dlon,
                         c.longitude])

        # transform all to raster CRS and sample once
        xs, ys = to_src.transform(lons, lats)
        samples = list(src.sample(zip(xs, ys)))
        elevs = np.array([s[0] for s in samples]).reshape(n, 5)

        # assign elevation & compute slope per coord
        for idx, c in enumerate(coords):
            elev_n, elev_e, elev_s, elev_w, elev_c = elevs[idx]
            c.elevation = float(elev_c)
            dz_dx = (elev_e - elev_w) / (2 * distance)
            dz_dy = (elev_n - elev_s) / (2 * distance)
            c.slope = math.degrees(math.atan(math.hypot(dz_dx, dz_dy)))


def set_distance_from_trail_vector(coords: List[Coord], trail_gdf: gpd.GeoDataFrame) -> None:
    """
    Calculate the distance from each Coord to the nearest trail using a vector layer.
    Updates the Coord.d_from_trail attribute in meters.
    
    Args:
        coords (List[Coord]): List of Coord objects.
        trail_gdf (gpd.GeoDataFrame): GeoDataFrame containing trail geometries.
    """
    if trail_gdf.empty:
        print("Trail layer is empty. Cannot calculate distances.")
        for c in coords:
            c.d_from_trail = float('inf')
        return

    # Combine all trail geometries into a single MultiLineString
    trails_union = trail_gdf.geometry.unary_union

    # Initialize geodesic calculator
    geod = Geod(ellps="WGS84")

    for c in coords:
        coord_point = Point(c.longitude, c.latitude)
        # Calculate the distance using shapely's distance method
        nearest_point = trails_union.interpolate(trails_union.project(coord_point))
        _, _, distance_meters = geod.inv(c.longitude, c.latitude, nearest_point.x, nearest_point.y)
        c.d_from_trail = abs(distance_meters)


if __name__ == "__main__":
    print("Generating Coord grid...")
    coords = create_coords_grid(CENTER_LON, CENTER_LAT, RADIUS_KM, SAMPLE_N)
    print(f"Generated {len(coords)} Coord objects")

    print("Loading trail layer...")
    trail_file = os.path.join(os.path.dirname(__file__), "data/gis_osm_roads_free_1.shp")
    trail_gdf = load_trail_layer(trail_file)

    print("Calculating distance to nearest trail for each Coord...")
    set_distance_from_trail_vector(coords, trail_gdf)

    # Print all coords and their distance from trail
    for c in coords:
        print(f"Coord ({c.latitude}, {c.longitude}) - Distance from trail: {c.d_from_trail:.2f} m")

    print("Getting NDVI for each Coord...")
    set_ndvi_from_local_tifs(coords, os.path.join(os.path.dirname(__file__), "data"))

    print("Getting elevation and slope for each Coord...")
    set_elevation_and_slope(coords)

    campsites = [c for c in coords if c.d_from_trail <= DISTANCE_FROM_TRAIL]
    print(f"Found {len(campsites)} potential campsites")

    # Create map
    m = create_map(CENTER_LAT, CENTER_LON)

    # Add layers
    add_layers_to_map(m)

    # Add bounding box
    add_bounding_box_to_map(m, CENTER_LAT, CENTER_LON, RADIUS_KM)

    # Add opacity slider
    add_opacity_slider(m)

    # Add markers
    add_markers_to_map(m, campsites)

    # Add layer control
    add_layer_control(m)

    # Add click-to-copy functionality
    add_click_to_copy(m)

    # Save map
    m.save("ADK_Campsites.html")
    print("Map saved to ADK_Campsites.html")
