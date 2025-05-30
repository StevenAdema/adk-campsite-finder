import folium
from folium import Element
from typing import List, Tuple
from dataclasses import dataclass
import math

@dataclass
class Coord:
    latitude: float
    longitude: float
    ndvi: float
    elevation: float = None
    slope: float = None

def create_map(center_lat: float, center_lon: float, zoom_start: int = 15) -> folium.Map:
    """Creates a Folium map centered at the given latitude and longitude with a custom tab title and logo."""
    m = folium.Map(location=[center_lat, center_lon], zoom_start=zoom_start, tiles=None)

    # Add custom tab title and logo
    title_script = Element("""
    <script>
    document.title = "⛺ ADK Campsite Finder";
    </script>
    """)
    m.get_root().html.add_child(title_script)

    return m

def add_layers_to_map(m: folium.Map) -> None:
    """Adds base layers and topographic layers to the map."""
    folium.TileLayer(
        tiles="https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}",
        attr="Esri World Imagery",
        name="Satellite"
    ).add_to(m)

    topo = folium.TileLayer(
        tiles="https://{s}.tile.opentopomap.org/{z}/{x}/{y}.png",
        attr="Map data: © OpenTopoMap (CC-BY-SA)",
        name="Topographic",
        overlay=True,
        control=True,
        show=True,
        opacity=1.0
    )
    topo.add_to(m)

def add_bounding_box_to_map(m: folium.Map, center_lat: float, center_lon: float, radius_km: float) -> None:
    """Adds a red bounding box to the map based on the center coordinates and radius."""
    # Calculate bounding box coordinates
    dlat = radius_km / 111.32  # Approximate degree offset for latitude
    dlon = radius_km / (111.32 * math.cos(math.radians(center_lat)))  # Approximate degree offset for longitude

    # Define the corners of the bounding box
    bounds = [
        [center_lat - dlat, center_lon - dlon],  # Bottom-left
        [center_lat - dlat, center_lon + dlon],  # Bottom-right
        [center_lat + dlat, center_lon + dlon],  # Top-right
        [center_lat + dlat, center_lon - dlon],  # Top-left
        [center_lat - dlat, center_lon - dlon],  # Closing the box
    ]

    # Add the bounding box to the map
    folium.PolyLine(
        bounds,
        color="red",
        weight=2,
        opacity=0.8,
        tooltip="Bounding Box"
    ).add_to(m)

    # Calculate the center of the top edge of the bounding box
    top_center_lat = center_lat + dlat
    top_center_lon = center_lon

def add_opacity_slider(m: folium.Map) -> None:
    """Adds an opacity slider for the topographic layer."""
    slider = Element("""
    <div id="opacity-slider" style="position: fixed; top: 80px; right: 20px; z-index:9999; background: white; padding: 8px; border-radius: 8px; box-shadow: 2px 2px 6px #888;">
      <label for="topo-opacity">Topo Opacity: <span id="opacity-value">100%</span></label><br>
      <input type="range" id="topo-opacity" min="0" max="1" value="1" step="0.01" style="width: 150px;">
    </div>
    <script>
    document.addEventListener('DOMContentLoaded', function() {
        setTimeout(function() {
            var slider = document.getElementById('topo-opacity');
            var valueDisplay = document.getElementById('opacity-value');
            var topoLayer = null;
            
            if (typeof window !== 'undefined') {
                for (var key in window) {
                    if (window[key] && typeof window[key] === 'object' && window[key]._layers) {
                        var map = window[key];
                        for (var layerId in map._layers) {
                            var layer = map._layers[layerId];
                            if (layer._url && layer._url.includes('opentopomap.org')) {
                                topoLayer = layer;
                                break;
                            }
                        }
                        if (topoLayer) break;
                    }
                }
            }
            
            if (slider && topoLayer) {
                slider.addEventListener('input', function(e) {
                    var opacity = parseFloat(e.target.value);
                    topoLayer.setOpacity(opacity);
                    valueDisplay.textContent = Math.round(opacity * 100) + '%';
                });
                console.log('Opacity slider connected successfully');
            } else {
                console.log('Could not find topographic layer or slider element');
            }
        }, 2000);
    });
    </script>
    """)
    m.get_root().html.add_child(slider)

def add_markers_to_map(m: folium.Map, coords: List[Coord], radius=3, color="lime") -> None:
    """Adds markers Coords."""
    for c in coords:
        folium.Marker(
            location=(c.latitude, c.longitude),
            icon=folium.DivIcon(html=f"""
            <div style="font-size: 24px; text-align: center; line-height: 1;">🏕️</div>
            """),
            tooltip=f"{c.latitude:.6f}, {c.longitude:.6f} | ndvi: {c.ndvi:.3f}, slope: {c.slope:.0f}°, elev: {c.elevation:.0f}"
        ).add_to(m)

def add_layer_control(m: folium.Map) -> None:
    """Adds a layer control to the map."""
    folium.LayerControl().add_to(m)

def add_click_to_copy(m: folium.Map) -> None:
    """Adds functionality to display latitude and longitude on click and add a marker with a label."""
    click_script = Element("""
    <script>
    document.addEventListener('DOMContentLoaded', function() {
        var map = null;
        for (var key in window) {
            if (window[key] && typeof window[key] === 'object' && window[key]._layers) {
                map = window[key];
                break;
            }
        }
        if (map) {
            map.on('click', function(e) {
                var lat = e.latlng.lat.toFixed(6);
                var lon = e.latlng.lng.toFixed(6);
                var coords = lat + ', ' + lon;

                // Remove the previous marker if it exists
                if (tempMarker) {
                    map.removeLayer(tempMarker);
                }

                // Add a new marker with a label
                tempMarker = L.marker([lat, lon]).addTo(map);
                tempMarker.bindPopup('Coordinates: ' + coords).openPopup();

                // Copy to clipboard
                navigator.clipboard.writeText(coords).then(function() {
                    console.log('Copied to clipboard: ' + coords);
                }).catch(function(err) {
                    console.error('Could not copy to clipboard', err);
                });
            });
        }
    });
    </script>
    """)
    m.get_root().html.add_child(click_script)