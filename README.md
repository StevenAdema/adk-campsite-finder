# 🏕️ ADK Campsite Finder

Find wilderness campsites in the Adirondacks Park (ADK) the meet DEC requiements as well as custom settings.

---
## 🥾 Background

My wife and I love hiking in the Adirondacks, but we often find ourselves seeking quieter, more secluded campsites away from the crowds. While the designated sites along the trails are convenient, they can fill up quickly, especially during weekends near popular peaks. The goal was to identify potential backcountry camping locations that meet both DEC wilderness camping regulations as well as our own personal preferences: higher elevation, near water source, not in dense forest, etc.

---

## ⚙️ How It Works

The ADK Campsite Finder uses geospatial tools to locate potential campsites in the Adirondacks. It creates a grid containing thousands of coordinates, all of which will be evaluated. The grid is analyzed using data from Microsoft's Planetary Computer catalog of high-resolution imagery, elevation data, and topographic features.

Coordinates that meet critera below are plotted on a map that can be loaded in any browser.

##### DEC Requirements
1. No camping on summits
2. No camping above 3,500 feet
3. No camping within 150 feet of any road, trail, or water feature
4. Tent must not be placed on vegetation.

##### Personal Requirements
5. Ground is relatively flat.
6. Camp above 2,000 feet elevation.
7. Camp within 0.5 miles of trail (minimize bushwacking)


### **Parameters**
`NDVI_THRESH`<br>
values: 0 - 1 (0: no vegetation, 1: dense vegetation)<br>
default: 0.2<br>
from the code:<br>
NDVI (Normalized Difference Vegetation Index) measures vegetation density using Near-Infrared (NIR) and Red light reflectance. It is calculated as (NIR - Red) / (NIR + Red).

NDVI is calculated for each coordinate using NAIP imagery or local TIF files. The set_ndvi function filters coordinates below a threshold (NDVI_THRESH), identifying open areas suitable for camping while avoiding dense forests or wetlands. This ensures compliance with DEC rules and personal preferences.

`MIN_ELEV_THRESH`, `MAX_ELEV_THRESH`  
values: Elevation in meters  
default: 800 (min), 1220 (max)  
from the code:  
Elevation is calculated using the 3DEP seamless dataset from Microsoft's Planetary Computer. For each coordinate, elevation is sampled at the center point and compared against the thresholds (`MIN_ELEV_THRESH` and `MAX_ELEV_THRESH`). Coordinates with elevation values outside this range are excluded. This ensures campsites are at comfortable elevations and comply with DEC rules.


`SLOPE_THRESH`  
values: Degrees (0 - 90)  
default: 8  
from the code:  
Slope is calculated using elevation values sampled at five points around each coordinate: north, south, east, west, and center. The slope is derived using the formula:<br>
dz_dx = (elev_e - elev_w) / (2 * distance)<br>
dz_dy = (elev_n - elev_s) / (2 * distance)<br>
slope = atan(hypot(dz_dx, dz_dy))

`DISTANCE_FROM_WATER`<br>
Work in progress.

`DISTANCE_FROM_TRAIL`<br>
Work in progress.

---

## 📸 Demo

![Map Screenshot](/img/screenshot.png)

- Red bounding box: Defines the area of interest.
- Tent Icon: Potential campsites.
- Red markers: Locations that do not meet the thresholds (if enabled).

---

## 🛠️ Installation

### **Prerequisites**
- Python 3.8+
- Pip
- Virtual environment (optional)

### **Steps**
1. Clone the repository:
   ```bash
   git clone https://github.com/stevenadema/adk-campsite-finder.git
   cd adk-campsite-finder
2. Set thresholds in main.py
3. Run main.py
```bash
   python main.py 
```
4. Open map ADK_Campsites.html in browser