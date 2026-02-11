# Gpr-Toolbox

**Automated hyperbola detection and GIS integration for ArcGIS Pro**

[![ArcGIS Pro](https://img.shields.io/badge/ArcGIS%20Pro-2.9%2B-0079C1?logo=esri&logoColor=white)](https://www.esri.com/arcgis-pro)
[![Python](https://img.shields.io/badge/Python-3.9%2B-3776AB?logo=python&logoColor=white)](https://www.python.org)
[![License](https://img.shields.io/badge/License-MIT-green)](LICENSE)
[![Esri Marketplace](https://img.shields.io/badge/Esri-Marketplace-0079C1?logo=esri)](https://www.esri.com/en-us/arcgis-marketplace)

> A Python Toolbox (`.pyt`) for ArcGIS Pro that automatically detects buried utilities from GSSI `.DZT` ground-penetrating radar files — no manual radargram interpretation required.

---

## Features

- **Automated hyperbola detection** — SciPy-powered multi-scale edge detection, adaptive thresholding, and `curve_fit()` hyperbola fitting with R² quality filtering
- **GPS interpolation** — Parses `.DZG` files, linearly interpolates GPS across all traces, and computes cumulative distance arrays
- **Null geometry on GPS absence** — No dummy `(0,0)` coordinates; detections without GPS are stored with `NULL` geometry and flagged fields for later joining
- **Selectable output coordinate system** — Accepts WKT, `.prj` path, or `SpatialReference` object; defaults to the active map's CRS
- **4-panel radargram export** — Enhanced, seismic contrast, edge-enhanced, and annotated detection panels at 400 dpi
- **Direct ArcGIS Pro integration** — Auto-adds point and polyline shapefiles to the active map on completion
- **Graceful degradation** — Runs in basic amplitude-detection mode if SciPy is unavailable; skips PNG export if matplotlib is missing

---

## Requirements

| Package | Status | Notes |
|---|---|---|
| `arcpy` | **Required** | Included with ArcGIS Pro |
| `numpy` | **Required** | Included with ArcGIS Pro |
| `scipy` | Recommended | Full hyperbola detection; falls back to basic mode without it |
| `matplotlib` | Recommended | Radargram PNG export; skipped if absent |

**ArcGIS Pro 2.9+** is required (tested on 3.x). The toolbox uses `arcpy.mp`, `arcpy.da.InsertCursor`, and `GPCoordinateSystem` parameter types.

---

## Installation

### 1. Install Optional Dependencies

Open the **Python Command Prompt** bundled with ArcGIS Pro (search "Python Command Prompt" in the Start menu):

```bash
pip install scipy matplotlib
```

### 2. Add the Toolbox

In ArcGIS Pro's **Catalog Pane**:

```
Toolboxes → right-click → Add Toolbox → select gpr_toolbox_.pyt
```

Or connect to the folder containing `gpr_toolbox_.pyt` and expand it directly.

---

## Usage

### Tool Parameters

| Parameter | Type | Default | Notes |
|---|---|---|---|
| Input Folder | `DEFolder` | — | Folder containing `.DZT` (and optionally `.DZG`) files |
| Output Folder | `DEFolder` | — | Created automatically if it doesn't exist |
| Survey Name | `GPString` | `GPR_Survey` | Prefix for output shapefile names |
| Velocity (m/ns) | `GPDouble` | auto | Leave blank to use dielectric from DZT header; try `0.09` for urban |
| Min Confidence | `GPDouble` | `0.58` | Lower to `0.55` for deep/narrow pipes; raise to `0.65` to reduce noise |
| Auto-add to map | `GPBoolean` | `True` | Loads outputs into the active ArcGIS Pro map on completion |
| Output Coordinate System | `GPCoordinateSystem` | Active map CRS | Accepts WKT, `.prj` path, or SR object; falls back to WGS 84 |

### Input File Structure

```
survey_folder/
├── line01.DZT       ← required
├── line01.DZG       ← optional (GPS)
├── line02.DZT
├── line02.DZG
└── ...
```

`.DZG` files must have the **same stem name** as their corresponding `.DZT` file. The tool scans for both `.DZG` and `.dzg` extensions.

### Output Structure

```
output_folder/
├── {survey}_Utilities.shp    ← point features, one per detected utility
├── {survey}_Utilities.dbf
├── {survey}_Utilities.shx
├── {survey}_Path.shp         ← polyline per survey file (NULL geom if no GPS)
├── Radargrams/
│   ├── line01_V9_1B.png      ← 4-panel radargram (if matplotlib installed)
│   └── line02_V9_1B.png
└── ...
```

---

## Output Schema

### `{survey}_Utilities.shp` — POINT

| Field | Type | Description |
|---|---|---|
| `UtilID` | TEXT | Unique identifier: `{stem}_U{n}` |
| `Depth_m` | DOUBLE | Estimated burial depth in metres |
| `Diam_m` | DOUBLE | Estimated diameter: `max(0.03, min(0.8, depth × 0.45))` m |
| `Confidence` | DOUBLE | Detection confidence 0–1 |
| `R_Squared` | DOUBLE | Hyperbola curve-fit R² |
| `Velocity` | DOUBLE | EM wave velocity used (m/ns) |
| `SourceFile` | TEXT | Source `.DZT` filename |
| `GPS_Avail` | TEXT | `'Yes'` if GPS data was parsed |
| `XY_Valid` | TEXT | `'Yes'` if geometry contains real coordinates |

### `{survey}_Path.shp` — POLYLINE

| Field | Type | Description |
|---|---|---|
| `FileName` | TEXT | Source `.DZT` filename |
| `GPS_Avail` | TEXT | `'Yes'` / `'No'` |

---

## Detection Algorithm

The `HyperbolaDetector` class processes each file through six stages:

1. **Preprocess** — Median filter (kernel=3) → DC removal (trace-mean subtraction) → adaptive gain ramp (1×→5× past sample 100) → ±98th percentile clip
2. **Edge Enhancement** — Multi-scale Sobel gradient (fine σ=0, coarse σ=1.5) combined with Laplacian
3. **Adaptive Threshold** — 68th-percentile threshold + binary closing (3 iterations, 8-connected structure)
4. **Cluster Extraction** — 8-connected flood-fill; clusters < 18 pixels are discarded
5. **Hyperbola Fitting** — `scipy.optimize.curve_fit()` with model `t = √(t₀² + ((x − x₀) / v)²)`, bounds enforced on velocity range 0.06–0.25 m/ns, R² filter ≥ 0.58
6. **Deduplication** — Detections within Δx < 0.35 m **and** Δdepth < 0.18 m of a higher-confidence result are removed

**Confidence score:** `0.6 × R² + 0.4 × min(1, cluster_size / 70)`

Without SciPy, a basic fallback scans amplitude variance across sliding windows and returns up to 35 candidates.

---

## Radargram Output

When matplotlib is available, each processed file generates a 4-panel 18×12 inch PNG at 400 dpi saved to `{output_folder}/Radargrams/`.

| Panel | Colormap | Description |
|---|---|---|
| Enhanced | Custom 9-stop GPR | Per-trace AGC (RMS), ±98th percentile contrast |
| Contrast | `seismic` | Polarity-sensitive, symmetric vmin/vmax |
| Edges | `hot` | Sobel + Laplace feature map normalised 0–1 |
| Detections | `gray` | Base image with confidence-coloured markers |

Detection markers: **Red** ≥ 0.80 · **Orange** 0.60–0.79 · **Yellow** < 0.60

---

## Velocity Reference

| Material | Velocity (m/ns) |
|---|---|
| Air | 0.30 |
| Dry sand | 0.15 |
| Concrete | 0.12 |
| Moist soil | 0.09–0.10 |
| Wet clay | 0.06–0.07 |
| Fresh water | 0.033 |

The tool reads dielectric permittivity (ε) from the DZT header and computes `v = 0.3 / √ε`. Override with the **Velocity** parameter if you know your site conditions.

---

## Known Limitations

- Trace cap of 3,000 per file (split long lines at the instrument level)
- Assumes GSSI `.DZT` format with 16-bit signed integer samples
- `.DZG` GPS parsing expects `$GSSIS` + `$GPGGA`/`$GNGGA` sentence pairs
- 3D (multi-channel) surveys are not currently supported

---

## Contributing

Issues, feature requests, and pull requests are welcome. Please open an Issue before starting large changes to confirm scope.

```bash
git clone https://github.com/micahwilliams/gpr-toolbox.git
cd gpr-toolbox
```

---

## License

MIT License — see [LICENSE](LICENSE) for details.

---

## Author

**Micah Williams**  
*Gpr-Toolbox*

---

*If this tool saves you time on a project, consider leaving a ⭐ on GitHub or a review on the Esri Marketplace.*
