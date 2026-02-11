# -*- coding: utf-8 -*-
"""
GPR Utility Detection  - Output SR + No Dummy Geometry + SR Fix
- Output Coordinate System parameter (GPCoordinateSystem)
- Removes dummy path/point geometries when GPS is absent (NULL geometry instead)
- FIX: Coerce Output Coordinate System to arcpy.SpatialReference even if WKT/.prj string is provided

Notes:
- PNG radargram export requires 'matplotlib' installed in your ArcGIS Pro Python env.
- If PNGs are not created, ensure you cloned Pro’s environment and added matplotlib.
"""

import arcpy
from pathlib import Path
import numpy as np
import struct

# Optional but highly recommended
try:
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    HAS_MPL = True
except:
    HAS_MPL = False

try:
    from scipy.optimize import curve_fit
    from scipy.ndimage import gaussian_filter1d, sobel, laplace, generate_binary_structure, binary_closing
    from scipy.interpolate import interp1d
    from scipy.signal import medfilt
    HAS_SCIPY = True
except:
    HAS_SCIPY = False
    arcpy.AddWarning("SciPy missing – basic mode only")

# --- helper to coerce parameter to SpatialReference ---
def _to_spatial_reference(sr_value):
    """Return arcpy.SpatialReference from a parameter value that may be SR, WKT, or .prj path."""
    if isinstance(sr_value, arcpy.SpatialReference):
        return sr_value
    if sr_value is None:
        return arcpy.SpatialReference(4326)
    sr = arcpy.SpatialReference()
    # Try interpreting the value as WKT
    try:
        sr.loadFromString(str(sr_value))
        if sr.factoryCode != 0 or sr.name:
            return sr
    except Exception:
        pass
    # Try if it's a path to a .prj file
    try:
        prj_path = Path(str(sr_value))
        if prj_path.exists() and prj_path.suffix.lower() == '.prj':
            sr_text = prj_path.read_text(encoding='utf-8', errors='ignore')
            sr.loadFromString(sr_text)
            return sr
    except Exception:
        pass
    # Fallback to WGS 84
    return arcpy.SpatialReference(4326)


class HyperbolaDetector:
    def __init__(self):
        self.velocity_range = (0.06, 0.25)

    def detect_hyperbolas(self, data, metadata):
        if not HAS_SCIPY:
            return self._basic_detection(data, metadata)

        enhanced = self._tuned_preprocess(data)
        edges = self._enhance_features(enhanced)
        mask = self._adaptive_threshold(edges)
        clusters = self._cluster_columns(mask)

        hyps = []
        for cluster in clusters:
            hyp = self._fit_hyperbola(cluster, enhanced, metadata)
            if hyp and hyp['r_squared'] > 0.58:
                hyps.append(hyp)

        return self._remove_duplicates(sorted(hyps, key=lambda x: -x['confidence']))

    def _tuned_preprocess(self, data):
        d = medfilt(data.astype(np.float64), kernel_size=3)
        d -= np.mean(d, axis=0)

        # Adaptive gain – safe for any number of samples
        nsamples = d.shape[0]
        gain = np.ones(nsamples)
        if nsamples > 100:
            deep_gain = np.linspace(1.0, 5.0, nsamples - 100)
            gain[100:] *= deep_gain
        d *= gain[:, np.newaxis]

        p2, p98 = np.percentile(d, [1, 99])
        d = np.clip(d, p2, p98)
        return d

    def _enhance_features(self, data):
        fine = np.sqrt(sobel(data, axis=0)**2 + sobel(data, axis=1)**2)
        smooth = gaussian_filter1d(data, sigma=1.5, axis=0)
        coarse = np.sqrt(sobel(smooth, axis=0)**2 + sobel(smooth, axis=1)**2)
        edges = 0.7*fine + 0.3*coarse
        lap = np.abs(laplace(data))
        return edges + 0.4*lap

    def _adaptive_threshold(self, edges):
        thresh = np.percentile(edges, 68)
        mask = edges > thresh
        return binary_closing(mask, generate_binary_structure(2, 2), iterations=3)

    def _cluster_columns(self, mask):
        h, w = mask.shape
        visited = np.zeros_like(mask, bool)
        clusters = []
        for col in range(w):
            for row in np.where(mask[:, col] & ~visited[:, col])[0]:
                if visited[row, col]:
                    continue
                cluster = []
                stack = [(row, col)]
                while stack:
                    y, x = stack.pop()
                    if not (0 <= y < h and 0 <= x < w) or visited[y, x] or not mask[y, x]:
                        continue
                    visited[y, x] = True
                    cluster.append((y, x))
                    for dy in [-1, 0, 1]:
                        for dx in [-1, 0, 1]:
                            if dy == dx == 0:
                                continue
                            stack.append((y + dy, x + dx))
                if len(cluster) >= 18:
                    clusters.append(cluster)
        return clusters

    def _fit_hyperbola(self, cluster, data, md):
        rows, cols = np.array(cluster).T
        x = cols * md['units_per_scan']
        t = rows * md['sample_rate_ns']

        def model(x, x0, t0, v):
            return np.sqrt(t0**2 + ((x - x0) / v)**2)

        try:
            p0 = [np.mean(x), np.min(t) + 0.5, md['velocity']]
            bounds = ([min(x) - 1, 0.1, self.velocity_range[0]],
                      [max(x) + 1, max(t), self.velocity_range[1]])
            popt, _ = curve_fit(model, x, t, p0=p0, bounds=bounds, maxfev=4000)
            x0, t0, v = popt
            pred = model(x, *popt)
            ss_res = np.sum((t - pred)**2)
            ss_tot = np.sum((t - np.mean(t))**2)
            r2 = 1 - ss_res/ss_tot if ss_tot > 0 else 0
            if r2 > 0.58:
                depth = v * t0 / 2.0
                conf = 0.6*r2 + 0.4*min(1, len(cluster)/70)

                # Store the trace index for accurate plotting
                trace_idx = int(round(x0 / md['units_per_scan']))

                return {'x0': float(x0), 'depth': float(depth), 'velocity': float(v),
                        'r_squared': float(r2), 'confidence': float(conf),
                        'diameter_estimate': max(0.03, min(0.8, depth*0.45)),
                        'trace_idx': trace_idx}
        except:
            pass
        return None

    def _remove_duplicates(self, hyps):
        good = []
        for h in hyps:
            if not any(abs(h['x0'] - e['x0']) < 0.35 and abs(h['depth'] - e['depth']) < 0.18 for e in good):
                good.append(h)
        return good

    def _basic_detection(self, data, md):
        hyps = []
        for col in range(10, data.shape[1] - 10, 6):
            win = data[:, col - 5:col + 6]
            if np.std(win) > np.mean(np.abs(win)) * 3:
                row = np.argmax(np.abs(win[:, 5])) + 15
                depth = row * md['sample_rate_ns'] * md['velocity'] / 2
                hyps.append({'x0': col * md['units_per_scan'], 'depth': depth,
                             'confidence': 0.68, 'r_squared': 0.62, 'velocity': md['velocity'],
                             'diameter_estimate': 0.18, 'trace_idx': col})
        return hyps[:35]


class DZTParser:
    def __init__(self):
        self.raw = self.processed = None
        self.metadata = {}
        self.gps_lat = self.gps_lon = []

    def read_dzt(self, path):
        with open(path, 'rb') as f:
            h = f.read(1024)
            extended = struct.unpack('<h', h[112:114])[0]
            hdr_size = 131072 if extended else 1024
            samples = struct.unpack('<h', h[4:6])[0]
            dielectric = struct.unpack('<f', h[54:58])[0]
            vel = 0.09 if dielectric <= 0 else 0.3/np.sqrt(max(dielectric, 1))
            f.seek(hdr_size)
            raw = np.frombuffer(f.read(), np.int16)[::2]
            traces = min(3000, len(raw)//samples)
            if len(raw) < traces * samples:
                return False
            self.raw = raw[:traces * samples].reshape((traces, samples)).T
            self.metadata = {'samples_per_trace': samples, 'traces': traces,
                             'sample_rate_ns': 0.111828, 'units_per_scan': 0.02,
                             'velocity': vel, 'has_gps': False}
            return True

    def read_dzg(self, path):
        if not Path(path).exists():
            return False
        try:
            text = Path(path).read_text(errors='ignore').splitlines()
            blocks = []
            i = 0
            while i < len(text) - 1:
                if text[i].startswith('$GSSIS'):
                    parts = text[i].split(',')
                    if len(parts) > 1:
                        try:
                            trace = int(parts[1])
                            gga = text[i + 1]
                            if gga.startswith(('$GPGGA', '$GNGGA')):
                                p = gga.split(',')
                                if len(p) >= 10:
                                    lat = float(p[2][:2]) + float(p[2][2:]) / 60
                                    lon = float(p[4][:3]) + float(p[4][3:]) / 60
                                    if p[3] in 'Ss': lat = -lat
                                    if p[5] in 'Ww': lon = -lon
                                    blocks.append((trace, lat, lon))
                            i += 2
                            continue
                        except:
                            pass
                i += 1
            if blocks:
                blocks.sort()
                tr, la, lo = zip(*blocks)
                interp_lat = interp1d(tr, la, kind='linear', fill_value='extrapolate')
                interp_lon = interp1d(tr, lo, kind='linear', fill_value='extrapolate')
                idx = np.arange(self.metadata['traces'])
                self.gps_lat = interp_lat(idx)
                self.gps_lon = interp_lon(idx)
                diffs = np.sqrt(np.diff(self.gps_lon)**2 + np.diff(self.gps_lat)**2) * 111194
                self.metadata['trace_distances'] = np.insert(np.cumsum(diffs), 0, 0)
                self.metadata['has_gps'] = True
                arcpy.AddMessage(f"DZG parsed → {len(blocks)} points → {self.metadata['traces']} traces")
                return True
        except Exception as e:
            arcpy.AddWarning(f"DZG error: {e}")
        return False

    def preprocess(self):
        if self.raw is None:
            return
        d = self.raw.astype(float)
        d -= np.mean(d, axis=0)
        gain = np.linspace(0.8, 4.0, d.shape[0])
        d *= gain[:, None]
        p2, p98 = np.percentile(d, [2, 98])
        d = np.clip(d, p2, p98)
        self.processed = d

    def depth_axis(self):
        return np.arange(self.metadata['samples_per_trace']) * self.metadata['sample_rate_ns'] * self.metadata['velocity'] / 2


def generate_radargram(parser, hyps, out_path, title):
    """
    Radargram export:
    - Crisp pixels (no interpolation/resampling)
    - Per-trace display AGC (RMS normalization)
    - Symmetric contrast (±98th percentile)
    - Edges panel if SciPy present
    - 400 dpi export
    """
    if not HAS_MPL:
        return
    try:
        cmap = LinearSegmentedColormap.from_list(
            'gpr',
            ['#000033', '#0000BB', '#0E4C92', '#2E8B57', '#90EE90',
             '#FFFF00', '#FFA500', '#FF4500', '#8B0000'],
            N=256
        )

        data = parser.processed
        depth = parser.depth_axis()

        # Distance axis
        if 'trace_distances' in parser.metadata:
            dist_array = parser.metadata['trace_distances']
            max_dist = float(dist_array[-1]) if len(dist_array) else 0.0
        else:
            max_dist = parser.metadata['traces'] * parser.metadata['units_per_scan']
            dist_array = np.linspace(0, max_dist, parser.metadata['traces'])

        # Display-only AGC
        display = data.copy()
        rms = np.sqrt(np.mean(display**2, axis=0)) + 1e-9
        display = display / rms

        vmax = np.percentile(display, 98)
        if vmax <= 0:
            vmax = np.max(np.abs(display)) + 1e-6
        vmin = -vmax

        # Edges panel
        edges_img = None
        if HAS_SCIPY:
            try:
                g0 = sobel(display, axis=0)
                g1 = sobel(display, axis=1)
                edges_img = np.sqrt(g0**2 + g1**2)
                if np.max(edges_img) > 0:
                    edges_img = edges_img / np.max(edges_img)
            except:
                edges_img = None

        fig, axs = plt.subplots(2, 2, figsize=(18, 12))
        panels = [
            ("Enhanced", cmap, display, dict(vmin=vmin, vmax=vmax)),
            ("Contrast", 'seismic', display, dict(vmin=vmin, vmax=vmax)),
            ("Edges", 'hot', edges_img if edges_img is not None else display, {}),
            ("Detections", 'gray', display, {})
        ]

        for i, (name, cm, img, extra_kwargs) in enumerate(panels):
            ax = axs[i // 2, i % 2]
            im = ax.imshow(
                img, cmap=cm, aspect='auto',
                extent=[0, max_dist, depth[-1], depth[0]],
                interpolation='none', resample=False, **extra_kwargs
            )
            ax.set_title(f"{name} – {len(hyps)} utilities")
            ax.set_xlabel("Distance (m)")
            ax.set_ylabel("Depth (m)")

            # Plot detections
            for h in hyps:
                c = 'red' if h['confidence'] > 0.8 else ('orange' if h['confidence'] > 0.6 else 'yellow')
                if 'trace_idx' in h and 0 <= h['trace_idx'] < len(dist_array):
                    plot_x = float(dist_array[h['trace_idx']])
                else:
                    plot_x = h['x0']
                ax.plot(plot_x, h['depth'], 'o', c=c, ms=12, mec='white', mew=2)
                ax.axvline(x=plot_x, color=c, alpha=0.4, linewidth=1.2, linestyle='--')
                ax.axhline(y=h['depth'], color=c, alpha=0.4, linewidth=1.2, linestyle='--')

            plt.colorbar(im, ax=ax, fraction=0.046)

        try:
            from matplotlib.lines import Line2D
            handles = [
                Line2D([0], [0], marker='o', color='w', label='≥ 0.80',
                       markerfacecolor='red', markeredgecolor='white', markersize=10),
                Line2D([0], [0], marker='o', color='w', label='0.60–0.79',
                       markerfacecolor='orange', markeredgecolor='white', markersize=10),
                Line2D([0], [0], marker='o', color='w', label='< 0.60',
                       markerfacecolor='yellow', markeredgecolor='white', markersize=10),
            ]
            axs[0, 1].legend(handles=handles, title='Confidence', loc='upper right')
        except:
            pass

        plt.suptitle(f"{title} - V9.1C (Crisp Plotting)", fontsize=16)
        plt.tight_layout()
        plt.savefig(out_path, dpi=400, bbox_inches='tight')
        plt.close()

    except Exception as e:
        arcpy.AddWarning(f"Radargram error: {e}")


class Toolbox(object):
    def __init__(self):
        self.label = "GPR Utility Detection PRO v9.1B (SR + NoDummyXY + SRFix)"
        self.tools = [GPRTool]


class GPRTool(object):
    def __init__(self):
        self.label = "GPR Utility Detection Tool"
        self.description = "Automated hyperbola detection • Output CRS selectable • No dummy XY when GPS is absent"

    def getParameterInfo(self):
        p0 = arcpy.Parameter(displayName="Input Folder (.DZT files)", name="in_folder",
                             datatype="DEFolder", parameterType="Required", direction="Input")
        p1 = arcpy.Parameter(displayName="Output Folder", name="out_folder",
                             datatype="DEFolder", parameterType="Required", direction="Output")
        p2 = arcpy.Parameter(displayName="Survey Name", name="survey_name",
                             datatype="GPString", parameterType="Optional", direction="Input"); p2.value = "GPR_Survey"
        p3 = arcpy.Parameter(displayName="Velocity (m/ns) – blank = auto (0.09 urban default)", name="velocity",
                             datatype="GPDouble", parameterType="Optional", direction="Input")
        p4 = arcpy.Parameter(displayName="Min Confidence (try 0.55 for very weak pipes)", name="min_conf",
                             datatype="GPDouble", parameterType="Optional", direction="Input"); p4.value = 0.58
        p5 = arcpy.Parameter(displayName="Auto-add to map", name="auto_add",
                             datatype="GPBoolean", parameterType="Optional", direction="Input"); p5.value = True
        p6 = arcpy.Parameter(displayName="Output Coordinate System", name="out_sr",
                             datatype="GPCoordinateSystem", parameterType="Optional", direction="Input")

        # Default to the active map's SR (fallback WGS 84)
        try:
            aprx = arcpy.mp.ArcGISProject("CURRENT"); m = aprx.activeMap or aprx.listMaps()[0]
            if m and m.spatialReference and m.spatialReference.factoryCode != 0:
                p6.value = m.spatialReference
            else:
                p6.value = arcpy.SpatialReference(4326)
        except:
            p6.value = arcpy.SpatialReference(4326)

        return [p0, p1, p2, p3, p4, p5, p6]

    def isLicensed(self): 
        return True

    def execute(self, parameters, messages):
        in_folder = Path(parameters[0].valueAsText)
        out_folder = Path(parameters[1].valueAsText)
        survey = parameters[2].valueAsText or "GPR_Survey"
        vel = parameters[3].value
        minconf = parameters[4].value or 0.58
        autoadd = parameters[5].value
        out_sr = parameters[6].value

        sr = _to_spatial_reference(out_sr)
        arcpy.AddMessage(f"Output SR: {sr.name} (WKID: {getattr(sr, 'factoryCode', 0)})")

        out_folder.mkdir(parents=True, exist_ok=True)
        (out_folder / "Radargrams").mkdir(exist_ok=True)

        utils_name = f"{survey}_Utilities"
        path_name  = f"{survey}_Path"

        utils_shp = out_folder / f"{utils_name}.shp"
        path_shp  = out_folder / f"{path_name}.shp"

        for s in (utils_shp, path_shp):
            if arcpy.Exists(str(s)):
                arcpy.Delete_management(str(s))

        # Create outputs in chosen SR
        arcpy.CreateFeatureclass_management(str(out_folder), utils_name, "POINT", spatial_reference=sr)
        for f in ["UtilID", "Depth_m", "Diam_m", "Confidence", "R_Squared", "Velocity", "SourceFile", "GPS_Avail", "XY_Valid"]:
            arcpy.AddField_management(str(utils_shp), f,
                                      "DOUBLE" if f not in ["UtilID", "SourceFile", "GPS_Avail", "XY_Valid"] else "TEXT",
                                      field_length=50)

        arcpy.CreateFeatureclass_management(str(out_folder), path_name, "POLYLINE", spatial_reference=sr)
        arcpy.AddField_management(str(path_shp), "FileName", "TEXT", field_length=100)
        arcpy.AddField_management(str(path_shp), "GPS_Avail", "TEXT", field_length=5)

        det = HyperbolaDetector()
        total = 0

        with arcpy.da.InsertCursor(str(utils_shp),
                                   ["SHAPE@XY", "UtilID", "Depth_m", "Diam_m", "Confidence", "R_Squared",
                                    "Velocity", "SourceFile", "GPS_Avail", "XY_Valid"]) as ucur, \
             arcpy.da.InsertCursor(str(path_shp),
                                   ["SHAPE@", "FileName", "GPS_Avail"]) as pcur:

            dzt_files = sorted(set(in_folder.glob("*.[Dd][Zz][Tt]")))
            if not dzt_files:
                messages.addMessage(f"ERROR: No .DZT files found in {in_folder}")
                return

            messages.addMessage(f"Found {len(dzt_files)} DZT files to process\n")

            for dzt in dzt_files:
                messages.addMessage(f"Processing {dzt.name}")
                p = DZTParser()
                if not p.read_dzt(str(dzt)):
                    continue
                p.read_dzg(dzt.with_suffix('.DZG') if dzt.with_suffix('.DZG').exists() else dzt.with_suffix('.dzg'))
                if vel:
                    p.metadata['velocity'] = vel
                p.preprocess()

                if 'trace_distances' in p.metadata:
                    diffs = np.diff(p.metadata['trace_distances'])
                    if len(diffs) > 0:
                        p.metadata['units_per_scan'] = float(np.mean(diffs))

                hyps = [h for h in det.detect_hyperbolas(p.processed, p.metadata) if h['confidence'] >= minconf]
                messages.addMessage(f"  → {len(hyps)} utilities found")

                # Radargram export
                if HAS_MPL:
                    png_dir = out_folder / "Radargrams"
                    generate_radargram(p, hyps, png_dir / f"{dzt.stem}_V9_1B.png", dzt.stem)

                has_gps = bool(p.metadata.get('has_gps', False) and len(getattr(p, 'gps_lat', [])) == p.metadata['traces'])

                # PATH: write geometry when GPS is available; else NULL geometry row
                if has_gps:
                    arr = arcpy.Array([arcpy.Point(lo, la) for la, lo in zip(p.gps_lat, p.gps_lon)])
                    pcur.insertRow([arcpy.Polyline(arr, sr), dzt.name, "Yes"])
                else:
                    pcur.insertRow([None, dzt.name, "No"])

                # POINTS: write XY when GPS is available; else NULL geometry with attributes
                for i, h in enumerate(hyps):
                    if 'trace_idx' in h:
                        idx = min(h['trace_idx'], p.metadata['traces']-1)
                    else:
                        idx = min(int(round(h['x0']/p.metadata['units_per_scan'])), p.metadata['traces']-1)

                    if has_gps:
                        lat = float(p.gps_lat[idx]); lon = float(p.gps_lon[idx])
                        ucur.insertRow([(lon, lat), f"{dzt.stem}_U{i+1}", h['depth'], h['diameter_estimate'],
                                        h['confidence'], h['r_squared'], h['velocity'], dzt.name, "Yes", "Yes"])
                        total += 1
                    else:
                        ucur.insertRow([None, f"{dzt.stem}_U{i+1}", h['depth'], h['diameter_estimate'],
                                        h['confidence'], h['r_squared'], h['velocity'], dzt.name, "No", "No"])

        if autoadd:
            try:
                aprx = arcpy.mp.ArcGISProject("CURRENT")
                m = aprx.activeMap or aprx.listMaps()[0]
                m.addDataFromPath(str(utils_shp))
                m.addDataFromPath(str(path_shp))
            except:
                pass

        messages.addMessage(f"DONE! {total} utilities mapped with valid XY; others stored as NULL geometry (no GPS)")
