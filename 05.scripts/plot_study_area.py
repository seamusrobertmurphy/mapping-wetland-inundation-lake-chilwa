#!/usr/bin/env python3
"""
Study-area and terrain figure for the Lake Chilwa manuscript.

Two panels from the self-derived terrain products (05.scripts/watershed-algorithms.qmd):
  (a) SRTM15+ elevation with hillshade, the basin and its sub-basins overlaid.
  (b) log10 D-infinity specific catchment area, the drainage structure that
      delineates the basin (Section 2.2.A), with the basin outline.

Local only, no Earth Engine: reads the DEM archives in 03.outputs/DEM and the
basin polygons in 03.outputs/SHP, matching plot_inundation_timeseries.py in style.

Output: 04.images/fig01_study_area.png / .pdf
Run: /opt/local/bin/python 05.scripts/plot_study_area.py
"""
import os
import tarfile
import numpy as np
import shapefile  # pyshp
import rasterio
from rasterio.windows import from_bounds
from pyproj import Geod
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource, LinearSegmentedColormap

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.normpath(os.path.join(HERE, "..", "03.outputs"))
IMG = os.path.normpath(os.path.join(HERE, "..", "04.images"))
DEM = os.path.join(OUT, "DEM")
EXT = os.path.join(DEM, "extracted")
SHP = os.path.join(OUT, "SHP")

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 10,
    "axes.titlesize": 10,
    "axes.titleweight": "bold",
    "axes.labelsize": 10,
    "axes.edgecolor": "#444444",
    "axes.linewidth": 0.8,
    "xtick.color": "#444444",
    "ytick.color": "#444444",
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})
BASIN_RED = "#b2182b"
SUBBASIN = "#33447066"

ARCHIVES = {"output_SRTM15Plus.tif": "rasters_SRTM15Plus",
            "Dinfareasca.tif": "Dinfarea"}


def _ensure_extracted():
    os.makedirs(EXT, exist_ok=True)
    for tif, arc in ARCHIVES.items():
        if os.path.exists(os.path.join(EXT, tif)):
            continue
        path = os.path.join(DEM, arc + ".tar.gz")
        if os.path.exists(path):
            with tarfile.open(path) as t:
                t.extractall(EXT)


def _read_window(tif, bounds):
    with rasterio.open(os.path.join(EXT, tif)) as s:
        win = from_bounds(*bounds, transform=s.transform)
        arr = s.read(1, window=win, masked=True).astype(float)
        w = s.window_transform(win)
        h, wd = arr.shape
        extent = [w.c, w.c + wd * w.a, w.f + h * w.e, w.f]  # l, r, b, t
    return np.ma.masked_invalid(arr), extent


def _rings(path):
    """Yield (xs, ys) for every ring of every shape in a shapefile."""
    r = shapefile.Reader(path)
    for sh in r.shapes():
        pts = sh.points
        parts = list(sh.parts) + [len(pts)]
        for i in range(len(parts) - 1):
            seg = pts[parts[i]:parts[i + 1]]
            if seg:
                xs, ys = zip(*seg)
                yield xs, ys


def _draw_outline(ax, path, **kw):
    for xs, ys in _rings(path):
        ax.plot(xs, ys, **kw)


def _scalebar(ax, lat, km=20):
    """A km scale bar sized for the local longitude degree length."""
    x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim()
    deg_km = 111.320 * np.cos(np.radians(lat))
    dx = km / deg_km
    sx = x0 + 0.06 * (x1 - x0); sy = y0 + 0.07 * (y1 - y0)
    ax.plot([sx, sx + dx], [sy, sy], color="#222222", lw=2.4,
            solid_capstyle="butt", zorder=6)
    ax.text(sx + dx / 2, sy + 0.012 * (y1 - y0), f"{km} km", ha="center",
            va="bottom", fontsize=7.5, color="#222222", zorder=6)


def build():
    _ensure_extracted()
    os.makedirs(IMG, exist_ok=True)
    basin = os.path.join(SHP, "chilwa_basin.shp")
    subs = os.path.join(SHP, "chilwa_subasins.shp")

    br = shapefile.Reader(basin).bbox  # xmin, ymin, xmax, ymax
    pad = 0.12
    bounds = (br[0] - pad, br[1] - pad, br[2] + pad, br[3] + pad)  # l,b,r,t
    lat_mid = (br[1] + br[3]) / 2

    elev, ext_e = _read_window("output_SRTM15Plus.tif", bounds)
    acc, ext_a = _read_window("Dinfareasca.tif", bounds)

    fig, (axa, axb) = plt.subplots(1, 2, figsize=(11, 5.4))

    # (a) elevation with hillshade. Truncate the terrain ramp to its land
    # portion so the low basin floor (~450 m) does not read as blue water.
    ls = LightSource(azdeg=315, altdeg=45)
    land = LinearSegmentedColormap.from_list(
        "land", plt.get_cmap("terrain")(np.linspace(0.25, 1.0, 256)))
    rgb = ls.shade(elev.filled(np.nan), cmap=land, blend_mode="soft",
                   vert_exag=0.0009, dx=1, dy=1,
                   vmin=float(elev.min()), vmax=float(elev.max()))
    axa.imshow(rgb, extent=ext_e, origin="upper")
    im_e = axa.imshow(elev, extent=ext_e, cmap=land, alpha=0,
                      vmin=float(elev.min()), vmax=float(elev.max()))
    _draw_outline(axa, subs, color=SUBBASIN, lw=0.4, zorder=3)
    _draw_outline(axa, basin, color=BASIN_RED, lw=1.8, zorder=4)
    cb = fig.colorbar(im_e, ax=axa, fraction=0.046, pad=0.02)
    cb.set_label("Elevation (m)", fontsize=8.5); cb.ax.tick_params(labelsize=7.5)
    axa.set_title("Lake Chilwa basin and terrain", loc="left")
    axa.set_xlabel("Longitude (°E)"); axa.set_ylabel("Latitude (°)")
    _scalebar(axa, lat_mid)
    axa.text(-0.02, 1.03, "a", transform=axa.transAxes, fontsize=11,
             fontweight="bold", va="bottom", ha="right")

    # (b) log10 D-infinity contributing area
    logacc = np.ma.masked_less_equal(acc, 0)
    logacc = np.ma.log10(logacc)
    blues = LinearSegmentedColormap.from_list(
        "acc", ["#f7fbff", "#9ecae1", "#4292c6", "#08306b"])
    im_a = axb.imshow(logacc, extent=ext_a, cmap=blues,
                      vmin=float(logacc.min()), vmax=float(logacc.max()))
    _draw_outline(axb, basin, color=BASIN_RED, lw=1.8, zorder=4)
    cb2 = fig.colorbar(im_a, ax=axb, fraction=0.046, pad=0.02)
    cb2.set_label("log₁₀ contributing area", fontsize=8.5)
    cb2.ax.tick_params(labelsize=7.5)
    axb.set_title("D-infinity drainage structure", loc="left")
    axb.set_xlabel("Longitude (°E)"); axb.set_ylabel("Latitude (°)")
    _scalebar(axb, lat_mid)
    axb.text(-0.02, 1.03, "b", transform=axb.transAxes, fontsize=11,
             fontweight="bold", va="bottom", ha="right")

    for ax in (axa, axb):
        ax.set_xlim(bounds[0], bounds[2]); ax.set_ylim(bounds[1], bounds[3])
        ax.tick_params(length=3)

    for e in ("png", "pdf"):
        fig.savefig(os.path.join(IMG, f"fig01_study_area.{e}"))
    plt.close(fig)

    geod = Geod(ellps="WGS84")
    area_km2 = 0.0
    for xs, ys in _rings(basin):
        a, _ = geod.polygon_area_perimeter(xs, ys)
        area_km2 += abs(a) / 1e6
    print("  wrote fig01_study_area.png / .pdf ->", IMG)
    print(f"  basin geodesic area: {area_km2:,.0f} km2  "
          f"(elev {elev.min():.0f}-{elev.max():.0f} m; 71 sub-basins)")


if __name__ == "__main__":
    build()
