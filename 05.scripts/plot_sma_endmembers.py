#!/usr/bin/env python3
"""
Spectral-mixture endmember figure for the Lake Chilwa manuscript.

Draws the four image-derived endmember reflectance spectra beside example
wet-year and dry-year open-water fraction maps. It is a post-processing step
over summaries exported by the `sma-endmembers`, `sma-timeseries`, and
`sma-export` chunks, so it runs locally without Earth Engine, matching
plot_inundation_timeseries.py.

Inputs
  03.outputs/sma_endmember_spectra.csv
      one row per class, columns: class, blue, green, red, nir, swir1, swir2
      surface reflectance (0-1), the per-class median spectrum.
  03.outputs/sma_fractions_wet.tif
  03.outputs/sma_fractions_dry.tif
      open-water fraction rasters (0-1) for a wet and a dry reference year,
      clipped to the basin, in geographic coordinates.
  03.outputs/chilwa_basin.geojson
      basin outline drawn over both maps (optional).

Output
  04.images/fig07_sma_endmembers.png / .pdf

Needs matplotlib and rasterio; on this machine that is /opt/local/bin/python.
Run: /opt/local/bin/python 05.scripts/plot_sma_endmembers.py
"""
import json
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.normpath(os.path.join(HERE, "..", "03.outputs"))
IMG = os.path.normpath(os.path.join(HERE, "..", "04.images"))

# ----------------------------------------------------------------------------
# Shared publication theme (identical to plot_inundation_timeseries.py)
# ----------------------------------------------------------------------------
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 10,
    "axes.titlesize": 10,
    "axes.titleweight": "bold",
    "axes.labelsize": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.edgecolor": "#444444",
    "axes.linewidth": 0.8,
    "xtick.color": "#444444",
    "ytick.color": "#444444",
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 8.5,
    "legend.frameon": False,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})

# Water blue and emergent green as in the sibling figures; dry vegetation and
# bare soil in the earth tones already carried by the L-band series.
CLASS_COLOUR = {
    "water": "#2166ac", "flooded_veg": "#5aae61",
    "dry_veg": "#b8860b", "bare_soil": "#8c510a",
}
CLASS_LABEL = {
    "water": "Open water",
    "flooded_veg": "Emergent / flooded vegetation",
    "dry_veg": "Dry vegetation",
    "bare_soil": "Bare soil",
}
CLASS_ORDER = ["water", "flooded_veg", "dry_veg", "bare_soil"]

# Landsat OLI band centres (micrometres) for the spectral-response axis.
BANDS = ["blue", "green", "red", "nir", "swir1", "swir2"]
WAVELENGTH = {"blue": 0.48, "green": 0.56, "red": 0.655,
              "nir": 0.865, "swir1": 1.61, "swir2": 2.20}

# Water-fraction ramp, matching the manuscript's own SMA visualisation.
FRAC_CMAP = LinearSegmentedColormap.from_list(
    "sma_water", ["white", "#00ffff", "#2166ac", "#08306b"])
GRID = "#e6e6e6"
OUTLINE = "#b2182b"


def _grid(ax, axis="y"):
    ax.grid(axis=axis, color=GRID, linewidth=0.7, zorder=0)
    ax.set_axisbelow(True)


def _panel_tag(ax, tag):
    ax.text(-0.02, 1.03, tag, transform=ax.transAxes, fontsize=11,
            fontweight="bold", va="bottom", ha="right")


def _geoms(path):
    """Polygon geometries of a GeoJSON layer."""
    if not os.path.exists(path):
        return []
    with open(path) as fh:
        gj = json.load(fh)
    return [f.get("geometry", f) for f in gj.get("features", [gj])]


def _rings(geoms):
    """Exterior rings of the geometries, as arrays of lon/lat."""
    out = []
    for geom in geoms:
        polys = ([geom["coordinates"]] if geom["type"] == "Polygon"
                 else geom["coordinates"])
        for poly in polys:
            out.append(np.asarray(poly[0], dtype=float))
    return out


def _read_fraction(path, geoms):
    """Open-water fraction of a raster, masked to the basin, with its extent.

    The Earth Engine download covers the bounding box rather than the clip
    geometry, so neighbouring water bodies outside the divide, Lake Chiuta to
    the north among them, are removed here.
    """
    import rasterio
    from rasterio.features import geometry_mask
    with rasterio.open(path) as src:
        names = list(src.descriptions or ())
        band = names.index("water_frac") + 1 if "water_frac" in names else 1
        arr = src.read(band, masked=True).astype(float)
        b = src.bounds
        outside = (geometry_mask(geoms, arr.shape, src.transform)
                   if geoms else np.zeros(arr.shape, bool))
    arr = np.ma.masked_invalid(arr)
    return np.ma.masked_where(outside, arr), (b.left, b.right, b.bottom, b.top)


def _map_panel(ax, arr, extent, rings, title, tag):
    im = ax.imshow(arr, cmap=FRAC_CMAP, vmin=0, vmax=1, extent=extent,
                   origin="upper", interpolation="nearest")
    for ring in rings:
        ax.plot(ring[:, 0], ring[:, 1], color=OUTLINE, lw=0.9, zorder=3)
    lat = 0.5 * (extent[2] + extent[3])
    ax.set_aspect(1.0 / np.cos(np.deg2rad(lat)))
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(True)
        s.set_edgecolor("#444444")
        s.set_linewidth(0.8)
    ax.set_title(title, loc="left", fontsize=9.5)
    _panel_tag(ax, tag)
    return im


def build(out_dir=OUT, img_dir=IMG,
          wet_label="Wet year, 2023", dry_label="Dry year, 2018"):
    os.makedirs(img_dir, exist_ok=True)
    spectra_csv = os.path.join(out_dir, "sma_endmember_spectra.csv")
    wet_tif = os.path.join(out_dir, "sma_fractions_wet.tif")
    dry_tif = os.path.join(out_dir, "sma_fractions_dry.tif")

    missing = [p for p in (spectra_csv, wet_tif, dry_tif) if not os.path.exists(p)]
    if missing:
        rel = [os.path.relpath(p, os.path.dirname(HERE)) for p in missing]
        raise SystemExit(
            "Missing SMA outputs, so the endmember figure cannot be built yet:\n  "
            + "\n  ".join(rel)
            + "\nRun the sma-export chunk under a live Earth Engine session, "
              "then rerun this script.")

    spec = pd.read_csv(spectra_csv)
    # Endmembers are derived per sensor era; the figure shows the Landsat 8 set
    # from the 2020 reference composite that the classification also uses.
    if "era" in spec.columns:
        spec = spec[spec["era"] == "L8"]
    spec = spec.set_index("class")
    geoms = _geoms(os.path.join(out_dir, "chilwa_basin.geojson"))
    wet, wet_ext = _read_fraction(wet_tif, geoms)
    dry, dry_ext = _read_fraction(dry_tif, geoms)
    rings = _rings(geoms)

    fig = plt.figure(figsize=(11.2, 4.1))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.55, 1, 1],
                          height_ratios=[1, 0.05], wspace=0.2, hspace=0.1)

    # (a) endmember reflectance spectra. Bands sit at even spacing rather than
    # true wavelength so the tightly grouped visible bands do not collide; the
    # centre wavelength is carried in each tick label.
    axs = fig.add_subplot(gs[0, 0])
    x = list(range(len(BANDS)))
    for cls in CLASS_ORDER:
        if cls not in spec.index:
            continue
        y = spec.loc[cls, BANDS].to_numpy(dtype=float)
        axs.plot(x, y, "-o", color=CLASS_COLOUR[cls], ms=4, lw=1.8,
                 label=CLASS_LABEL[cls], zorder=3)
    _grid(axs)
    for xv in x:
        axs.axvline(xv, color=GRID, lw=0.6, zorder=0)
    axs.set_xlabel("Band (centre wavelength, µm)")
    axs.set_ylabel("Surface reflectance")
    axs.set_ylim(bottom=0)
    axs.set_xlim(-0.3, len(BANDS) - 0.7)
    axs.set_xticks(x)
    axs.set_xticklabels([f"{b.upper()}\n{WAVELENGTH[b]:.2f}" for b in BANDS],
                        fontsize=8)
    axs.legend(loc="upper left")
    axs.set_title("Image-derived endmember spectra", loc="left")
    _panel_tag(axs, "a")

    # (b, c) wet-year and dry-year open-water fraction, on a shared scale
    axw = fig.add_subplot(gs[0, 1])
    axd = fig.add_subplot(gs[0, 2])
    _map_panel(axw, wet, wet_ext, rings, wet_label, "b")
    im = _map_panel(axd, dry, dry_ext, rings, dry_label, "c")

    cax = fig.add_subplot(gs[1, 1:])
    cb = fig.colorbar(im, cax=cax, orientation="horizontal")
    cb.set_label("Open-water fraction", fontsize=8.5)
    cb.outline.set_linewidth(0.6)
    cb.ax.tick_params(labelsize=8)

    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(img_dir, f"fig07_sma_endmembers.{ext}"))
    plt.close(fig)
    print("  wrote fig07_sma_endmembers.png / .pdf ->", img_dir)

    print("\n================ ENDMEMBER SPECTRA ================")
    print("  " + f"{'class':<30}" + "  ".join(f"{b:>6}" for b in BANDS))
    for cls in CLASS_ORDER:
        if cls in spec.index:
            vals = "  ".join(f"{spec.loc[cls, b]:6.3f}" for b in BANDS)
            print(f"  {CLASS_LABEL[cls]:<30}{vals}")
    for arr, lab in ((wet, wet_label), (dry, dry_label)):
        print(f"  {lab}: mean water fraction {arr.mean():.3f}, "
              f"fraction of basin above 0.5 = {(arr > 0.5).sum() / arr.count():.3f}")
    print("==================================================\n")


if __name__ == "__main__":
    build()
