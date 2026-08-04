#!/usr/bin/env python3
"""High-resolution RGB for locating the field photographs.

Two jobs need different imagery and conflating them wastes effort.

  Locating a photograph. Needs the finest resolution available; the date
  barely matters, because the things a photograph can be matched against
  (village layout, river mouths, roads, the hill skyline, Nchisi Island)
  persist across decades. Sentinel-2 at 10 m is 50 times finer than the
  best 2012 optical and the basin has over a thousand near-cloudless scenes.

  Validating a 2012 classification. Needs contemporaneous imagery, where
  500 m MODIS is the free ceiling (see export_rgb_photo_dates.py).

This script serves the first job, and adds the nearest 30 m Landsat either
side of the fieldwork, which is the finest contemporaneous-ish optical there is.

Exports are tiled because a 10 m mosaic over the lake exceeds Earth Engine's
50 MB single-request ceiling, then mosaicked locally.

Writes 03.outputs/TIF/RGB_highres/
"""

import io
import json
import os
import zipfile

import ee
import requests

ee.Initialize(project="murphys-deforisk")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "03.outputs", "TIF", "RGB_highres")
os.makedirs(OUT, exist_ok=True)

gj = json.load(open(os.path.join(ROOT, "03.outputs", "JSON", "chilwa_basin.geojson")))
geom = gj["features"][0]["geometry"] if gj.get("type") == "FeatureCollection" else gj
basin = ee.Geometry(geom)

# The photographs were all taken on or near the water, so the lake and its
# margin is the placement area. The full basin at 10 m is needlessly large:
# a 45 km buffer runs to 6,257 km2 and over 2 GB of tiles once the download
# grid is counted, against 2,434 km2 here, which still clears the reed belt
# by a wide margin on every shore.
lake_area = basin.centroid(1).buffer(28000).intersection(basin, 1)


def tiled_download(image, bands, scale, region, name, n=4):
    """Download in an n x n grid and report each tile. Returns paths written."""
    b = region.bounds().getInfo()["coordinates"][0]
    xs = [p[0] for p in b]
    ys = [p[1] for p in b]
    x0, x1, y0, y1 = min(xs), max(xs), min(ys), max(ys)
    paths = []
    for i in range(n):
        for j in range(n):
            gx0 = x0 + (x1 - x0) * i / n
            gx1 = x0 + (x1 - x0) * (i + 1) / n
            gy0 = y0 + (y1 - y0) * j / n
            gy1 = y0 + (y1 - y0) * (j + 1) / n
            reg = ee.Geometry.Rectangle([gx0, gy0, gx1, gy1])
            p = os.path.join(OUT, f"{name}_t{i}{j}.tif")
            if os.path.exists(p) and os.path.getsize(p) > 2000:
                paths.append(p)
                continue
            try:
                url = image.select(bands).clip(reg).getDownloadURL(
                    {"scale": scale, "region": reg, "format": "GEO_TIFF", "crs": "EPSG:4326"}
                )
                r = requests.get(url, timeout=900)
                if r.status_code == 200 and len(r.content) > 2000:
                    with open(p, "wb") as fh:
                        fh.write(r.content)
                    paths.append(p)
                    print(f"    {os.path.basename(p)}: {len(r.content)/1e6:.1f} MB", flush=True)
                else:
                    print(f"    {os.path.basename(p)}: HTTP {r.status_code} {r.content[:80]}", flush=True)
            except Exception as e:
                print(f"    {os.path.basename(p)}: ERROR {str(e)[:100]}", flush=True)
    return paths


def s2_mask(img):
    scl = img.select("SCL")
    good = scl.neq(3).And(scl.neq(8)).And(scl.neq(9)).And(scl.neq(10)).And(scl.neq(11))
    return img.updateMask(good)


def main():
    print("=== Sentinel-2, 10 m true colour, near-cloudless composite ===")
    s2 = (
        ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
        .filterBounds(lake_area)
        .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", 5))
        .map(s2_mask)
    )
    n = s2.size().getInfo()
    print(f"  scenes under 5% cloud: {n}")
    img = s2.median()
    # 8 x 8 keeps each tile near 35 MB, inside Earth Engine's 50 MB ceiling.
    tiled_download(img, ["B4", "B3", "B2"], 10, lake_area, "S2_10m_lake_composite", n=8)

    print("\n=== Landsat 5, 30 m, 2011 (nearest before the fieldwork) ===")

    def mask_l(img):
        qa = img.select("QA_PIXEL")
        return img.updateMask(qa.bitwiseAnd(1 << 3).eq(0).And(qa.bitwiseAnd(1 << 4).eq(0)))

    l5 = (
        ee.ImageCollection("LANDSAT/LT05/C02/T1_L2")
        .filterBounds(basin)
        .filterDate("2011-01-01", "2011-10-19")
        .map(mask_l)
    )
    print(f"  scenes: {l5.size().getInfo()}")
    tiled_download(l5.median(), ["SR_B3", "SR_B2", "SR_B1"], 30, basin, "L5_30m_2011", n=4)

    print("\n=== Landsat 8, 30 m, 2013 (nearest after the fieldwork) ===")
    l8 = (
        ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
        .filterBounds(basin)
        .filterDate("2013-03-25", "2013-12-31")
        .map(mask_l)
    )
    print(f"  scenes: {l8.size().getInfo()}")
    tiled_download(l8.median(), ["SR_B4", "SR_B3", "SR_B2"], 30, basin, "L8_30m_2013", n=4)

    print("\ndone; mosaic the tiles with 05.scripts/mosaic_highres.R")


if __name__ == "__main__":
    main()
