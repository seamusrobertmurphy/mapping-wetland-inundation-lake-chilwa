#!/usr/bin/env python3
"""Landsat 7 true colour sharpened to 15 m, on the field photograph dates.

Landsat 7 is the only optical sensor imaging the basin between October 2011 and
March 2013, and it carries a panchromatic band at 15 m beside its 30 m colour
bands. Sharpening lays the brightness of the 15 m band over the colour of the
30 m bands, which halves the pixel size while keeping the hue. That is the
finest imagery taken inside the fieldwork window, and its job is to let a
photograph be placed against a shoreline, channel or reed edge seen on the same
fortnight.

The method is high-pass modulation. Each surface-reflectance colour band is
multiplied by the panchromatic band divided by a 45 m mean of itself. That ratio
averages one, so the 30 m colours are kept and only the 15 m edges are added. A
hue, saturation and value substitution was tried first and rejected, because the
panchromatic band spans 0.52 to 0.90 micrometres, takes in near infrared, and
turned vegetation pink once it replaced the brightness of the colour image.

Dates, area and the two windows match 05.scripts/export_ground_truth_aid.py so
the layers overlay the 30 m versions. The near composite is within 24 days of
the photograph and keeps its scan-line wedges; the gap-filled composite is
within 75 days and fills them, because the wedges move between passes. This is
a placement aid only and enters no measured quantity.

Writes 03.outputs/TIF/pansharpened_15m/ and quicklooks to
03.outputs/PNG/pansharpened_15m/
"""

import io
import json
import os
import subprocess
import zipfile

import ee
import geopandas as gpd
import requests

ee.Initialize(project="murphys-deforisk")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "03.outputs", "TIF", "pansharpened_15m")
PNG = os.path.join(ROOT, "03.outputs", "PNG", "pansharpened_15m")
os.makedirs(OUT, exist_ok=True)
os.makedirs(PNG, exist_ok=True)

_lk = gpd.read_file(os.path.join(ROOT, "02.inputs", "SHP", "lakes_site.shp")).to_crs(4326)
LAKE = ee.Geometry(json.loads(json.dumps(_lk.geometry.union_all().__geo_interface__)))
AOI = LAKE.buffer(12000).bounds()
SCALE = 15

DATES = [
    ("2011-12-21", "dry vegetation and exposed lakebed, 2 frames"),
    ("2012-02-20", "open water at Phaloni and the Sombani inflow, 10 frames"),
    ("2012-02-22", "flooded vegetation and open water, 43 frames"),
    ("2012-03-29", "flooded vegetation, 3 frames"),
    ("2012-04-19", "open water, 5 frames"),
    ("2012-05-20", "flooded vegetation at a landing site, 4 frames"),
    ("2012-09-26", "open water, 4 frames"),
]
WINDOWS = ((24, "near"), (75, "gapfilled"))


def clear(img):
    qa = img.select("QA_PIXEL")
    bad = (qa.bitwiseAnd(1 << 1).Or(qa.bitwiseAnd(1 << 2)).Or(qa.bitwiseAnd(1 << 3))
           .Or(qa.bitwiseAnd(1 << 4)).Or(qa.bitwiseAnd(1 << 5)))
    return img.updateMask(bad.eq(0)).resample("bilinear")


def composite(date, days):
    d = ee.Date(date)
    # Landsat 7 ETM+ Collection 2 Tier 1 Level 2 surface reflectance, Earth Engine catalogue
    # https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LE07_C02_T1_L2
    # USGS landing page https://doi.org/10.5066/P9TU80IG ; no inputs manifest row exists for it.
    sr = (ee.ImageCollection("LANDSAT/LE07/C02/T1_L2").filterBounds(AOI)
          .filterDate(d.advance(-days, "day"), d.advance(days, "day")).map(clear))
    # Landsat 7 ETM+ Collection 2 Tier 1 top-of-atmosphere reflectance, for the 15 m band B8
    # https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LE07_C02_T1_TOA
    # USGS landing page https://doi.org/10.5066/P9TU80IG ; no inputs manifest row exists for it.
    toa = (ee.ImageCollection("LANDSAT/LE07/C02/T1_TOA").filterBounds(AOI)
           .filterDate(d.advance(-days, "day"), d.advance(days, "day")).map(clear))
    rgb = (sr.select(["SR_B3", "SR_B2", "SR_B1"], ["red", "green", "blue"]).median()
           .multiply(0.0000275).add(-0.2))
    return sr, rgb, toa.select("B8").median()


def sharpen(rgb, pan):
    detail = pan.divide(pan.focal_mean(1, "square", "pixels"))
    sharp = rgb.multiply(detail).updateMask(rgb.select("red").mask())
    cover = rgb.select("red").mask().gt(0).reduceRegion(
        ee.Reducer.mean(), AOI, 150, bestEffort=True).getInfo()["red"]
    return sharp.multiply(2200).clamp(1, 255).toByte(), cover


def download(image, name, tiles=3):
    out = os.path.join(OUT, f"{name}.tif")
    if os.path.exists(out):
        print(f"  exists {os.path.basename(out)}", flush=True)
        return out
    b = AOI.bounds().getInfo()["coordinates"][0]
    xs = [c[0] for c in b]; ys = [c[1] for c in b]
    x0, x1, y0, y1 = min(xs), max(xs), min(ys), max(ys)
    parts = []
    for i in range(tiles):
        for j in range(tiles):
            box = ee.Geometry.Rectangle([
                x0 + (x1 - x0) * i / tiles, y0 + (y1 - y0) * j / tiles,
                x0 + (x1 - x0) * (i + 1) / tiles, y0 + (y1 - y0) * (j + 1) / tiles])
            p = os.path.join(OUT, f"_{name}_{i}{j}.tif")
            if not (os.path.exists(p) and os.path.getsize(p) > 0):
                url = image.getDownloadURL({"scale": SCALE, "crs": "EPSG:4326",
                                            "region": box.getInfo()["coordinates"],
                                            "filePerBand": False, "format": "ZIPPED_GEO_TIFF"})
                r = requests.get(url, timeout=900)
                if not r.ok:
                    raise RuntimeError(f"{name} {i}{j}: {r.status_code} {r.text[:200]}")
                with zipfile.ZipFile(io.BytesIO(r.content)) as z:
                    open(p, "wb").write(z.read(z.namelist()[0]))
            parts.append(p)
    vrt = out.replace(".tif", ".vrt")
    subprocess.run(["gdalbuildvrt", vrt] + parts, check=True, capture_output=True)
    subprocess.run(["gdal_translate", "-of", "GTiff", "-co", "COMPRESS=DEFLATE", "-co", "TILED=YES",
                    "-a_nodata", "0", vrt, out], check=True, capture_output=True)
    for f in parts + [vrt]:
        os.remove(f)
    print(f"  wrote {os.path.basename(out)} ({os.path.getsize(out)/1e6:.1f} MB)", flush=True)
    return out


def quicklook(tif, name, width=1400):
    png = os.path.join(PNG, f"{name}.png")
    subprocess.run(["gdal_translate", "-of", "PNG", "-outsize", str(width), "0", tif, png],
                   check=True, capture_output=True)


def main():
    for date, what in DATES:
        print(f"{date}: {what}", flush=True)
        for days, tag in WINDOWS:
            col, rgb, pan = composite(date, days)
            n = col.size().getInfo()
            if n == 0:
                print(f"  {tag}: no Landsat 7 within {days} days", flush=True)
                continue
            img, cover = sharpen(rgb, pan)
            t = download(img, f"{date}_rgb15_{tag}")
            quicklook(t, f"{date}_rgb15_{tag}")
            print(f"  {tag}: {n} scenes within {days} days, {cover*100:.1f} per cent of the area observed",
                  flush=True)


if __name__ == "__main__":
    main()
