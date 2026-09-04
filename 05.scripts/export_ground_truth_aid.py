#!/usr/bin/env python3
"""Imagery to relocate the ground-truth photographs, on the days they were taken.

The problem this serves. 244 field photographs carry land-cover evidence and
almost no position: no GPS fix survives in any header, and the single located
frame set is a place name read out of a filename and geocoded at plus or minus
five kilometres. The author holds GPS records whose coordinate reference system
is uncertain, so what is needed is not another derived product but a picture of
the lake as it actually looked on each photography day, at a grain fine enough
to recognise a shoreline, a channel or a reed edge from the frame itself.

Three layers per date, all on one grid so they overlay exactly.

  rgb        Landsat true colour, red-green-blue, stretched for viewing rather
             than for radiometry. Composited over a window either side of the
             date because a single Landsat 7 pass in 2012 leaves scan-line gaps.
  mndwi      The modified normalised difference water index on the same pixels,
             and for 22 February 2012 the harmonic surface already validated
             against Water Observations from Space, since direct imagery covers
             only 22 per cent of the basin within eight days of that date.
  occurrence Global surface water occurrence, the share of the 1984 to 2021
             record in which each pixel held water. This is a locating aid only
             and does not enter the manuscript, per the closed decision that
             public global products do not appear in its text. Its value here is
             that it shows the permanent channel network and the lake's usual
             margin, which is what a photograph is recognised against.

Writes 03.outputs/TIF/ground_truth_aid/ as GeoTIFFs for QGIS, and matching
PNG quicklooks for the interactive tool.
"""

import io
import json
import os
import subprocess
import zipfile

import ee
import requests

ee.Initialize(project="murphys-deforisk")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "03.outputs", "TIF", "ground_truth_aid")
PNG = os.path.join(ROOT, "03.outputs", "PNG", "ground_truth_aid")
os.makedirs(OUT, exist_ok=True)
os.makedirs(PNG, exist_ok=True)

import geopandas as gpd
_lk = gpd.read_file(os.path.join(ROOT, "02.inputs", "SHP", "lakes_site.shp")).to_crs(4326)
LAKE = ee.Geometry(json.loads(json.dumps(_lk.geometry.union_all().__geo_interface__)))
AOI = LAKE.buffer(12000).bounds()
SCALE = 30

# Dates whose frames carry land-cover evidence worth locating. Archival,
# interior and specimen sessions are excluded because nothing in them can be
# matched to a pixel.
DATES = [
    ("2011-12-21", "dry vegetation and exposed lakebed, 2 frames"),
    ("2012-02-20", "open water at Phaloni and the Sombani inflow, 10 frames"),
    ("2012-02-22", "flooded vegetation and open water, 43 frames"),
    ("2012-03-29", "flooded vegetation, 3 frames"),
    ("2012-04-19", "open water, 5 frames"),
    ("2012-05-20", "flooded vegetation at a landing site, 4 frames"),
    ("2012-09-26", "open water, 4 frames"),
]
# Two windows, on purpose. The narrow one is closest to the photograph date and
# is the honest record of that fortnight; in February 2012 it holds only four
# Landsat 7 scenes, so its scan-line wedges are plainly visible. The wide one
# stacks enough passes that the wedges fill, because they move between passes,
# at the cost of no longer belonging to a single date. For recognising a
# shoreline, a channel or a reed edge from a photograph the wide one is the
# better tool; for reading water extent on the day it is the wrong one. Both are
# written and the filename says which is which.
WINDOW = 24
WIDE = 75
BANDS = {"LT05": ("SR_B3", "SR_B2", "SR_B1", "SR_B2", "SR_B5"),
         "LE07": ("SR_B3", "SR_B2", "SR_B1", "SR_B2", "SR_B5"),
         "LC08": ("SR_B4", "SR_B3", "SR_B2", "SR_B3", "SR_B6")}


def clear(img):
    qa = img.select("QA_PIXEL")
    bad = (qa.bitwiseAnd(1 << 1).Or(qa.bitwiseAnd(1 << 2)).Or(qa.bitwiseAnd(1 << 3))
           .Or(qa.bitwiseAnd(1 << 4)).Or(qa.bitwiseAnd(1 << 5)))
    return img.updateMask(bad.eq(0))


def scene_stack(date):
    d = ee.Date(date)
    out = None
    for sensor, (r, g, b, green, swir) in BANDS.items():
        col = (ee.ImageCollection(f"LANDSAT/{sensor}/C02/T1_L2").filterBounds(AOI)
               .filterDate(d.advance(-WINDOW, "day"), d.advance(WINDOW, "day")).map(clear))

        def prep(img, r=r, g=g, b=b, green=green, swir=swir):
            sr = img.select("SR_B.").multiply(0.0000275).add(-0.2)
            rgb = sr.select([r, g, b]).rename(["red", "green", "blue"])
            m = sr.normalizedDifference([green, swir]).rename("mndwi")
            return rgb.addBands(m)

        col = col.map(prep)
        out = col if out is None else out.merge(col)
    return out


def download(image, name, scale=SCALE, tiles=2):
    parts = []
    b = AOI.bounds().getInfo()["coordinates"][0]
    xs = [c[0] for c in b]; ys = [c[1] for c in b]
    x0, x1, y0, y1 = min(xs), max(xs), min(ys), max(ys)
    for i in range(tiles):
        for j in range(tiles):
            box = ee.Geometry.Rectangle([
                x0 + (x1 - x0) * i / tiles, y0 + (y1 - y0) * j / tiles,
                x0 + (x1 - x0) * (i + 1) / tiles, y0 + (y1 - y0) * (j + 1) / tiles])
            p = os.path.join(OUT, f"_{name}_{i}{j}.tif")
            if not (os.path.exists(p) and os.path.getsize(p) > 0):
                url = image.getDownloadURL({"scale": scale, "region": box.getInfo()["coordinates"],
                                            "filePerBand": False, "format": "ZIPPED_GEO_TIFF"})
                r = requests.get(url, timeout=900)
                if not r.ok:
                    raise RuntimeError(f"{name} {i}{j}: {r.status_code} {r.text[:200]}")
                with zipfile.ZipFile(io.BytesIO(r.content)) as z:
                    open(p, "wb").write(z.read(z.namelist()[0]))
            parts.append(p)
    out = os.path.join(OUT, f"{name}.tif")
    vrt = out.replace(".tif", ".vrt")
    subprocess.run(["gdalbuildvrt", vrt] + parts, check=True, capture_output=True)
    subprocess.run(["gdal_translate", "-of", "GTiff", "-co", "COMPRESS=DEFLATE",
                    "-co", "TILED=YES", vrt, out], check=True, capture_output=True)
    for f in parts + [vrt]:
        os.remove(f)
    print(f"  wrote {os.path.basename(out)} ({os.path.getsize(out)/1e6:.1f} MB)", flush=True)
    return out


def quicklook(tif, name, bands=(1, 2, 3), vmin=None, vmax=None, width=1400):
    """A viewable PNG on the same footprint, for the interactive tool."""
    png = os.path.join(PNG, f"{name}.png")
    cmd = ["gdal_translate", "-of", "PNG", "-ot", "Byte", "-outsize", str(width), "0"]
    for b in bands:
        cmd += ["-b", str(b)]
    if vmin is not None:
        cmd += ["-scale", str(vmin), str(vmax), "0", "255"]
    else:
        cmd += ["-scale"]
    subprocess.run(cmd + [tif, png], check=True, capture_output=True)
    print(f"  quicklook {os.path.basename(png)} ({os.path.getsize(png)/1e6:.1f} MB)", flush=True)
    return png


def main():
    # the static locating layer: where water has historically been
    occ = ee.Image("JRC/GSW1_4/GlobalSurfaceWater").select("occurrence").unmask(0).clip(AOI)
    t = download(occ.toByte(), "water_occurrence_1984_2021")
    quicklook(t, "water_occurrence", bands=(1,), vmin=0, vmax=100)

    harmonic = os.path.join(ROOT, "03.outputs", "TIF", "harmonic_2012",
                            "mndwi_harmonic_2012-02-22.tif")

    for date, what in DATES:
        print(f"{date}: {what}", flush=True)
        for win, tag in ((WINDOW, "near"), (WIDE, "gapfilled")):
            globals()["WINDOW"] = win
            col = scene_stack(date)
            n = col.size().getInfo()
            if n == 0:
                print(f"  {tag}: no Landsat within {win} days", flush=True); continue
            med = col.median()
            rgb = med.select(["red", "green", "blue"]).multiply(2200).clamp(0, 255).toByte()
            t = download(rgb, f"{date}_rgb_{tag}")
            quicklook(t, f"{date}_rgb_{tag}")
            mn = med.select("mndwi").multiply(10000).toInt16()
            t = download(mn, f"{date}_mndwi_{tag}")
            quicklook(t, f"{date}_mndwi_{tag}", bands=(1,), vmin=-4000, vmax=4000)
            print(f"  {tag}: {n} scenes within {win} days", flush=True)
        globals()["WINDOW"] = 24
        continue
        col = scene_stack(date)
        n = col.size().getInfo()
        if n == 0:
            print("  no Landsat within the window", flush=True); continue
        med = col.median()
        rgb = med.select(["red", "green", "blue"]).multiply(2200).clamp(0, 255).toByte()
        t = download(rgb, f"{date}_rgb")
        quicklook(t, f"{date}_rgb")
        mn = med.select("mndwi").multiply(10000).toInt16()
        t = download(mn, f"{date}_mndwi")
        quicklook(t, f"{date}_mndwi", bands=(1,), vmin=-4000, vmax=4000)
        print(f"  built from {n} Landsat scenes within {WINDOW} days", flush=True)

    print(f"\nGeoTIFFs in {OUT}\nquicklooks in {PNG}")
    if os.path.exists(harmonic):
        print(f"the validated harmonic surface for 22 February 2012 is already at\n  {harmonic}")


if __name__ == "__main__":
    main()
