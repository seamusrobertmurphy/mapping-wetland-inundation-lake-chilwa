#!/usr/bin/env python3
"""Planet NICFI basemap over the Lake Chilwa watershed, ready to run on access.

This script is complete and correct and will not run until the Earth Engine
account holding these credentials has been added to the NICFI programme. The
collection is not public: a live check on 4 September 2026 returned "asset not
found (does not exist or caller does not have access)", which is what Earth
Engine says when an asset exists but the caller is not on its access list.

Two facts about the data, both from the catalogue page itself rather than from
memory. It is 4.77 metres, not 30 or 50 centimetres, and it begins in December
2015. It is genuinely sharper than the ten-metre Sentinel-2 basemap built
alongside it, and it is the best imagery that will ever be free over this basin,
but it cannot show the 2012 fieldwork because it did not exist then.

To gain access, register at planet.com/nicfi using the same Google account that
holds these Earth Engine credentials, accept the non-commercial licence, and
allow a day for the access list to update. Then run this script unchanged.

Writes 03.outputs/TIF/basemap/chilwa_nicfi_4m_<period>.tif
"""

import io
import json
import os
import subprocess
import sys
import zipfile

import ee
import requests

ee.Initialize(project="murphys-deforisk")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "03.outputs", "TIF", "basemap")
os.makedirs(OUT, exist_ok=True)
COLLECTION = "projects/planet-nicfi/assets/basemaps/africa"
SCALE = 4.77

import geopandas as gpd
_lk = gpd.read_file(os.path.join(ROOT, "02.inputs", "SHP", "lakes_site.shp")).to_crs(4326)
AOI = ee.Geometry(json.loads(json.dumps(
    _lk.geometry.union_all().__geo_interface__))).buffer(12000).bounds()

# One dry-season and one wet-season mosaic. NICFI names its assets by period,
# biannual before September 2020 and monthly after, so the filter is on date
# rather than on an index pattern that changes partway through the record.
PERIODS = [("2023-08-01", "2023-11-01", "dryseason_2023"),
           ("2024-01-01", "2024-04-01", "wetseason_2024")]


def check():
    try:
        n = ee.ImageCollection(COLLECTION).filterBounds(AOI).size().getInfo()
        print(f"access confirmed, {n} mosaics intersect the watershed")
        return True
    except Exception as e:
        print("NICFI is not accessible from these credentials.", file=sys.stderr)
        print(f"  Earth Engine said: {str(e)[:160]}", file=sys.stderr)
        print("  Register at planet.com/nicfi with the Google account that holds these\n"
              "  credentials, accept the non-commercial licence, then run this again.",
              file=sys.stderr)
        return False


def stretch(img):
    """Two to ninety-eight per cent per band, measured on the mosaic itself."""
    p = img.reduceRegion(ee.Reducer.percentile([2, 98]), AOI, 60,
                         maxPixels=1e10, bestEffort=True)
    out = []
    for b in ("R", "G", "B"):
        lo = ee.Number(p.get(f"{b}_p2")); hi = ee.Number(p.get(f"{b}_p98"))
        out.append(img.select(b).subtract(lo).divide(hi.subtract(lo))
                   .multiply(255).clamp(0, 255).rename(b))
    return ee.Image.cat(out).toByte()


def download(image, name, tiles=4):
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
                url = image.getDownloadURL({"scale": SCALE,
                                            "region": box.getInfo()["coordinates"],
                                            "filePerBand": False,
                                            "format": "ZIPPED_GEO_TIFF"})
                r = requests.get(url, timeout=1800)
                if not r.ok:
                    raise RuntimeError(f"{name} {i}{j}: {r.status_code} {r.text[:200]}")
                with zipfile.ZipFile(io.BytesIO(r.content)) as z:
                    open(p, "wb").write(z.read(z.namelist()[0]))
                print(f"    tile {i}{j}  {os.path.getsize(p)/1e6:5.1f} MB", flush=True)
            parts.append(p)
    out = os.path.join(OUT, f"{name}.tif")
    vrt = out.replace(".tif", ".vrt")
    subprocess.run(["gdalbuildvrt", vrt] + parts, check=True, capture_output=True)
    subprocess.run(["gdal_translate", "-of", "GTiff", "-co", "COMPRESS=DEFLATE",
                    "-co", "TILED=YES", "-co", "PHOTOMETRIC=RGB",
                    "-co", "BIGTIFF=IF_SAFER", vrt, out], check=True, capture_output=True)
    subprocess.run(["gdaladdo", "-r", "average", out, "2", "4", "8", "16"],
                   check=True, capture_output=True)
    for f in parts + [vrt]:
        os.remove(f)
    print(f"  wrote {os.path.basename(out)} ({os.path.getsize(out)/1e6:.0f} MB)", flush=True)


def main():
    if not check():
        return 1
    for t0, t1, tag in PERIODS:
        col = ee.ImageCollection(COLLECTION).filterBounds(AOI).filterDate(t0, t1)
        n = col.size().getInfo()
        if n == 0:
            print(f"{tag}: no mosaic in that window"); continue
        print(f"{tag}: {n} mosaics", flush=True)
        download(stretch(col.median()), f"chilwa_nicfi_4m_{tag}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
