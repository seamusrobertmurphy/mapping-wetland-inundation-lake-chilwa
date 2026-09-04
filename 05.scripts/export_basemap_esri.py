#!/usr/bin/env python3
"""A sub-five-metre basemap of the Lake Chilwa basin from Esri World Imagery.

Why this and not a satellite collection. Nothing under five metres is available
free over this basin as downloadable imagery: Maxar is commercial, SkySat, NAIP
and ALOS AVNIR-2 return zero scenes here, and NICFI at 4.77 m needs a
registration. What does exist is a web tile service. Esri World Imagery serves
roughly 0.6 metres over Lake Chilwa, dated 2025, and a tile fetched at zoom 18
over Phaloni resolves individual trees, field boundaries and footpaths.

This is a positioning aid and nothing else. It has no calibrated radiometry, no
stated acquisition date per pixel, and no band values that mean anything
physical. It must never enter a water index, a classification or any measured
quantity. Its only job is to let a photograph be matched to a place.

Zoom 16 gives 2.3 metres at this latitude, which is inside the five-metre
requirement and keeps the mosaic to a manageable size over a basin this large.
Zoom 17 doubles the detail and quadruples the tile count; pass --zoom 17 if the
shoreline needs it, and expect roughly four times the download.

Licence. Esri World Imagery is free to view and is standard practice as a
backdrop in QGIS. Bulk download for redistribution is not permitted. This builds
a local copy for one researcher's positioning work, which is the same use as the
live layer, and it is not to be published or shared as imagery.

Writes 03.outputs/TIF/basemap/chilwa_esri_<res>m_z<zoom>.tif
"""

import argparse
import math
import os
import subprocess
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor

import geopandas as gpd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "03.outputs", "TIF", "basemap")
TMP = os.path.join(OUT, "_tiles")
os.makedirs(TMP, exist_ok=True)
URL = ("https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery"
       "/MapServer/tile/{z}/{y}/{x}")
UA = {"User-Agent": "Mozilla/5.0 (research positioning aid)"}


def deg2tile(lat, lon, z):
    n = 2 ** z
    x = (lon + 180.0) / 360.0 * n
    r = math.radians(lat)
    y = (1.0 - math.log(math.tan(r) + 1 / math.cos(r)) / math.pi) / 2.0 * n
    return x, y


def tile2deg(x, y, z):
    n = 2 ** z
    lon = x / n * 360.0 - 180.0
    lat = math.degrees(math.atan(math.sinh(math.pi * (1 - 2 * y / n))))
    return lat, lon


def fetch(args):
    z, x, y = args
    p = os.path.join(TMP, f"{z}_{x}_{y}.jpg")
    if os.path.exists(p) and os.path.getsize(p) > 500:
        return p
    for attempt in range(3):
        try:
            d = urllib.request.urlopen(
                urllib.request.Request(URL.format(z=z, x=x, y=y), headers=UA),
                timeout=40).read()
            if len(d) > 500:
                open(p, "wb").write(d)
                return p
        except Exception:
            time.sleep(1.5 * (attempt + 1))
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--zoom", type=int, default=16)
    ap.add_argument("--buffer-km", type=float, default=12.0)
    a = ap.parse_args()
    z = a.zoom

    lk = gpd.read_file(os.path.join(ROOT, "02.inputs", "SHP", "lakes_site.shp")).to_crs(4326)
    geom = lk.to_crs(32736).buffer(a.buffer_km * 1000).to_crs(4326).union_all()
    w, s, e, n = geom.bounds

    x0, y0 = deg2tile(n, w, z)
    x1, y1 = deg2tile(s, e, z)
    x0, x1 = int(math.floor(x0)), int(math.ceil(x1))
    y0, y1 = int(math.floor(y0)), int(math.ceil(y1))
    res = 156543.03 * math.cos(math.radians((n + s) / 2)) / (2 ** z)
    jobs = [(z, x, y) for x in range(x0, x1) for y in range(y0, y1)]
    print(f"zoom {z}, {res:.2f} m per pixel, {x1-x0} by {y1-y0} tiles = {len(jobs):,}", flush=True)
    print(f"extent {w:.4f} to {e:.4f} east, {s:.4f} to {n:.4f} north", flush=True)

    got, miss = [], 0
    with ThreadPoolExecutor(max_workers=12) as ex:
        for i, p in enumerate(ex.map(fetch, jobs)):
            if p:
                got.append(p)
            else:
                miss += 1
            if (i + 1) % 500 == 0:
                print(f"  {i+1:,}/{len(jobs):,} fetched, {miss} missing", flush=True)
    print(f"  {len(got):,} tiles on disk, {miss} could not be fetched", flush=True)

    # georeference each tile, then build one mosaic
    vrts = []
    for p in got:
        zz, xx, yy = (int(v) for v in os.path.basename(p)[:-4].split("_"))
        top, left = tile2deg(xx, yy, zz)
        bot, right = tile2deg(xx + 1, yy + 1, zz)
        v = p.replace(".jpg", ".vrt")
        if not os.path.exists(v):
            subprocess.run(["gdal_translate", "-of", "VRT", "-a_srs", "EPSG:4326",
                            "-a_ullr", str(left), str(top), str(right), str(bot), p, v],
                           check=True, capture_output=True)
        vrts.append(v)

    lst = os.path.join(TMP, "tiles.txt")
    open(lst, "w").write("\n".join(vrts))
    big = os.path.join(TMP, "mosaic.vrt")
    subprocess.run(["gdalbuildvrt", "-input_file_list", lst, big], check=True, capture_output=True)
    out = os.path.join(OUT, f"chilwa_esri_{res:.1f}m_z{z}.tif")
    subprocess.run(["gdal_translate", "-of", "GTiff", "-co", "COMPRESS=JPEG",
                    "-co", "PHOTOMETRIC=YCBCR", "-co", "JPEG_QUALITY=88",
                    "-co", "TILED=YES", "-co", "BIGTIFF=IF_SAFER", big, out],
                   check=True, capture_output=True)
    subprocess.run(["gdaladdo", "-r", "average", out, "2", "4", "8", "16", "32"],
                   check=True, capture_output=True)
    print(f"\nwrote {out}  ({os.path.getsize(out)/1e6:.0f} MB at {res:.2f} m)", flush=True)


if __name__ == "__main__":
    sys.exit(main())
