#!/usr/bin/env python3
"""A readable basemap of the Lake Chilwa basin, for recognising ground from photographs.

Why this exists. The near-date Landsat composites are the honest record of the
fortnight a photograph was taken, and they are close to unreadable: four
Landsat 7 scenes in February 2012, scan-line wedges through half the frame,
cloud holes through the rest, at thirty metres. That is fine for measuring water
and useless for recognising a hill, a village, a road or the mouth of a channel.

Those are two different jobs and they need two different images. This one
abandons the date entirely. It composites Sentinel-2 across nine dry seasons at
ten metres, so every pixel is chosen from dozens of clear looks and the result
has no gaps at all. It cannot tell you where the water was in 2012. It can tell
you what the ground looks like, and the ground has not moved.

Two seasons are built. The dry-season image, August to October, shows the
exposed lakebed, the channel network and the settlement pattern at their
clearest. The wet-season image, January to March, shows the flood extent the
basin reaches in a normal year, which is what the reed belt sits inside.

The stretch is computed from the image rather than assumed. An earlier version
multiplied reflectance by a fixed constant, which crushed the water to black and
washed the land out. Here the second and ninety-eighth percentiles of each band
are measured over the scene and mapped to the byte range, separately per band,
which is what makes a satellite image look like a photograph.

Writes 03.outputs/TIF/basemap/ at 10 m, ready to drag into QGIS.
"""

import io
import json
import os
import subprocess
import time
import zipfile

import ee
import requests

ee.Initialize(project="murphys-deforisk")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "03.outputs", "TIF", "basemap")
os.makedirs(OUT, exist_ok=True)

import geopandas as gpd
_lk = gpd.read_file(os.path.join(ROOT, "02.inputs", "SHP", "lakes_site.shp")).to_crs(4326)
LAKE = ee.Geometry(json.loads(json.dumps(_lk.geometry.union_all().__geo_interface__)))
AOI = LAKE.buffer(12000).bounds()
SCALE = 10

SEASONS = {
    "dry": ((8, 10), "August to October, the lakebed and channels at their clearest"),
    "wet": ((1, 3), "January to March, the flood extent of a normal year"),
}
YEARS = (2016, 2024)


def masked(img):
    """Sentinel-2 scene classification: drop cloud, shadow, cirrus and snow."""
    scl = img.select("SCL")
    good = (scl.eq(4).Or(scl.eq(5)).Or(scl.eq(6)).Or(scl.eq(7))
            .Or(scl.eq(2)).Or(scl.eq(11)))
    return img.updateMask(good).divide(10000)


def composite(months):
    m0, m1 = months
    col = (ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED").filterBounds(AOI)
           .filterDate(f"{YEARS[0]}-01-01", f"{YEARS[1] + 1}-01-01")
           .filter(ee.Filter.calendarRange(m0, m1, "month"))
           .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", 35))
           .map(masked))
    return col.select(["B4", "B3", "B2"]).median().rename(["red", "green", "blue"]), col.size()


def stretch(img):
    """Two to ninety-eight per cent per band, measured on the image itself."""
    p = img.reduceRegion(ee.Reducer.percentile([2, 98]), AOI, 120,
                         maxPixels=1e10, bestEffort=True)
    out = []
    for b in ("red", "green", "blue"):
        lo = ee.Number(p.get(f"{b}_p2")); hi = ee.Number(p.get(f"{b}_p98"))
        out.append(img.select(b).subtract(lo).divide(hi.subtract(lo))
                   .multiply(255).clamp(0, 255).rename(b))
    return ee.Image.cat(out).toByte(), p


def download(image, name, tiles=3):
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
    print(f"  wrote {os.path.basename(out)} ({os.path.getsize(out)/1e6:.0f} MB, "
          f"with overviews so QGIS pans smoothly)", flush=True)


ASSET_ROOT = "projects/murphys-deforisk/assets"


def to_asset(image, name):
    """Store the mosaic first, then read tiles from the stored copy.

    A median over six hundred Sentinel-2 scenes is recomputed for every tile
    the interactive download requests, which exceeds the user memory limit. The
    batch exporter has no such limit, so the composite is computed once, written
    to an asset, and the tiles are then pulled from stored pixels where no
    computation happens at all."""
    aid = f"{ASSET_ROOT}/{name}"
    try:
        ee.data.getAsset(aid)
        print(f"  asset already exists: {aid}", flush=True)
        return ee.Image(aid)
    except ee.EEException:
        pass
    task = ee.batch.Export.image.toAsset(
        image=image, description=name[:99], assetId=aid,
        region=AOI, scale=SCALE, crs="EPSG:4326", maxPixels=1e10)
    task.start()
    print(f"  asset export started ({task.id}); this is the slow step", flush=True)
    while task.active():
        time.sleep(45)
    st = task.status()
    if st["state"] != "COMPLETED":
        raise SystemExit(f"asset export failed: {st.get('error_message', st)}")
    print("  asset export complete", flush=True)
    return ee.Image(aid)


def main():
    for key, (months, what) in SEASONS.items():
        print(f"{key} season: {what}", flush=True)
        img, n = composite(months)
        st, p = stretch(img)
        info = ee.Dictionary({"n": n, "p": p}).getInfo()
        print(f"  {info['n']} Sentinel-2 scenes under 35 per cent cloud, "
              f"{YEARS[0]} to {YEARS[1]}", flush=True)
        print("  stretch: " + ", ".join(
            f"{b} {info['p'][b+'_p2']:.3f} to {info['p'][b+'_p98']:.3f}"
            for b in ("red", "green", "blue")), flush=True)
        name = f"chilwa_s2_10m_{key}season_{YEARS[0]}_{YEARS[1]}"
        stored = to_asset(st, name)
        download(stored, name)


if __name__ == "__main__":
    main()
