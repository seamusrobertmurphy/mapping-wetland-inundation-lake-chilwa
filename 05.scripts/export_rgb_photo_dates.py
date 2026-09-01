#!/usr/bin/env python3
"""Export RGB GeoTIFFs matching the field-photograph dates.

What is actually available over this basin, verified by query rather than
assumed:

  2011-12 to 2012-12   Landsat 7 only. Thematic Mapper imaging ends
                       18 October 2011 and the Operational Land Imager opens
                       25 March 2013, so the Enhanced Thematic Mapper Plus is
                       the sole 30 m sensor across the fieldwork window. It
                       holds 41 scenes over the basin in 2012, and although
                       the scan-line corrector gaps hole any single scene, the
                       wedges move between passes, so 15 scenes over February
                       to June leave 18 km2 of the 8,752 km2 basin unobserved.
                       ASTER returns zero scenes for every photo month, because
                       it acquires on request and this basin was never tasked.
                       This script still writes MODIS at 500 m, which remains
                       the better choice for a single named day, since Landsat
                       7 near one date is thin, in that two scenes within
                       eight days of 22 February 2012 reach only 22 per cent
                       of the basin. For a
                       date rather than a day, see the harmonic fit in
                       05.scripts/harmonic_fit_2012.py, and for the coverage
                       measurement, 05.scripts/l7_2012_coverage_test.py.
  2014-11              Landsat 8 has six scenes, so that month is exported at
                       30 m as well.

MOD09GA (Terra) and MYD09GA (Aqua) are daily at 500 m and carry blue, green
and red, so true colour is possible. MOD09Q1/MOD09GQ reach 250 m but hold only
red and near-infrared, so they cannot make a true-colour image; they are not
used here.

Per-date images are attempted first, since a photograph belongs to a day. Where
cloud defeats the day, the monthly median stands in. Both are written, and the
clear-pixel fraction is reported so a scene is never presented as clean when it
is not.

Writes 03.outputs/TIF/RGB_photo_dates/
"""

import csv
import io
import os
import zipfile

import ee
import requests

ee.Initialize(project="murphys-deforisk")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "03.outputs", "TIF", "RGB_photo_dates")
os.makedirs(OUT, exist_ok=True)

import json

gj = json.load(open(os.path.join(ROOT, "03.outputs", "JSON", "chilwa_basin.geojson")))
geom = gj["features"][0]["geometry"] if gj.get("type") == "FeatureCollection" else gj
aoi = ee.Geometry(geom)
region = aoi.bounds()

MODIS_RGB = ["sur_refl_b01", "sur_refl_b04", "sur_refl_b03"]  # red, green, blue


def modis_clear(img):
    """Mask cloud and cloud shadow from the 500 m state flags."""
    s = img.select("state_1km")
    clear = s.bitwiseAnd(3).eq(0).And(s.bitwiseAnd(1 << 2).eq(0))
    return img.updateMask(clear)


def fetch(image, name, scale, bands):
    """Download a clipped GeoTIFF. Returns the written path or None."""
    try:
        url = image.select(bands).clip(aoi).getDownloadURL(
            {"scale": scale, "region": region, "format": "GEO_TIFF", "crs": "EPSG:4326"}
        )
        r = requests.get(url, timeout=600)
        if r.status_code != 200 or len(r.content) < 2000:
            print(f"    {name}: download failed ({r.status_code}, {len(r.content)} bytes)")
            return None
        path = os.path.join(OUT, name + ".tif")
        if r.content[:2] == b"PK":  # zipped multiband
            z = zipfile.ZipFile(io.BytesIO(r.content))
            z.extractall(OUT)
            print(f"    {name}: extracted {len(z.namelist())} band files")
            return OUT
        with open(path, "wb") as fh:
            fh.write(r.content)
        print(f"    {name}: {len(r.content)/1e6:.1f} MB")
        return path
    except Exception as e:
        print(f"    {name}: ERROR {str(e)[:110]}")
        return None


def clear_fraction(img, scale):
    try:
        v = (
            img.select(MODIS_RGB[0])
            .mask()
            .reduceRegion(ee.Reducer.mean(), aoi, scale, maxPixels=int(1e9))
            .get(MODIS_RGB[0])
            .getInfo()
        )
        return float(v) if v is not None else 0.0
    except Exception:
        return 0.0


def main():
    meta = list(csv.DictReader(open(os.path.join(ROOT, "03.outputs", "CSV", "field_photo_metadata.csv"))))
    dates = sorted({m["date"] for m in meta if m["date"].strip()})
    months = sorted({d[:7] for d in dates})
    counts = {}
    for m in meta:
        if m["date"].strip():
            counts[m["date"]] = counts.get(m["date"], 0) + 1

    log = []

    print(f"{len(dates)} distinct photo dates across {len(months)} months\n")

    print("=== monthly MODIS median, 500 m true colour ===")
    for mo in months:
        y, mm = int(mo[:4]), int(mo[5:7])
        start = f"{y}-{mm:02d}-01"
        end = f"{y+1}-01-01" if mm == 12 else f"{y}-{mm:02d}".rsplit("-", 1)[0] + f"-{mm+1:02d}-01"
        col = (
            ee.ImageCollection("MODIS/061/MOD09GA")
            .filterBounds(aoi)
            .filterDate(start, end)
            .map(modis_clear)
        )
        n = col.size().getInfo()
        img = col.median()
        cf = clear_fraction(img, 500)
        print(f"  {mo}  n={n:3d}  clear={cf:.2f}")
        fetch(img, f"MODIS500_monthly_{mo}", 500, MODIS_RGB)
        log.append({"target": mo, "kind": "monthly median", "sensor": "MOD09GA", "scale_m": 500,
                    "n_images": n, "clear_fraction": round(cf, 3), "n_photos": sum(counts[d] for d in dates if d.startswith(mo))})

    print("\n=== per-date MODIS, 500 m true colour ===")
    for d in dates:
        col = (
            ee.ImageCollection("MODIS/061/MOD09GA")
            .filterBounds(aoi)
            .filterDate(d, ee.Date(d).advance(1, "day"))
            .map(modis_clear)
        )
        n = col.size().getInfo()
        if n == 0:
            print(f"  {d}: no image")
            continue
        img = col.mosaic()
        cf = clear_fraction(img, 500)
        flag = "" if cf > 0.5 else "  (mostly cloud)"
        print(f"  {d}  photos={counts[d]:3d}  clear={cf:.2f}{flag}")
        fetch(img, f"MODIS500_date_{d}", 500, MODIS_RGB)
        log.append({"target": d, "kind": "single date", "sensor": "MOD09GA", "scale_m": 500,
                    "n_images": n, "clear_fraction": round(cf, 3), "n_photos": counts[d]})

    print("\n=== Landsat 8, 30 m true colour, November 2014 ===")
    l8 = (
        ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
        .filterBounds(aoi)
        .filterDate("2014-11-01", "2014-12-01")
    )
    n = l8.size().getInfo()
    if n:
        def mask_l8(img):
            qa = img.select("QA_PIXEL")
            m = qa.bitwiseAnd(1 << 3).eq(0).And(qa.bitwiseAnd(1 << 4).eq(0))
            return img.updateMask(m)

        img = l8.map(mask_l8).median()
        print(f"  2014-11  n={n}")
        fetch(img, "L8_30m_monthly_2014-11", 30, ["SR_B4", "SR_B3", "SR_B2"])
        log.append({"target": "2014-11", "kind": "monthly median", "sensor": "L8 OLI", "scale_m": 30,
                    "n_images": n, "clear_fraction": "", "n_photos": sum(counts[d] for d in dates if d.startswith("2014-11"))})

    with open(os.path.join(OUT, "manifest.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(log[0].keys()))
        w.writeheader()
        w.writerows(log)
    print(f"\nmanifest written: {len(log)} entries")


if __name__ == "__main__":
    main()
