#!/usr/bin/env python3
"""ASTER VNIR at 15 m over the Lake Chilwa basin.

ASTER is the finest free optical sensor that ever imaged this basin. It
acquires on request rather than systematically, and the request record decides
everything: across all seven ASTER collections in NASA CMR (L1A, L1B, L1T and
the L2 surface products including AST_09), the entire 2012 holding over this
basin is 15, 22 and 31 July and 19 October. Nothing between February and June,
so the fieldwork window itself cannot be covered at 15 m by any product.

What does exist is cloud-free coverage in the neighbouring years, seasonally
matched to the photo sessions: 12 May 2011 at zero per cent cloud and
17 April 2011 at four per cent, against photo sessions through April and May
2012. For placing ground-truth points, where shorelines and villages are
stable between years, that is the most useful imagery available anywhere.

VNIR carries green (B01), red (B02) and near-infrared (B3N) only. There is no
blue band, so the composite is the standard false-colour NIR-red-green, not
true colour. SWIR failed in January 2009 and is unusable for any of these
dates.

Writes 03.outputs/TIF/RGB_highres/ASTER15m_<date>_t*.tif (mosaic with
05.scripts/mosaic_highres.R)
"""

import json
import os

import ee
import requests

ee.Initialize(project="murphys-deforisk")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "03.outputs", "TIF", "RGB_highres")
os.makedirs(OUT, exist_ok=True)

gj = json.load(open(os.path.join(ROOT, "03.outputs", "JSON", "chilwa_basin.geojson")))
geom = gj["features"][0]["geometry"] if gj.get("type") == "FeatureCollection" else gj
aoi = ee.Geometry(geom)

BANDS = ["B3N", "B02", "B01"]  # NIR, red, green -> false colour
SCALE = 15

# Dates whose footprint contains the LAKE, not merely the basin. The
# distinction is decisive: the cloud-free May 2011 pass covers the eastern
# basin and reaches 0.1 per cent of lake pixels. 6 July 2011 (1 per cent
# cloud) and 12 August 2013 (0 per cent) are the nearest near-cloud-free
# scenes over the water, bracketing the fieldwork.
DATES = ["2011-07-06", "2013-08-12", "2013-10-22"]


def tiles(image, name, n=6):
    b = aoi.bounds().getInfo()["coordinates"][0]
    xs = [p[0] for p in b]
    ys = [p[1] for p in b]
    x0, x1, y0, y1 = min(xs), max(xs), min(ys), max(ys)
    ok = 0
    for i in range(n):
        for j in range(n):
            reg = ee.Geometry.Rectangle([
                x0 + (x1 - x0) * i / n, y0 + (y1 - y0) * j / n,
                x0 + (x1 - x0) * (i + 1) / n, y0 + (y1 - y0) * (j + 1) / n])
            p = os.path.join(OUT, f"{name}_t{i}{j}.tif")
            if os.path.exists(p) and os.path.getsize(p) > 2000:
                ok += 1
                continue
            try:
                url = image.clip(reg).getDownloadURL(
                    {"scale": SCALE, "region": reg, "format": "GEO_TIFF", "crs": "EPSG:4326"})
                r = requests.get(url, timeout=900)
                if r.status_code == 200 and len(r.content) > 2000:
                    with open(p, "wb") as fh:
                        fh.write(r.content)
                    ok += 1
                    print(f"    t{i}{j}: {len(r.content)/1e6:.1f} MB", flush=True)
                else:
                    print(f"    t{i}{j}: empty/HTTP {r.status_code}", flush=True)
            except Exception as e:
                print(f"    t{i}{j}: {str(e)[:90]}", flush=True)
    return ok


for d in DATES:
    col = (ee.ImageCollection("ASTER/AST_L1T_003")
           .filterBounds(aoi)
           .filterDate(d, ee.Date(d).advance(1, "day"))
           .filter(ee.Filter.lt("CLOUDCOVER", 20)))
    n = col.size().getInfo()
    if n == 0:
        print(f"{d}: no usable scene")
        continue
    print(f"{d}: {n} scenes, mosaicking at {SCALE} m")
    img = col.select(BANDS).mosaic()
    got = tiles(img, f"ASTER15m_{d}")
    print(f"  tiles written: {got}")

print("done; run 05.scripts/mosaic_highres.R to join")
