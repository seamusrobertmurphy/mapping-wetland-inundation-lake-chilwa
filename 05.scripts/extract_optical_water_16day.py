#!/usr/bin/env python3
"""Observed water area over the Lake Chilwa basin on the sixteen-day grid.

One row per sixteen-day step per sensor family, holding the water area actually
seen in that window and how much of the basin was clear enough to see. Nothing
is modelled here and nothing is filled in: a step with no clear observation
returns nothing, and the fusion downstream is what carries the record across
those gaps. Keeping observation and model apart is the point, because a
harmonic curve read every sixteen days is interpolation wearing the costume of
data, and only a file that says which steps were observed can tell the two
apart later.

Landsat is the anchor. Its thirty-metre grid resolves the reed fringe where the
water boundary actually sits, and every coarser stream is calibrated against it
rather than against another coarse stream. Water is taken as a modified
normalised difference water index above zero, the threshold the 2012 work
validated against Water Observations from Space.

MODIS is the dense stream from February 2000, at five hundred metres from Terra
and Aqua together. It carries the same index, because the MODIS surface
reflectance product has both a green and a shortwave infrared band, so the
water rule is identical and only the grain changes.

Each row also records the clear area behind it. A step where cloud left a
quarter of the region visible is a much weaker observation than one where all
of it was clear, and the fusion needs that difference as a variance rather than
as a footnote.

Clear cover is judged over the lake and a fifteen-kilometre belt around it, not
over the whole basin. Cloud sitting on the hills forty kilometres away does not
affect a measurement of the lake, and judging on the basin discarded steps that
were in fact fully usable: measured over the basin only 31 of 366 steps before
2000 reached eighty per cent clear, which for the era with no other dense sensor
is the difference between a record and a gap. The basin figure is kept in its
own column so nothing is lost.

Writes 03.outputs/CSV/optical_water_16day.csv
"""

import argparse
import csv
import datetime
import json
import os
import time

import ee

ee.Initialize(project="murphys-deforisk")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "03.outputs", "CSV", "optical_water_16day.csv")
gj = json.load(open(os.path.join(ROOT, "03.outputs", "JSON", "chilwa_basin.geojson")))
AOI = ee.Geometry(gj["features"][0]["geometry"] if gj.get("type") == "FeatureCollection" else gj)

import geopandas as _gpd
_lk = _gpd.read_file(os.path.join(ROOT, "02.inputs", "SHP", "lakes_site.shp")).to_crs(4326)
LAKE = ee.Geometry(json.loads(json.dumps(_lk.geometry.union_all().__geo_interface__)))
WET = LAKE.buffer(15000).intersection(AOI, maxError=100)   # everywhere the lake can reach
WET_KM2 = 2895.0 + 1073.0    # measured: the ten-kilometre belt plus the lake outline

START, END, STEP = datetime.date(1984, 1, 1), datetime.date(2024, 12, 31), 16
AREA = ee.Image.pixelArea().divide(1e6).rename("km2")
THRESH = 0.0

LANDSAT = {"LT04": ("SR_B2", "SR_B5"), "LT05": ("SR_B2", "SR_B5"),
           "LE07": ("SR_B2", "SR_B5"), "LC08": ("SR_B3", "SR_B6"),
           "LC09": ("SR_B3", "SR_B6")}
LS_SPAN = {"LT04": (1982, 1993), "LT05": (1984, 2012), "LE07": (1999, 2025),
           "LC08": (2013, 2025), "LC09": (2021, 2025)}
MODIS = [("MODIS/061/MOD09A1", 2000), ("MODIS/061/MYD09A1", 2002)]


def steps():
    n = (END - START).days // STEP + 1
    return [START + datetime.timedelta(days=STEP * i) for i in range(n)]


def ls_mndwi(t0, t1):
    """Clear Landsat observations in the window, as an MNDWI collection."""
    out = None
    for sensor, (green, swir) in LANDSAT.items():
        lo, hi = LS_SPAN[sensor]
        if t1.year < lo or t0.year > hi:
            continue
        col = (ee.ImageCollection(f"LANDSAT/{sensor}/C02/T1_L2")
               .filterBounds(AOI).filterDate(t0.isoformat(), t1.isoformat()))

        def prep(img, green=green, swir=swir):
            qa = img.select("QA_PIXEL")
            bad = (qa.bitwiseAnd(1 << 1).Or(qa.bitwiseAnd(1 << 2))
                   .Or(qa.bitwiseAnd(1 << 3)).Or(qa.bitwiseAnd(1 << 4))
                   .Or(qa.bitwiseAnd(1 << 5)))
            sr = img.select("SR_B.").multiply(0.0000275).add(-0.2)
            return sr.normalizedDifference([green, swir]).rename("m").updateMask(bad.eq(0))

        col = col.map(prep)
        out = col if out is None else out.merge(col)
    return out


def modis_mndwi(t0, t1):
    """MODIS surface reflectance carries green (b4) and shortwave infrared (b6)."""
    out = None
    for cid, first in MODIS:
        if t1.year < first:
            continue
        col = (ee.ImageCollection(cid).filterBounds(AOI)
               .filterDate(t0.isoformat(), t1.isoformat()))

        def prep(img):
            qa = img.select("StateQA")
            cloud = qa.bitwiseAnd(3).neq(0).Or(qa.bitwiseAnd(1 << 2).neq(0))
            g = img.select("sur_refl_b04").multiply(0.0001)
            s = img.select("sur_refl_b06").multiply(0.0001)
            return (g.subtract(s).divide(g.add(s)).rename("m")
                    .updateMask(cloud.Not()))

        col = col.map(prep)
        out = col if out is None else out.merge(col)
    return out


def stats(col, scale):
    """Water area and clear area in the window, from the median composite."""
    if col is None:
        return None
    # a window with no scene at all gives a median carrying no bands, which the
    # threshold cannot be applied to, so an empty fully-masked band stands in
    # and the row is dropped client-side on its zero scene count
    med = ee.Image(ee.Algorithms.If(col.size().gt(0), col.median(),
                                    ee.Image(0).rename("m").selfMask()))
    # both areas come from one reduction over a two-band image, because Earth
    # Engine caps how many aggregations a single request may carry
    two = (med.gt(THRESH).multiply(AREA).rename("water")
           .addBands(med.mask().multiply(AREA).rename("clear")))
    return ee.Dictionary({
        "n": col.size(),
        "area": two.reduceRegion(ee.Reducer.sum(), AOI, scale,
                                 maxPixels=1e10, bestEffort=True),
        "wet": two.reduceRegion(ee.Reducer.sum(), WET, scale,
                                maxPixels=1e10, bestEffort=True),
    })


def run(family, years, retries=3):
    scale = 90 if family == "landsat" else 500
    maker = ls_mndwi if family == "landsat" else modis_mndwi
    done = {}
    if os.path.exists(OUT):
        for r in csv.DictReader(open(OUT)):
            done[(r["sensor_family"], r["date"])] = r
        print(f"{len(done)} rows already on disk")

    todo = [s for s in steps() if s.year in years and (family, s.isoformat()) not in done]
    CHUNK = 8   # aggregations per request, well inside the concurrency cap
    chunks = [todo[i:i + CHUNK] for i in range(0, len(todo), CHUNK)]
    print(f"{family}: {len(todo)} steps to compute in {len(chunks)} requests")

    for batch in chunks:
        y = f"{batch[0].isoformat()}..{batch[-1].isoformat()}"
        req = {}
        for s in batch:
            d = stats(maker(s, s + datetime.timedelta(days=STEP)), scale)
            if d is not None:
                req[s.isoformat()] = d
        if not req:
            continue
        last = None
        for a in range(retries):
            try:
                got = ee.Dictionary(req).getInfo(); break
            except Exception as e:
                last = e; time.sleep(30 * (a + 1))
        else:
            print(f"  {y}: FAILED {last}", flush=True); continue

        vals = []
        for s in batch:
            g = got.get(s.isoformat())
            if g is None or not g.get("n"):
                continue
            clear = (g.get("area") or {}).get("clear") or 0.0
            water = (g.get("area") or {}).get("water") or 0.0
            wclear = (g.get("wet") or {}).get("clear") or 0.0
            wwater = (g.get("wet") or {}).get("water") or 0.0
            wf = wclear / WET_KM2
            done[(family, s.isoformat())] = {
                "date": s.isoformat(), "year": s.year, "sensor_family": family,
                "n_scenes": g["n"],
                "water_km2": round(wwater, 1), "clear_km2": round(wclear, 1),
                "clear_fraction": round(wf, 4),
                "basin_water_km2": round(water, 1),
                "basin_clear_fraction": round(clear / 8752.0, 4),
                "scale_m": scale,
            }
            if wf >= 0.5:
                vals.append(wwater)
        print(f"  {y}: {sum(1 for s in batch if (family, s.isoformat()) in done)}/{len(batch)} steps"
              + (f", water {min(vals):,.0f} to {max(vals):,.0f} km2 where over half the lake area was clear"
                 if vals else ""), flush=True)

        fields = ["date", "year", "sensor_family", "n_scenes", "water_km2", "clear_km2",
                  "clear_fraction", "basin_water_km2", "basin_clear_fraction", "scale_m"]
        with open(OUT, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields)
            w.writeheader()
            for k in sorted(done, key=lambda k: (k[1], k[0])):
                w.writerow(done[k])
    print(f"wrote {OUT}, {len(done)} rows")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--family", choices=["landsat", "modis"], required=True)
    ap.add_argument("--start", type=int, default=1984)
    ap.add_argument("--end", type=int, default=2024)
    a = ap.parse_args()
    run(a.family, set(range(a.start, a.end + 1)))
