#!/usr/bin/env python3
"""AVHRR over Lake Chilwa on the sixteen-day grid, 1984 to 2022.

Why AVHRR at all. No other sensor sees this basin frequently before 2000.
Landsat gives four to thirty clear scenes a year in that era and MODIS does not
launch until February 2000, so the sixteen-day record for 1984 to 1999 has no
dense stream unless AVHRR supplies it. The catalogue holds 14,903 AVHRR images
over the basin from June 1981, daily at 5.6 km.

What AVHRR can and cannot do here. Five-kilometre pixels put roughly thirty-five
whole pixels across a thousand square kilometres of lake, so a shoreline cannot
be mapped and a water area cannot be measured by counting pixels; the
shoreline pixels are all mixtures. What it can do is carry a basin-scale index
that moves with the water, and be calibrated against a sensor that does measure
area. AVHRR has no shortwave infrared on the early platforms, so a modified
water index is impossible and the normalised difference vegetation index is
used instead: water is dark in the near infrared and so reads low, vegetation
and dry ground read high. Checked over three years the separation between the
lake and the land beyond ten kilometres was 0.292 in the wet year 1990, 0.175
in 2012 and 0.040 in the 1995 recession, when much of the lake bed was exposed
and genuinely was land, so the index tracks the hydrology rather than failing.

The index is deliberately a contrast, the lake against its own surroundings on
the same day, not an absolute reflectance. That cancels the sensor drift,
orbital decay and changing overpass time that make raw AVHRR reflectance
unusable across a forty-year record spanning seven satellites.

A median over the lake polygon was tried first and does not work. Tested against
451 well-observed optical steps it reached a correlation of only 0.379 and
compressed a true range of 7 to 1,570 square kilometres into 739 to 1,283,
because a median saturates: a lake half emptied still shows a water median while
its area has halved. Area needs a fraction, not a central value. Every pixel is
therefore unmixed against two endmembers, open water at an index of -0.05 and
the same-date land reference, and the resulting fractions are integrated over
the lake and a fifteen-kilometre belt around it. That is linear in area by
construction and it keeps the same-date referencing that cancels the drift.

Calibration to area happens elsewhere, against MODIS over their 2000 to 2022
overlap, which is twenty-two years of coincident observation.

Writes 03.outputs/CSV/avhrr_index_16day.csv
"""

import argparse
import csv
import datetime
import json
import os
import time

import ee
import geopandas as gpd

ee.Initialize(project="murphys-deforisk")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "03.outputs", "CSV", "avhrr_index_16day.csv")

gj = json.load(open(os.path.join(ROOT, "03.outputs", "JSON", "chilwa_basin.geojson")))
AOI = ee.Geometry(gj["features"][0]["geometry"] if gj.get("type") == "FeatureCollection" else gj)
_g = gpd.read_file(os.path.join(ROOT, "02.inputs", "SHP", "lakes_site.shp")).to_crs(4326)
LAKE = ee.Geometry(json.loads(json.dumps(_g.geometry.union_all().__geo_interface__)))
REF = AOI.difference(LAKE.buffer(10000), maxError=100)   # land beyond 10 km of the lake
WET = LAKE.buffer(15000)          # everywhere the lake can reach
NDVI_WATER = -0.05                # open water endmember, from the observed lake minima

START = datetime.date(1984, 1, 1)
END = datetime.date(2024, 12, 31)
STEP = 16
SCALE = 5566          # the collection's own grid


def steps():
    n = (END - START).days // STEP + 1
    return [START + datetime.timedelta(days=STEP * i) for i in range(n)]


def window_stats(t0, t1):
    """Median NDVI in the lake and on the far reference, plus how much was clear."""
    col = ee.ImageCollection("NOAA/CDR/AVHRR/SR/V5").filterBounds(AOI).filterDate(
        t0.isoformat(), t1.isoformat())

    def prep(im):
        qa = im.select("QA")
        # NOAA CDR QA: bit 1 cloudy, bit 2 cloud shadow, bit 6 night
        bad = qa.bitwiseAnd(1 << 1).Or(qa.bitwiseAnd(1 << 2)).Or(qa.bitwiseAnd(1 << 6))
        nd = im.normalizedDifference(["SREFL_CH2", "SREFL_CH1"]).rename("ndvi")
        return nd.updateMask(bad.eq(0))

    clear = col.map(prep)
    med = clear.median()
    cnt = clear.count().rename("n")

    # two-endmember unmixing against the same-day land reference, integrated
    # sub-pixel over everywhere the lake can reach
    ref_ndvi = ee.Number(med.reduceRegion(ee.Reducer.median(), REF, SCALE,
                                          maxPixels=1e9, bestEffort=True).get("ndvi"))
    span = ref_ndvi.subtract(NDVI_WATER).max(0.05)
    frac = (ee.Image.constant(ref_ndvi).subtract(med)
            .divide(ee.Image.constant(span)).clamp(0, 1).rename("f"))
    area = ee.Number(frac.multiply(ee.Image.pixelArea().divide(1e6)).rename("km2")
                     .reduceRegion(ee.Reducer.sum(), WET, SCALE,
                                   maxPixels=1e9, bestEffort=True).get("km2"))
    seen = ee.Number(med.mask().multiply(ee.Image.pixelArea().divide(1e6)).rename("km2")
                     .reduceRegion(ee.Reducer.sum(), WET, SCALE,
                                   maxPixels=1e9, bestEffort=True).get("km2"))
    return ee.Dictionary({
        "water_km2": area, "seen_km2": seen,
        "n_images": col.size(),
        "lake": med.reduceRegion(ee.Reducer.median(), LAKE, SCALE,
                                 maxPixels=1e9, bestEffort=True).get("ndvi"),
        "ref": med.reduceRegion(ee.Reducer.median(), REF, SCALE,
                                maxPixels=1e9, bestEffort=True).get("ndvi"),
        "lake_obs": cnt.reduceRegion(ee.Reducer.mean(), LAKE, SCALE,
                                     maxPixels=1e9, bestEffort=True).get("n"),
    })


def run(years, retries=3):
    done = {}
    if os.path.exists(OUT):
        for r in csv.DictReader(open(OUT)):
            done[r["date"]] = r
        print(f"{len(done)} steps already on disk")

    todo = [s for s in steps() if s.year in years and s.isoformat() not in done]
    print(f"{len(todo)} steps to compute")
    # one request per year keeps the aggregation count well inside the limit
    by_year = {}
    for s in todo:
        by_year.setdefault(s.year, []).append(s)

    for y in sorted(by_year):
        batch = by_year[y]
        req = ee.Dictionary({s.isoformat(): window_stats(s, s + datetime.timedelta(days=STEP))
                             for s in batch})
        last = None
        for a in range(retries):
            try:
                got = req.getInfo(); break
            except Exception as e:
                last = e; time.sleep(30 * (a + 1))
        else:
            print(f"  {y}: FAILED {last}", flush=True); continue

        for s in batch:
            g = got[s.isoformat()]
            lake, ref = g.get("lake"), g.get("ref")
            done[s.isoformat()] = {
                "date": s.isoformat(), "year": s.year,
                "n_images": g.get("n_images") or 0,
                "water_km2": round(g["water_km2"], 1) if g.get("water_km2") is not None else "",
                "seen_km2": round(g["seen_km2"], 1) if g.get("seen_km2") is not None else "",
                "mean_clear_obs_in_lake": round(g.get("lake_obs") or 0, 2),
                "ndvi_lake": round(lake, 5) if lake is not None else "",
                "ndvi_reference": round(ref, 5) if ref is not None else "",
                "contrast": round(ref - lake, 5) if (lake is not None and ref is not None) else "",
            }
        ok = sum(1 for s in batch if done[s.isoformat()]["contrast"] != "")
        vals = [done[s.isoformat()]["water_km2"] for s in batch
                if done[s.isoformat()].get("water_km2") not in ("", None)]
        print(f"  {y}: {ok}/{len(batch)} steps, water {min(vals):,.0f} to {max(vals):,.0f} km2"
              if vals else f"  {y}: no values", flush=True)

        fields = ["date", "year", "n_images", "water_km2", "seen_km2",
                  "mean_clear_obs_in_lake", "ndvi_lake", "ndvi_reference", "contrast"]
        with open(OUT, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields)
            w.writeheader()
            for k in sorted(done):
                w.writerow(done[k])
    print(f"wrote {OUT}, {len(done)} steps")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", type=int, default=1984)
    ap.add_argument("--end", type=int, default=2022)
    a = ap.parse_args()
    run(set(range(a.start, a.end + 1)))
