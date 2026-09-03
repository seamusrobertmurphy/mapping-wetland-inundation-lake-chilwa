#!/usr/bin/env python3
"""Rebuild the annual spectral-index series with Landsat 7 included.

03.outputs/CSV/landsat_annual_indices.csv holds 37 years and is missing 1985,
1988, 2002 and 2012, because it was built under the Landsat 5 and Landsat 8
rule. Landsat 7 flew from 1999, so it can close 2002 and 2012. The two 1980s
gaps can only be closed by Landsat 5 and are reported here rather than assumed.

The method follows the committed series: cloud, shadow, cirrus and snow masked
from QA_PIXEL, scenes over 30 per cent cloud dropped, an annual median composite
per index, then the spatial mean over the basin. Overlapping years are
recomputed as a check, so a method that has drifted shows up as a mismatch
against the committed file rather than passing silently.

Writes 03.outputs/CSV/landsat_annual_indices_rebuilt.csv
"""

import csv
import json
import os

import ee

ee.Initialize(project="murphys-deforisk")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "03.outputs", "CSV")

gj = json.load(open(os.path.join(ROOT, "03.outputs", "JSON", "chilwa_basin.geojson")))
geom = gj["features"][0]["geometry"] if gj.get("type") == "FeatureCollection" else gj["geometry"]
AOI = ee.Geometry(geom)

MAX_CLOUD = 30
SENSORS = {                       # blue, green, red, nir, swir1, swir2
    "LT05": ("SR_B1", "SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B7"),
    "LE07": ("SR_B1", "SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B7"),
    "LC08": ("SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B6", "SR_B7"),
}
NAMES = ["blue", "green", "red", "nir", "swir1", "swir2"]


def harmonise(img, bands):
    """Scale to reflectance, mask cloud and shadow, rename to a common schema."""
    qa = img.select("QA_PIXEL")
    bad = (qa.bitwiseAnd(1 << 1).Or(qa.bitwiseAnd(1 << 2)).Or(qa.bitwiseAnd(1 << 3))
             .Or(qa.bitwiseAnd(1 << 4)).Or(qa.bitwiseAnd(1 << 5)))
    sr = img.select(list(bands)).multiply(0.0000275).add(-0.2).rename(NAMES)
    return sr.updateMask(bad.eq(0)).copyProperties(img, ["system:time_start"])


def indices(img):
    """The five indices the committed series carries."""
    b = {n: img.select(n) for n in NAMES}
    ndwi = b["green"].subtract(b["nir"]).divide(b["green"].add(b["nir"]))
    mndwi = b["green"].subtract(b["swir1"]).divide(b["green"].add(b["swir1"]))
    aweish = (b["blue"].add(b["green"].multiply(2.5))
              .subtract(b["nir"].add(b["swir1"]).multiply(1.5))
              .subtract(b["swir2"].multiply(0.25)))
    wri = b["green"].add(b["red"]).divide(b["nir"].add(b["swir1"]))
    ndpi = b["swir1"].subtract(b["green"]).divide(b["swir1"].add(b["green"]))
    return (ndwi.rename("NDWI").addBands(mndwi.rename("MNDWI"))
            .addBands(aweish.rename("AWEIsh")).addBands(wri.rename("WRI"))
            .addBands(ndpi.rename("NDPI")))


def year_row(year, sensors):
    """One year of the series, or None where no scene survives the cloud screen."""
    col = None
    counts = {}
    for s in sensors:
        c = (ee.ImageCollection(f"LANDSAT/{s}/C02/T1_L2")
             .filterBounds(AOI).filterDate(f"{year}-01-01", f"{year+1}-01-01")
             .filter(ee.Filter.lt("CLOUD_COVER", MAX_CLOUD)))
        n = c.size().getInfo()
        counts[s] = n
        if n == 0:
            continue
        c = c.map(lambda i, b=SENSORS[s]: harmonise(i, b)).map(indices)
        col = c if col is None else col.merge(c)
    if col is None:
        return None, counts
    stats = (ee.ImageCollection(col).median()
             .reduceRegion(ee.Reducer.mean(), AOI, 30, maxPixels=1e10).getInfo())
    return stats, counts


def main():
    committed = {int(float(r["year"])): r for r in
                 csv.DictReader(open(os.path.join(OUT, "landsat_annual_indices.csv")))}
    rows, checks = [], []
    for year in range(1984, 2025):
        sensors = [s for s in ("LT05", "LE07", "LC08")
                   if not (s == "LE07" and year < 1999)
                   and not (s == "LC08" and year < 2013)
                   and not (s == "LT05" and year > 2011)]
        stats, counts = year_row(year, sensors)
        n = sum(counts.values())
        if stats is None or stats.get("MNDWI") is None:
            print(f"  {year}  no usable scene under {MAX_CLOUD} % cloud  {counts}")
            continue
        row = {"year": year, "n_scenes": n,
               **{k: stats.get(k) for k in ("NDWI", "MNDWI", "AWEIsh", "WRI", "NDPI")},
               "sensors": "+".join(s for s, c in counts.items() if c)}
        rows.append(row)
        flag = ""
        if year in committed:
            old = float(committed[year]["MNDWI"])
            d = abs(old - row["MNDWI"])
            checks.append(d)
            flag = f"  committed MNDWI {old:+.4f}, diff {d:.4f}"
        else:
            flag = "  NEW, absent from the committed series"
        print(f"  {year}  n={n:3d}  MNDWI {row['MNDWI']:+.4f}  [{row['sensors']}]{flag}")

    with open(os.path.join(OUT, "landsat_annual_indices_rebuilt.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\n{len(rows)} years written, was {len(committed)}")
    if checks:
        print(f"agreement on shared years: median |diff| in MNDWI "
              f"{sorted(checks)[len(checks)//2]:.4f}, worst {max(checks):.4f}")


if __name__ == "__main__":
    main()
