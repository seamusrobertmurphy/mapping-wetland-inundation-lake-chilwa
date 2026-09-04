#!/usr/bin/env python3
"""What optical imagery actually exists over the Lake Chilwa basin, by year.

The committed annual index record is empty for 1985 and 1988, and the note
attached to it says no Landsat 5 scene over the basin clears a 30 per cent
cloud screen in those years. That claim was only ever tested against Landsat 5
Thematic Mapper. Landsat 4 carried the same instrument until 1993 and both
satellites also carried the Multispectral Scanner, so the claim has to be
re-tested against every collection Earth Engine holds before a year is called
empty. The same survey then supplies the scene counts that decide which years
can be composited directly and which need the harmonic model.

Writes 03.outputs/CSV/optical_archive_survey.csv, one row per sensor and year,
with the scene count at three cloud screens.
"""

import csv
import os

import ee

ee.Initialize(project="murphys-deforisk")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import json
gj = json.load(open(os.path.join(ROOT, "03.outputs", "JSON", "chilwa_basin.geojson")))
geom = gj["features"][0]["geometry"] if gj.get("type") == "FeatureCollection" else gj
AOI = ee.Geometry(geom)

COLLECTIONS = {
    "L4_TM":  "LANDSAT/LT04/C02/T1_L2",
    "L5_TM":  "LANDSAT/LT05/C02/T1_L2",
    "L7_ETM": "LANDSAT/LE07/C02/T1_L2",
    "L8_OLI": "LANDSAT/LC08/C02/T1_L2",
    "L9_OLI": "LANDSAT/LC09/C02/T1_L2",
    "L1_MSS": "LANDSAT/LM01/C02/T1",
    "L2_MSS": "LANDSAT/LM02/C02/T1",
    "L3_MSS": "LANDSAT/LM03/C02/T1",
    "L4_MSS": "LANDSAT/LM04/C02/T1",
    "L5_MSS": "LANDSAT/LM05/C02/T1",
    "S2_SR":  "COPERNICUS/S2_SR_HARMONIZED",
}
SCREENS = [30, 60, 101]   # 101 means every scene, screened or not


def main():
    out = os.path.join(ROOT, "03.outputs", "CSV", "optical_archive_survey.csv")
    rows = []
    for name, cid in COLLECTIONS.items():
        try:
            base = ee.ImageCollection(cid).filterBounds(AOI)
            first = base.limit(1).size().getInfo()
        except ee.EEException as e:
            print(f"{name}: collection unavailable, {e}")
            continue
        if first == 0:
            print(f"{name}: no scene over the basin at all")
            continue
        cloud = "CLOUD_COVER" if "T1_L2" in cid else (
            "CLOUDY_PIXEL_PERCENTAGE" if "S2" in cid else None)
        for year in range(1972, 2026):
            col = base.filterDate(f"{year}-01-01", f"{year + 1}-01-01")
            n_all = col.size().getInfo()
            if n_all == 0:
                continue
            row = {"sensor": name, "collection": cid, "year": year}
            for s in SCREENS:
                if cloud is None or s > 100:
                    row[f"n_cloud_lt_{s}"] = n_all
                else:
                    row[f"n_cloud_lt_{s}"] = col.filter(
                        ee.Filter.lt(cloud, s)).size().getInfo()
            rows.append(row)
            print(f"  {name} {year}: {row['n_cloud_lt_101']} scenes, "
                  f"{row['n_cloud_lt_60']} under 60 per cent cloud, "
                  f"{row['n_cloud_lt_30']} under 30", flush=True)

    fields = ["sensor", "collection", "year"] + [f"n_cloud_lt_{s}" for s in SCREENS]
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {out}, {len(rows)} sensor-years")

    print("\nyears 1984 to 2024 with no scene at all under 60 per cent cloud, any sensor:")
    by_year = {}
    for r in rows:
        by_year.setdefault(r["year"], 0)
        by_year[r["year"]] += r["n_cloud_lt_60"]
    empty = [y for y in range(1984, 2025) if by_year.get(y, 0) == 0]
    print(" ", empty if empty else "none")


if __name__ == "__main__":
    main()
