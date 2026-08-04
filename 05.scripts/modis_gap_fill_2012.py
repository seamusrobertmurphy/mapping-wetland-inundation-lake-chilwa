#!/usr/bin/env python3
"""MODIS open-water series over the Lake Chilwa basin, 2000 to 2024.

Purpose: the Landsat record has no 2012, the year that carries both the
recession and the fieldwork. MOD09A1 spans it without interruption at 500 m,
which is also the grain at which imprecisely located field points can validate
a class, because the Typha fringe is one to three kilometres wide.

MOD09A1 rather than the 250 m MOD09Q1: only the 500 m product carries the
green (b04) and shortwave-infrared (b06) pair that MNDWI needs.

Writes 03.outputs/CSV/modis_water_fraction_8day.csv
"""

import csv
import os

import ee

ee.Initialize(project="murphys-deforisk")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASIN = os.path.join(ROOT, "03.outputs", "JSON", "chilwa_basin.geojson")

with open(BASIN) as fh:
    import json

    gj = json.load(fh)
geom = gj["features"][0]["geometry"] if gj.get("type") == "FeatureCollection" else gj
aoi = ee.Geometry(geom)

SCALE = 500


def add_indices(img):
    """MNDWI from green and SWIR1, cloud-screened on StateQA bits 0-1."""
    g = img.select("sur_refl_b04").multiply(0.0001)
    s = img.select("sur_refl_b06").multiply(0.0001)
    mndwi = g.subtract(s).divide(g.add(s)).rename("MNDWI")
    state = img.select("StateQA")
    clear = state.bitwiseAnd(3).eq(0)
    return img.addBands(mndwi).updateMask(clear)


def water_feature(img):
    water = img.select("MNDWI").gt(0).rename("water")
    stats = water.reduceRegion(
        reducer=ee.Reducer.mean().combine(ee.Reducer.count(), sharedInputs=True),
        geometry=aoi,
        scale=SCALE,
        maxPixels=int(1e9),
    )
    return ee.Feature(
        None,
        {
            "date": img.date().format("YYYY-MM-dd"),
            "year": img.date().get("year"),
            "water_fraction": stats.get("water_mean"),
            "n_clear_cells": stats.get("water_count"),
        },
    )


col = (
    ee.ImageCollection("MODIS/061/MOD09A1")
    .filterDate("2000-01-01", "2025-01-01")
    .filterBounds(aoi)
    .map(add_indices)
)

print("images in collection:", col.size().getInfo())

rows = []
for yr in range(2000, 2025):
    sub = col.filterDate(f"{yr}-01-01", f"{yr+1}-01-01")
    fc = ee.FeatureCollection(sub.map(water_feature))
    got = fc.getInfo()["features"]
    for f in got:
        rows.append(f["properties"])
    print(f"  {yr}: {len(got)} steps", flush=True)

out = os.path.join(ROOT, "03.outputs", "CSV", "modis_water_fraction_8day.csv")
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["date", "year", "water_fraction", "n_clear_cells"])
    w.writeheader()
    for r in sorted(rows, key=lambda x: x["date"]):
        w.writerow(r)

print("wrote", out, len(rows), "rows")
