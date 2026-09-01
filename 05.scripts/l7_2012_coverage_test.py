"""Measure what Landsat 7 actually delivers over the Lake Chilwa basin in 2012.

The question this answers: if every valid Landsat 7 observation in a window is
stacked, how much of the basin is left with no observation at all? SLC-off gap
wedges are fixed in scene geometry, so whether stacking multiple passes fills
them is an empirical matter, not something to assume either way.

Counts, per 30 m pixel, the number of cloud-free non-gap observations, then
reports the fraction of the basin at each count.
"""
import json
import ee

ee.Initialize(project="murphys-deforisk")

with open("03.outputs/JSON/chilwa_basin.geojson") as fh:
    gj = json.load(fh)
geom = gj["features"][0]["geometry"] if gj.get("type") == "FeatureCollection" else gj["geometry"]
basin = ee.Geometry(geom)

DILATE = 1  # pixels of buffer applied to the cloud mask


def clear_mask(img):
    """Valid = observed by the sensor and not flagged cloud, shadow, cirrus or snow."""
    qa = img.select("QA_PIXEL")
    bad = (qa.bitwiseAnd(1 << 1)          # dilated cloud
           .Or(qa.bitwiseAnd(1 << 2))     # cirrus
           .Or(qa.bitwiseAnd(1 << 3))     # cloud
           .Or(qa.bitwiseAnd(1 << 4))     # cloud shadow
           .Or(qa.bitwiseAnd(1 << 5)))    # snow
    observed = img.select("SR_B4").mask()  # SLC-off gaps and fill are masked here
    return observed.And(bad.eq(0)).rename("valid").unmask(0)


def coverage(start, end, label):
    col = (ee.ImageCollection("LANDSAT/LE07/C02/T1_L2")
           .filterBounds(basin)
           .filterDate(start, end))
    n = col.size().getInfo()
    if n == 0:
        print(f"{label}: no scenes")
        return
    counts = ee.ImageCollection(col.map(clear_mask)).sum().rename("n").clip(basin)
    area = ee.Image.pixelArea().divide(1e6)
    total = area.updateMask(ee.Image(1).clip(basin)).reduceRegion(
        ee.Reducer.sum(), basin, 30, maxPixels=1e10).get("area").getInfo()
    rows = []
    for k in (0, 1, 2, 3, 5):
        m = counts.eq(0) if k == 0 else counts.gte(k)
        a = area.updateMask(m).reduceRegion(
            ee.Reducer.sum(), basin, 30, maxPixels=1e10).get("area").getInfo()
        rows.append((k, a, 100.0 * a / total))
    print(f"\n{label}: {n} scenes, basin {total:,.0f} km2")
    print(f"  never observed        {rows[0][1]:9,.0f} km2  {rows[0][2]:5.1f} %")
    for k, a, pct in rows[1:]:
        print(f"  >= {k} clear looks      {a:9,.0f} km2  {pct:5.1f} %")


if __name__ == "__main__":
    coverage("2012-02-01", "2012-07-01", "Fieldwork window, Feb to Jun 2012")
    coverage("2012-01-01", "2013-01-01", "Whole of 2012")


def by_path(start, end, label):
    """Does a single path fill its own gaps, or is it the 166/167 sidelap that does it?"""
    base = (ee.ImageCollection("LANDSAT/LE07/C02/T1_L2")
            .filterBounds(basin).filterDate(start, end))
    area = ee.Image.pixelArea().divide(1e6)
    total = area.updateMask(ee.Image(1).clip(basin)).reduceRegion(
        ee.Reducer.sum(), basin, 30, maxPixels=1e10).get("area").getInfo()
    print(f"\n{label}, by WRS-2 path")
    for path in (166, 167):
        col = base.filter(ee.Filter.eq("WRS_PATH", path))
        n = col.size().getInfo()
        if n == 0:
            print(f"  path {path}: no scenes")
            continue
        counts = ee.ImageCollection(col.map(clear_mask)).sum().clip(basin)
        seen = area.updateMask(counts.gte(1)).reduceRegion(
            ee.Reducer.sum(), basin, 30, maxPixels=1e10).get("area").getInfo()
        print(f"  path {path}: {n:2d} scenes, {seen:7,.0f} km2 seen at least once"
              f" ({100*seen/total:5.1f} % of basin)")


def near_date(centre, days, label):
    """Coverage in a tight window, which is what a single-date map actually needs."""
    c = ee.Date(centre)
    coverage(c.advance(-days, "day").format("YYYY-MM-dd").getInfo(),
             c.advance(days, "day").format("YYYY-MM-dd").getInfo(),
             f"{label} ({centre} +/- {days} d)")
