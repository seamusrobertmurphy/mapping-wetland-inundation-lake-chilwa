#!/usr/bin/env python3
"""Water standing within vegetation at Lake Chilwa from ALOS PALSAR, by year.

The same measurement the Envisat work makes for 2011 and 2012, extended to the
years an L-band record exists. PALSAR transmits horizontally, so a wave passing
between vertical reed stems reflects off the water, bounces into the stem base
and returns doubled, brightening HH while HV, which comes from volume
scattering inside the stand, does not share the gain. The difference HH minus
HV is therefore high over flooded vegetation and the incidence-angle term
cancels in it, because both channels share one geometry.

Sentinel-1 cannot do this over this basin. It transmits vertically, which the
stems absorb, and a check of the archive on 3 September 2026 found exactly one
horizontally polarised scene, an Extra Wide swath frame in 2023 covering 537 of
the basin's 8,752 square kilometres, so the horizontal record is PALSAR and
Envisat alone.

The decision rule is the one settled on the Envisat frames. Both thresholds sit
at the ninety-fifth percentile of dry land more than ten kilometres from the
lake outline, so the rule admits at most five per cent of dry land on either
channel, and the rate it actually admits there is reported as the noise floor.
Open water comes from the optical harmonic model read on the mosaic's own
median acquisition date, never from the radar, and the flooded class is taken
only outside it. That date matters and is not assumed: the 2016 mosaic was
acquired on 27 June, in the dry season, so a radar figure compared against an
annual optical mean would be comparing two different times of year. The mosaic
carries its acquisition time in its own epoch band, that time is read back
before anything else is computed, and both the open-water layer and the optical
estimate of the vegetated class are evaluated on it, so the radar and the
optical measurement being set against it belong to the same week.

Calibration follows Shimada et al. (2009): gamma0 in decibels is
10*log10(DN^2) - 83.0 for the yearly mosaic.

Writes 03.outputs/CSV/palsar_flooded_vegetation.csv
"""

import csv
import datetime
import json
import math
import os
import time

import ee

import harmonic_open_water_series as H
import sma_harmonic_flooded_veg as S

COMPOSITE_DAYS = 45   # Landsat looks either side of the mosaic date

ee.Initialize(project="murphys-deforisk")

ROOT = H.ROOT
OUT = os.path.join(ROOT, "03.outputs", "CSV", "palsar_flooded_vegetation.csv")
AOI = H.AOI
AREA = H.AREA_KM2

COLLECTION = "JAXA/ALOS/PALSAR/YEARLY/SAR_EPOCH"
CAL_OFFSET = -83.0        # Shimada et al. (2009), yearly mosaic
WATER_MNDWI = 0.10        # above this the optical model calls it open water
BELT_KM = 10.0            # the reed fringe lies within this of the lake outline
FAR_KM = 10.0             # beyond this is the dry-land reference
DRY_PCTL = 95


def db(img, band):
    return img.select(band).pow(2).log10().multiply(10).add(CAL_OFFSET).rename(band)


def _unused_window_for(year):
    """The window and model order the open-water series already chose for this year.

    Read from 03.outputs/CSV/harmonic_open_water_annual.csv rather than probed
    again, so the optical layer the radar is judged against is the same surface
    that the open-water column reports, and so the expensive probe runs once."""
    path = os.path.join(ROOT, "03.outputs", "CSV", "harmonic_open_water_annual.csv")
    if not os.path.exists(path):
        raise SystemExit("run harmonic_open_water_series.py first; its windows are needed here")
    for r in csv.DictReader(open(path)):
        if int(r["year"]) == year:
            harm = H.FULL_HARMONICS if r["model_order"] == "annual+semiannual" else H.SIMPLE_HARMONICS
            return int(r["half_window"]), harm
    raise SystemExit(f"{year} is not in the open-water series yet")


def year_row(year, lake, belt, far, retries=3):
    col = (ee.ImageCollection(COLLECTION).filterBounds(AOI)
           .filterDate(f"{year}-01-01", f"{year + 1}-01-01"))
    if col.size().getInfo() == 0:
        return None
    img = ee.Image(col.mosaic())

    # read the mosaic's own acquisition time first, so every optical layer it
    # is compared against is evaluated on the day the radar looked
    epoch_ms = ee.Number(img.select("epoch").reduceRegion(
        ee.Reducer.median(), AOI, 300, maxPixels=1e10, bestEffort=True).get("epoch")).getInfo()
    mdate = datetime.datetime.fromtimestamp(epoch_ms / 1000.0, datetime.UTC).date()
    mid = ee.Date(mdate.isoformat())

    hh = db(img, "HH")
    hv = db(img, "HV")
    diff = hh.subtract(hv).rename("d")
    valid = hh.mask().And(hv.mask())

    # Open water and the optical estimate of the vegetated class both come from
    # the Landsat scenes nearest the mosaic date, not from a fitted annual
    # surface. That keeps the comparison observational on both sides and
    # removes any dependence on a model built for a different product.
    t0 = mid.advance(-COMPOSITE_DAYS, "day")
    t1 = mid.advance(COMPOSITE_DAYS, "day")
    sens = [k for k, (lo, hi) in H.LS_SPAN.items()] if hasattr(H, "LS_SPAN") else \
        H.sensors_for(year, year)

    ls = None
    for sensor in H.sensors_for(year - 1, year + 1):
        green, swir = H.BANDS[sensor]
        c = (ee.ImageCollection(f"LANDSAT/{sensor}/C02/T1_L2")
             .filterBounds(AOI).filterDate(t0, t1).map(H.clear))

        def prep(img, green=green, swir=swir, sensor=sensor):
            sr = img.select("SR_B.").multiply(0.0000275).add(-0.2)
            m = sr.normalizedDifference([green, swir]).rename("m")
            spec = [S.EM[S.ERA[sensor]][cl] for cl in S.CLASSES]
            ref = img.select(S.BANDS6[sensor]).multiply(0.0000275).add(-0.2).rename(S.SIX)
            v = ref.unmix(spec, True, True).rename(S.CLASSES).select("flooded_veg").rename("v")
            return m.addBands(v)

        c = c.map(prep)
        ls = c if ls is None else ls.merge(c)

    n_ls = ls.size()
    comp = ee.Image(ee.Algorithms.If(n_ls.gt(0), ls.median(),
                                     ee.Image.cat(ee.Image(0).rename("m"),
                                                  ee.Image(0).rename("v")).selfMask()))
    mndwi = comp.select("m")
    veg_frac = comp.select("v").clamp(0, 1).rename("v")

    water = mndwi.gt(WATER_MNDWI)
    nonwater = mndwi.lte(WATER_MNDWI).And(valid)
    dryland = nonwater.And(far)

    def pct(image, mask, p):
        return ee.Number(image.updateMask(mask).reduceRegion(
            ee.Reducer.percentile([p]), AOI, 100, maxPixels=1e10,
            bestEffort=True).values().get(0))

    t_diff = pct(diff, dryland, DRY_PCTL)
    t_hh = pct(hh, dryland, DRY_PCTL)
    bright = diff.gt(t_diff).And(hh.gt(t_hh)).And(nonwater)

    def km2(mask):
        return ee.Number(mask.multiply(AREA).rename("km2").reduceRegion(
            ee.Reducer.sum(), AOI, 100, maxPixels=1e10, bestEffort=True).get("km2"))

    stats = ee.Dictionary({
        "n_landsat": n_ls, "t_diff": t_diff, "t_hh": t_hh,
        "flooded_belt": km2(bright.And(belt)),
        "flooded_basin": km2(bright),
        "nonwater_belt": km2(nonwater.And(belt)),
        "dry_bright": km2(bright.And(far)),
        "dry_total": km2(dryland),
        "water_km2": km2(water.And(valid)),
        "valid_km2": km2(valid),
        # the optical estimate of the vegetated class on the same day, so the
        # two methods are compared like for like rather than across seasons
        "optical_veg_same_date": ee.Number(veg_frac.multiply(AREA).rename("km2")
                                           .reduceRegion(ee.Reducer.sum(), AOI, H.SCALE,
                                                         maxPixels=1e10,
                                                         bestEffort=True).get("km2")),
        "hh_water": pct(hh, water.And(valid), 50),
        "hh_dry": pct(hh, dryland, 50),
        "hv_dry": pct(hv, dryland, 50),
    })
    last = None
    for a in range(retries):
        try:
            g = stats.getInfo(); break
        except Exception as e:
            last = e; time.sleep(20 * (a + 1))
    else:
        raise RuntimeError(f"{year}: {last}")

    fpr = (g["dry_bright"] / g["dry_total"]) if g["dry_total"] else float("nan")
    share = (g["flooded_belt"] / g["nonwater_belt"]) if g["nonwater_belt"] else float("nan")
    return {
        "year": year, "sensor": "ALOS PALSAR yearly mosaic",
        "mosaic_date": mdate.isoformat(),
        "n_landsat_in_composite": g.get("n_landsat"),
        "valid_km2": round(g["valid_km2"], 1),
        "optical_open_water_km2": round(g["water_km2"], 1),
        "optical_flooded_same_date_km2": round(g["optical_veg_same_date"], 1),
        "flooded_veg_belt_km2": round(g["flooded_belt"], 1),
        "flooded_veg_basin_km2": round(g["flooded_basin"], 1),
        "nonwater_belt_km2": round(g["nonwater_belt"], 1),
        "flooded_share_of_belt": round(share, 4),
        "dry_land_false_positive": round(fpr, 4),
        "t_hh_minus_hv_dB": round(g["t_diff"], 2),
        "t_hh_dB": round(g["t_hh"], 2),
        "hh_water_median_dB": round(g["hh_water"], 2),
        "hh_dry_median_dB": round(g["hh_dry"], 2),
        "hv_dry_median_dB": round(g["hv_dry"], 2),
        "method": f"HH minus HV above the dry-land {DRY_PCTL}th percentile and HH above the "
                  f"same, outside optical water at MNDWI > {WATER_MNDWI}, within "
                  f"{BELT_KM:.0f} km of the lake outline, all layers read on {mdate.isoformat()}",
    }


def main():
    import geopandas as gpd
    g = gpd.read_file(os.path.join(ROOT, "02.inputs", "SHP", "lakes_site.shp")).to_crs(4326)
    lake = ee.Geometry(json.loads(json.dumps(g.geometry.union_all().__geo_interface__)))
    belt = ee.Image.constant(1).clip(lake.buffer(BELT_KM * 1000)).mask().gt(0)
    far = ee.Image.constant(1).clip(lake.buffer(FAR_KM * 1000)).mask().gt(0).Not()

    done = {}
    if os.path.exists(OUT):
        for r in csv.DictReader(open(OUT)):
            done[int(r["year"])] = r
    rows = []
    for year in list(range(2007, 2011)) + list(range(2015, 2025)):
        if year in done:
            rows.append(done[year]); continue
        r = year_row(year, lake, belt, far)
        if r is None:
            print(f"  {year}: no mosaic"); continue
        rows.append(r); done[year] = r
        print(f"  {year} ({r['mosaic_date']}): flooded {r['flooded_veg_belt_km2']:7.0f} km2 in the belt "
              f"({100 * r['flooded_share_of_belt']:.1f}% of non-water there), open water "
              f"{r['optical_open_water_km2']:7.0f} km2, false positive on dry land "
              f"{100 * r['dry_land_false_positive']:.1f}%, thresholds "
              f"{r['t_hh_minus_hv_dB']:.1f} and {r['t_hh_dB']:.1f} dB; optical says "
              f"{r['optical_flooded_same_date_km2']:.0f} km2 on the same day", flush=True)
        with open(OUT, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            for y in sorted(done):
                w.writerow(done[y])
    print(f"wrote {OUT}, {len(rows)} years")


if __name__ == "__main__":
    main()
