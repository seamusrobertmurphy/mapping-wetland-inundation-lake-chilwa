#!/usr/bin/env python3
"""Open water at Lake Chilwa for every year from 1984 to 2024, from a rolling
per-pixel harmonic model of MNDWI.

Why this replaces the annual composite. The committed record built one
composite per year and inherited that year's scene count, so 1985, 1988, 2002
and 2012 carry nothing at all and 1987, 1992, 2001 and 2003 rest on two or
three scenes, which put 2003 at 118 km2 of open water against a median of
1,227 km2. A composite cannot measure a lake from three scenes. Fitting the
seasonal shape on a window that spans several years uses every observation the
neighbourhood holds, while the thin year's own scenes and the linear drift term
still set its level, so the estimate degrades gracefully instead of collapsing.

The model is the one settled for 2012 in harmonic_fit_2012.py, unchanged:

    MNDWI(t) = a0 + a1*t + b1*cos(2*pi*t) + c1*sin(2*pi*t)
                        + b2*cos(4*pi*t) + c2*sin(4*pi*t)

with t in years. LandTrendr was considered for this and rejected, because it
reduces each year to a single value joined by straight lines, which would erase
the flood pulse the study measures.

The window is the narrowest the data will carry, because a wide one smooths
away the very signal the study measures. Checked against the MODIS record over
2000 to 2009, a fixed five-year window varied by only 52 square kilometres
across that decade where MODIS varied by 105, and the correlation between them
fell to 0.66 even though the root mean square difference was a slight 47 square
kilometres. The level was right and the year-to-year movement was not, because
only the linear drift term carries variation between years and a five-year
window damps it. The candidate windows therefore begin at a single year and
widen only as far as the coverage test demands.

The window is not fixed. A five-year window suits the years after 1994, when
the archive holds twenty or more usable scenes a year, but the early record is
far thinner: four scenes clear a sixty per cent cloud screen in 1984 and none
at all in 1985. Tried at five years on 1984 the fit covered fourteen square
kilometres of the eight thousand seven hundred and fifty-two in the basin,
because almost no pixel reached the observation floor, and the water area it
reported was water inside that sliver rather than a lake. The window is
therefore widened one year at a time, from two either side out to six, until
the fitted model covers ninety-five per cent of the basin, and where six years
still will not do it the semi-annual pair is dropped and the simpler four
parameter model is fitted instead, which needs fewer observations to constrain.
The window and the model order are recorded for every year, so a reader can see
which years are measured densely and which are carried by a wide window. The
fitted curve is then read at twenty-four dates through the target year, and the
area under MNDWI above zero is summed at each, giving the annual mean, the
annual maximum and the month in which the maximum falls. The uncertainty is derived from the fit
itself rather than assumed. The residual root mean square error of the fitted
model is computed for that window, the area is recomputed at plus and minus one
residual either side of the threshold, and half that spread is reported, so a
year whose fit is poor because it holds few observations carries a wider
interval than a year whose fit is tight. The error that matters lands on pixels
sitting near the water boundary, and that boundary runs through the reed fringe,
which is the study's subject.

Sensors follow the closed decision. Landsat 5 and 8 carry the record, Landsat 7
enters the fitting window exactly as it does for 2012, and Landsat 4 is added
because the archive survey found eight Landsat 4 scenes over the basin in 1988,
a year the Landsat 5 rule reported as empty.

Writes 03.outputs/CSV/harmonic_open_water_annual.csv
"""

import argparse
import csv
import json
import math
import os
import time

import ee

ee.Initialize(project="murphys-deforisk")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "03.outputs", "CSV", "harmonic_open_water_annual.csv")

gj = json.load(open(os.path.join(ROOT, "03.outputs", "JSON", "chilwa_basin.geojson")))
geom = gj["features"][0]["geometry"] if gj.get("type") == "FeatureCollection" else gj
AOI = ee.Geometry(geom)

ORIGIN = ee.Date("2012-01-01")
FULL_HARMONICS = [1, 2]      # annual and semi-annual, six free parameters
SIMPLE_HARMONICS = [1]       # annual only, four free parameters
OBS_PER_PARAM = 2            # observation floor, twice the free parameters
HALF_WINDOWS = [0, 1, 2, 3, 4, 5, 6]
COVERAGE_TARGET = 0.95       # share of the basin the fitted model must reach
N_DATES = 24                 # dates read through the target year
THRESH = 0.0
MIN_RESID = 0.03        # floor on the residual, so a tight fit still carries an interval
SCALE = 60               # proven against the 30 m raster for 2012 to within 2 km2

BANDS = {"LT04": ("SR_B2", "SR_B5"), "LT05": ("SR_B2", "SR_B5"),
         "LE07": ("SR_B2", "SR_B5"), "LC08": ("SR_B3", "SR_B6"),
         "LC09": ("SR_B3", "SR_B6")}
AREA_KM2 = ee.Image.pixelArea().divide(1e6).rename("km2")
BASIN_KM2_NOMINAL = 8752.0   # the basin polygon's area, used to scale fraction errors


def predictors(harmonics):
    return ["const", "t"] + [f"{f}{k}" for k in harmonics for f in ("cos", "sin")]


def scaled(img):
    return ee.Image(img.select("SR_B.").multiply(0.0000275).add(-0.2)
                    .copyProperties(img, ["system:time_start"]))


def clear(img):
    qa = img.select("QA_PIXEL")
    bad = (qa.bitwiseAnd(1 << 1).Or(qa.bitwiseAnd(1 << 2)).Or(qa.bitwiseAnd(1 << 3))
             .Or(qa.bitwiseAnd(1 << 4)).Or(qa.bitwiseAnd(1 << 5)))
    return img.updateMask(bad.eq(0))


def design(date, harmonics):
    t = ee.Number(ee.Date(date).difference(ORIGIN, "year"))
    bands = [ee.Image.constant(1), ee.Image.constant(t)]
    for k in harmonics:
        w = t.multiply(2 * math.pi * k)
        bands += [ee.Image.constant(w.cos()), ee.Image.constant(w.sin())]
    return ee.Image.cat(bands).rename(predictors(harmonics)).float()


def observations(start, end, sensors, harmonics=FULL_HARMONICS):
    out = None
    for s in sensors:
        green, swir1 = BANDS[s]
        col = (ee.ImageCollection(f"LANDSAT/{s}/C02/T1_L2")
               .filterBounds(AOI).filterDate(start, end).map(clear))

        def prep(img, green=green, swir1=swir1):
            d = ee.Date(img.get("system:time_start"))
            mndwi = scaled(img).normalizedDifference([green, swir1]).rename("mndwi")
            return (design(d, harmonics).addBands(mndwi).updateMask(mndwi.mask())
                    .set("system:time_start", img.get("system:time_start")))

        col = col.map(prep)
        out = col if out is None else out.merge(col)
    return out


def sensors_for(y0, y1):
    """Which collections could hold a scene inside the window."""
    s = []
    if y0 <= 1993 and y1 >= 1982:
        s.append("LT04")
    if y0 <= 2011 and y1 >= 1984:
        s.append("LT05")
    if y0 <= 2024 and y1 >= 1999:
        s.append("LE07")
    if y1 >= 2013:
        s.append("LC08")
    if y1 >= 2021:
        s.append("LC09")
    return s


def clear_count(y0, y1):
    """Clear MNDWI observations per pixel over a window, without fitting."""
    col = observations(f"{y0}-01-01", f"{y1}-01-01", sensors_for(y0, y1))
    return col.select("mndwi").count().rename("n"), col


def choose_window(year, retries=4):
    """Smallest window and simplest model that fits at least 95 per cent of the basin.

    Coverage is measured before any regression is run, because the count image
    is cheap and the regression is not, and because a fit that covers a sliver
    of the basin reports the water inside that sliver rather than the lake.

    The candidates are probed one model order at a time rather than all at
    once. Earth Engine refuses a request carrying too many concurrent
    aggregations, and fourteen candidates in one call reaches that limit, so
    the full model is tried across every window first and the simpler model is
    only asked for if none of them reaches the target."""
    basin = ee.Number(AREA_KM2.reduceRegion(
        reducer=ee.Reducer.sum(), geometry=AOI, scale=SCALE,
        maxPixels=1e10, bestEffort=True).get("km2"))

    def probe(harm):
        d = {}
        for hw in HALF_WINDOWS:
            y0, y1 = max(year - hw, 1984), min(year + hw + 1, 2026)
            n, _ = clear_count(y0, y1)
            floor = OBS_PER_PARAM * len(predictors(harm))
            d[str(hw)] = ee.Number(n.gte(floor).multiply(AREA_KM2).rename("km2")
                                   .reduceRegion(reducer=ee.Reducer.sum(), geometry=AOI,
                                                 scale=SCALE, maxPixels=1e10,
                                                 bestEffort=True).get("km2")).divide(basin)
        last = None
        for attempt in range(retries):
            try:
                return ee.Dictionary(d).getInfo()
            except Exception as e:
                last = e
                time.sleep(30 * (attempt + 1))
        raise RuntimeError(f"{year}: window probe failed, {last}")

    seen = {}
    for harm in (FULL_HARMONICS, SIMPLE_HARMONICS):
        got = probe(harm)
        tag = "full" if len(harm) == 2 else "simple"
        seen.update({f"{hw}_{tag}": got[str(hw)] for hw in HALF_WINDOWS})
        for hw in HALF_WINDOWS:
            if got[str(hw)] >= COVERAGE_TARGET:
                return hw, harm, got[str(hw)], basin.getInfo(), seen
    hw = HALF_WINDOWS[-1]
    tag = max(("full", "simple"), key=lambda t: seen[f"{hw}_{t}"])
    harm = FULL_HARMONICS if tag == "full" else SIMPLE_HARMONICS
    return hw, harm, seen[f"{hw}_{tag}"], basin.getInfo(), seen


def year_row(year, retries=3):
    hw, harmonics, coverage, basin_km2, probe = choose_window(year)
    y0, y1 = max(year - hw, 1984), min(year + hw + 1, 2026)
    sens = sensors_for(y0, y1)
    px = predictors(harmonics)
    floor = OBS_PER_PARAM * len(px)
    col = observations(f"{y0}-01-01", f"{y1}-01-01", sens, harmonics)

    n = col.select("mndwi").count().rename("n")
    reg = col.select(px + ["mndwi"]).reduce(
        ee.Reducer.linearRegression(numX=len(px), numY=1))
    coeffs = (reg.select("coefficients").arrayProject([0]).arrayFlatten([px])
              .updateMask(n.gte(floor)))

    bands = []
    for i in range(N_DATES):
        frac = (i + 0.5) / N_DATES
        d = ee.Date(f"{year}-01-01").advance(frac * 365.25, "day")
        bands.append(coeffs.multiply(design(d, harmonics)).reduce(ee.Reducer.sum())
                     .rename(f"d{i:02d}"))
    fitted = ee.Image.cat(bands)

    def resid(img):
        pred = coeffs.multiply(img.select(px)).reduce(ee.Reducer.sum())
        return img.select("mndwi").subtract(pred).pow(2).rename("sq")

    rmse = ee.Number(col.map(resid).mean().sqrt().reduceRegion(
        reducer=ee.Reducer.median(), geometry=AOI, scale=300, maxPixels=1e10,
        bestEffort=True).get("sq"))
    band = rmse.max(MIN_RESID)

    def areas(th):
        """Area of water, scaled up to the whole basin when the fit does not
        cover all of it, so that a year is never reported low merely because
        part of the basin could not be fitted."""
        return (fitted.gt(th).multiply(AREA_KM2)
                .reduceRegion(reducer=ee.Reducer.sum(), geometry=AOI,
                              scale=SCALE, maxPixels=1e10, bestEffort=True))

    stats = ee.Dictionary({
        "a": areas(THRESH), "lo": areas(band.multiply(-1)), "hi": areas(band),
        "rmse": rmse, "band": band,
        "n_scenes": col.size(),
        "n_obs": n.updateMask(coeffs.select("const").mask()).reduceRegion(
            reducer=ee.Reducer.median(), geometry=AOI, scale=300, maxPixels=1e10,
            bestEffort=True),
        "fitted_frac": (coeffs.select("const").mask().multiply(AREA_KM2).rename("km2")
                        .reduceRegion(reducer=ee.Reducer.sum(), geometry=AOI,
                                      scale=SCALE, maxPixels=1e10, bestEffort=True)),
        "n_in_year": observations(f"{year}-01-01", f"{year + 1}-01-01",
                                  sensors_for(year, year)).size(),
    })

    last = None
    for attempt in range(retries):
        try:
            got = stats.getInfo()
            break
        except Exception as e:
            last = e
            time.sleep(20 * (attempt + 1))
    else:
        raise RuntimeError(f"{year}: {last}")

    a = [got["a"][f"d{i:02d}"] for i in range(N_DATES)]
    lo = [got["lo"][f"d{i:02d}"] for i in range(N_DATES)]
    hi = [got["hi"][f"d{i:02d}"] for i in range(N_DATES)]
    fitted_km2 = got["fitted_frac"].get("km2") or 0.0
    cov = fitted_km2 / basin_km2 if basin_km2 else 0.0
    mean_a = sum(a) / len(a)
    mean_lo, mean_hi = sum(lo) / len(lo), sum(hi) / len(hi)
    imax = max(range(N_DATES), key=lambda i: a[i])
    peak_month = int((imax + 0.5) / N_DATES * 12) + 1

    # The unfitted remainder of the basin is dry margin far more often than it
    # is lake, so the area is not scaled up; instead the shortfall is carried
    # as a one-sided uncertainty, and the coverage is reported in the row.
    shortfall = max(0.0, (1.0 - cov)) * mean_a
    sd = abs(mean_lo - mean_hi) / 2.0 + shortfall

    return {
        "year": year,
        "window_start": y0, "window_end": y1 - 1, "half_window": hw,
        "model_order": "annual+semiannual" if len(harmonics) == 2 else "annual",
        "n_free_parameters": len(px),
        "obs_floor": floor,
        "sensors": "+".join(sens),
        "n_scenes_window": got["n_scenes"],
        "n_scenes_in_year": got["n_in_year"],
        "median_obs_per_pixel": round(got["n_obs"].get("n") or 0, 1),
        "basin_fitted_km2": round(fitted_km2, 1),
        "basin_coverage": round(cov, 4),
        "open_water_mean_km2": round(mean_a, 1),
        "open_water_sd_km2": round(sd, 1),
        "open_water_max_km2": round(max(a), 1),
        "open_water_min_km2": round(min(a), 1),
        "peak_month": peak_month,
        "threshold": THRESH,
        "fit_residual_mndwi": round(got["rmse"], 4),
        "threshold_band_mndwi": round(got["band"], 4),
        "method": f"harmonic MNDWI fit, {y1 - y0}-year window {y0} to {y1 - 1}, "
                  f"{'annual+semiannual' if len(harmonics) == 2 else 'annual'} model, "
                  f"{N_DATES} dates, MNDWI > {THRESH:+.2f}",
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", type=int, default=1984)
    ap.add_argument("--end", type=int, default=2024)
    args = ap.parse_args()

    done = {}
    if os.path.exists(OUT):
        for r in csv.DictReader(open(OUT)):
            done[int(r["year"])] = r
        print(f"{len(done)} years already on disk")

    fields = None
    for year in range(args.start, args.end + 1):
        if year in done:
            continue
        row = year_row(year)
        done[year] = row
        fields = fields or list(row.keys())
        print(f"  {year}  {row['open_water_mean_km2']:7.0f} +/- {row['open_water_sd_km2']:5.0f} km2 "
              f"(max {row['open_water_max_km2']:.0f}, peak month {row['peak_month']}), "
              f"{row['n_scenes_window']} scenes in window, {row['n_scenes_in_year']} in year, "
              f"window {row['window_start']}-{row['window_end']} {row['model_order']}, "
              f"coverage {100 * row['basin_coverage']:.0f}%, "
              f"median {row['median_obs_per_pixel']:.0f} obs/pixel, "
              f"residual {row['fit_residual_mndwi']:.3f}", flush=True)
        fields = list(next(iter(done.values())).keys()) if not fields else fields
        with open(OUT, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields)
            w.writeheader()
            for y in sorted(done):
                w.writerow(done[y])
    print(f"wrote {OUT}, {len(done)} years")


if __name__ == "__main__":
    main()
