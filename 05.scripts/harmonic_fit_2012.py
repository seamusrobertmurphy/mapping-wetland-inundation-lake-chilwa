#!/usr/bin/env python3
"""Harmonic fit of the 2012 water signal, to date the Lake Chilwa recession.

Why a harmonic fit rather than LandTrendr. LandTrendr reduces each year to one
value and joins those values with straight lines, which suits forests that
change slowly in one direction. Lake Chilwa rises and falls every year, so an
annual straight-line fit would remove the flood pulse that is being measured.
A harmonic fit keeps it, by modelling the repeating shape of the year.

Model, fitted independently at every 30 m pixel:

    MNDWI(t) = a0 + a1*t + b1*cos(2*pi*t) + c1*sin(2*pi*t)
                        + b2*cos(4*pi*t) + c2*sin(4*pi*t)

with t in years. The annual pair carries the single wet peak, the semi-annual
pair lets the rise and the recession differ in steepness, and the linear term
absorbs drift across the window. Six coefficients are fitted against the clear
observations available at that pixel, and the model is then evaluated on
22 February 2012, the date of the 43-photograph flooded-vegetation session,
for which the direct imagery covers only 22 per cent of the basin.

Landsat 7 supplies the observations. Its SLC-off gaps move between passes, so
stacking the 2012 scenes leaves no part of the basin unobserved; this was
measured in 05.scripts/l7_2012_coverage_test.py.

Reports the fit residual, the number of observations behind every pixel, and a
holdout test in which whole dates are withheld and predicted.

Writes 03.outputs/TIF/harmonic_2012/ and 03.outputs/CSV/harmonic_fit_2012_*.csv
"""

import csv
import io
import json
import math
import os
import zipfile

import ee
import requests

ee.Initialize(project="murphys-deforisk")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_TIF = os.path.join(ROOT, "03.outputs", "TIF", "harmonic_2012")
OUT_CSV = os.path.join(ROOT, "03.outputs", "CSV")
os.makedirs(OUT_TIF, exist_ok=True)

gj = json.load(open(os.path.join(ROOT, "03.outputs", "JSON", "chilwa_basin.geojson")))
geom = gj["features"][0]["geometry"] if gj.get("type") == "FeatureCollection" else gj
AOI = ee.Geometry(geom)

TARGET = "2012-02-22"          # the 43-photograph flooded-vegetation session
ORIGIN = ee.Date("2012-01-01")  # t = 0
HARMONICS = [1, 2]              # annual and semi-annual
MIN_OBS = 12                    # at least twice the free parameters


def predictors(trend):
    """Design columns. The linear term is near-collinear with the annual pair
    over a window of exactly one period, so it is carried only across years."""
    return (["const"] + (["t"] if trend else [])
            + [f"{f}{k}" for k in HARMONICS for f in ("cos", "sin")])
MNDWI_THRESHOLDS = (-0.10, 0.00, 0.10)

AREA_KM2 = ee.Image.pixelArea().divide(1e6)


def scaled(img):
    """Collection 2 Level-2 digital numbers to surface reflectance."""
    return ee.Image(img.select("SR_B.").multiply(0.0000275).add(-0.2)
                    .copyProperties(img, ["system:time_start"]))


def clear(img):
    """Drop cloud, shadow, cirrus and snow. SLC-off gaps are already unobserved."""
    qa = img.select("QA_PIXEL")
    bad = (qa.bitwiseAnd(1 << 1).Or(qa.bitwiseAnd(1 << 2)).Or(qa.bitwiseAnd(1 << 3))
             .Or(qa.bitwiseAnd(1 << 4)).Or(qa.bitwiseAnd(1 << 5)))
    return img.updateMask(bad.eq(0))


def design(date, trend):
    """The predictor values at one instant."""
    t = ee.Number(ee.Date(date).difference(ORIGIN, "year"))
    bands = [ee.Image.constant(1)] + ([ee.Image.constant(t)] if trend else [])
    for k in HARMONICS:
        w = t.multiply(2 * math.pi * k)
        bands += [ee.Image.constant(w.cos()), ee.Image.constant(w.sin())]
    return ee.Image.cat(bands).rename(predictors(trend)).float()


def observations(start, end, sensors=("LE07",), trend=False):
    """Clear MNDWI observations with their design row attached.

    Band positions differ between the Thematic Mapper family and the Operational
    Land Imager, so green and SWIR1 are named per sensor rather than assumed."""
    BANDS = {"LT05": ("SR_B2", "SR_B5"), "LE07": ("SR_B2", "SR_B5"),
             "LC08": ("SR_B3", "SR_B6"), "LC09": ("SR_B3", "SR_B6")}
    out = None
    for sensor in sensors:
        green, swir1 = BANDS[sensor]
        col = (ee.ImageCollection(f"LANDSAT/{sensor}/C02/T1_L2")
               .filterBounds(AOI).filterDate(start, end).map(clear))

        def prep(img, green=green, swir1=swir1):
            d = ee.Date(img.get("system:time_start"))
            mndwi = scaled(img).normalizedDifference([green, swir1]).rename("mndwi")
            return (design(d, trend).addBands(mndwi)
                    .updateMask(mndwi.mask())
                    .set("system:time_start", img.get("system:time_start")))

        col = col.map(prep)
        out = col if out is None else out.merge(col)
    return out.sort("system:time_start")


def fit(col, trend=False, min_obs=MIN_OBS):
    """Least squares at every pixel.

    Pixels carrying fewer than min_obs clear observations are dropped. Without
    that mask a pixel with as many observations as free parameters fits exactly,
    reporting zero residual and extrapolating far outside the -1 to 1 range that
    MNDWI can occupy."""
    px = predictors(trend)
    reg = col.select(px + ["mndwi"]).reduce(
        ee.Reducer.linearRegression(numX=len(px), numY=1))
    coeffs = reg.select("coefficients").arrayProject([0]).arrayFlatten([px])
    n = col.select("mndwi").count().rename("n")
    return coeffs.updateMask(n.gte(min_obs)), n


def evaluate(coeffs, date, trend=False):
    """Fitted MNDWI on one date."""
    return (coeffs.multiply(design(date, trend))
            .reduce(ee.Reducer.sum()).rename("mndwi_fit"))


def residual_rmse(coeffs, col, trend=False):
    """Root mean square difference between fitted and observed.

    In-fit only, so it flatters the model; use holdout() for an honest figure."""
    px = predictors(trend)

    def resid(img):
        pred = coeffs.multiply(img.select(px)).reduce(ee.Reducer.sum())
        return img.select("mndwi").subtract(pred).pow(2).rename("se")

    return ee.ImageCollection(col.map(resid)).mean().sqrt().rename("rmse")


def holdout(col, date, trend=False, min_obs=MIN_OBS):
    """Refit with one whole date withheld, then predict it.

    Withholding the date rather than scattered pixels is the test that matters
    here, because the question is whether the model can stand in for a day on
    which the sensor saw little."""
    d = ee.Date(date)
    kept = col.filter(ee.Filter.date(d.advance(-1, "day"), d.advance(1, "day")).Not())
    held = col.filter(ee.Filter.date(d.advance(-1, "day"), d.advance(1, "day")))
    coeffs, _ = fit(kept, trend, min_obs)
    px = predictors(trend)

    def resid(img):
        pred = coeffs.multiply(img.select(px)).reduce(ee.Reducer.sum())
        return img.select("mndwi").subtract(pred).rename("err")

    err = ee.ImageCollection(held.map(resid)).mean()
    return err


def basin_stat(img, reducer=ee.Reducer.mean(), scale=30):
    return img.reduceRegion(reducer, AOI, scale, maxPixels=1e10).getInfo()


def water_area(mndwi_img, threshold):
    m = mndwi_img.gt(threshold)
    v = AREA_KM2.updateMask(m).reduceRegion(ee.Reducer.sum(), AOI, 30, maxPixels=1e10)
    return v.get("area").getInfo()


# ---------------------------------------------------------------------------
# The window. A fit inside 2012 alone was tried first and abandoned: the median
# pixel carries 8 clear Landsat 7 observations that year, against 5 free
# parameters, and requiring a safe 12 observations left only 20.7 per cent of
# the basin. Widening to the three years either side brings in Landsat 5 (to
# October 2011) and Landsat 8 (from March 2013) and raises the count enough to
# fit, at the cost of assuming the shape of the year is stable across the
# window. The linear term absorbs the drift the recession puts into that
# assumption, and holdout() is what tests whether the assumption survives.
# ---------------------------------------------------------------------------

WINDOW = ("2010-06-01", "2014-06-01")
SENSORS = ("LT05", "LE07", "LC08")


def build():
    col = observations(*WINDOW, sensors=SENSORS, trend=True)
    coeffs, n = fit(col, trend=True)
    return col, coeffs, n


def report(scale=60):
    col, coeffs, n = build()
    print(f"window {WINDOW[0]} to {WINDOW[1]}, sensors {', '.join(SENSORS)}")
    print(f"scenes {col.size().getInfo()}, statistics at {scale} m\n")

    tot = AREA_KM2.updateMask(ee.Image(1).clip(AOI)).reduceRegion(
        ee.Reducer.sum(), AOI, scale, maxPixels=1e10).get("area").getInfo()
    print("clear observations per pixel")
    print("  percentiles", n.clip(AOI).reduceRegion(
        ee.Reducer.percentile([5, 25, 50, 75, 95]), AOI, scale, maxPixels=1e10).getInfo())
    fitted_area = AREA_KM2.updateMask(coeffs.select("const").mask()).reduceRegion(
        ee.Reducer.sum(), AOI, scale, maxPixels=1e10).get("area").getInfo()
    print(f"  pixels meeting the {MIN_OBS}-observation floor: "
          f"{fitted_area:,.0f} km2, {100 * fitted_area / tot:.1f} % of basin\n")

    img = evaluate(coeffs, TARGET, trend=True)
    print(f"fitted MNDWI on {TARGET}",
          img.reduceRegion(ee.Reducer.percentile([5, 50, 95]), AOI, scale,
                           maxPixels=1e10).getInfo())
    rmse = residual_rmse(coeffs, col, trend=True)
    print("in-fit RMSE (flatters the model)",
          rmse.reduceRegion(ee.Reducer.mean(), AOI, scale, maxPixels=1e10).getInfo())
    return col, coeffs, img


def validate(col, dates, scale=60):
    """Withhold each date in turn, predict it, and report the error."""
    print("\nholdout, whole date withheld then predicted")
    for d in dates:
        err = holdout(col, d, trend=True)
        st = err.abs().reduceRegion(
            ee.Reducer.mean().combine(ee.Reducer.percentile([50, 90]), sharedInputs=True),
            AOI, scale, maxPixels=1e10).getInfo()
        print(f"  {d}  mean |error| {st.get('err_mean'):.4f} MNDWI, "
              f"median {st.get('err_p50'):.4f}, p90 {st.get('err_p90'):.4f}")


def export(image, name, scale=30, tiles=2):
    """Pull a GeoTIFF down in tiles and merge them.

    Values are stored as int16 MNDWI x 10000, because a float band is larger
    than the 50 MB direct-download ceiling; at 30 m even int16 over the whole
    basin is 56 MB, hence the tiling."""
    import subprocess

    scaled_img = image.multiply(10000).toInt16().clip(AOI)
    b = AOI.bounds().getInfo()["coordinates"][0]
    xs = [c[0] for c in b]
    ys = [c[1] for c in b]
    x0, x1, y0, y1 = min(xs), max(xs), min(ys), max(ys)
    parts = []
    for i in range(tiles):
        for j in range(tiles):
            box = ee.Geometry.Rectangle([
                x0 + (x1 - x0) * i / tiles, y0 + (y1 - y0) * j / tiles,
                x0 + (x1 - x0) * (i + 1) / tiles, y0 + (y1 - y0) * (j + 1) / tiles])
            part = os.path.join(OUT_TIF, f"{name}_{i}{j}.tif")
            if os.path.exists(part) and os.path.getsize(part) > 0:
                parts.append(part)
                print(f"  tile {i}{j}  {os.path.getsize(part) / 1e6:5.1f} MB (on disk)")
                continue
            url = scaled_img.getDownloadURL(
                {"scale": scale, "region": box.getInfo()["coordinates"],
                 "filePerBand": False, "format": "ZIPPED_GEO_TIFF"})
            r = requests.get(url, timeout=900)
            if not r.ok:
                raise RuntimeError(f"tile {i}{j}: {r.status_code} {r.text[:400]}")
            with zipfile.ZipFile(io.BytesIO(r.content)) as z:
                with open(part, "wb") as fh:
                    fh.write(z.read(z.namelist()[0]))
            parts.append(part)
            print(f"  tile {i}{j}  {os.path.getsize(part) / 1e6:5.1f} MB")

    out = os.path.join(OUT_TIF, f"{name}.tif")
    vrt = out.replace(".tif", ".vrt")
    subprocess.run(["gdalbuildvrt", "-srcnodata", "-32768", "-vrtnodata", "-32768",
                    vrt] + parts, check=True, capture_output=True)
    subprocess.run(["gdal_translate", "-of", "GTiff", "-co", "COMPRESS=DEFLATE",
                    "-co", "TILED=YES", "-a_nodata", "-32768", vrt, out],
                   check=True, capture_output=True)
    for f in parts + [vrt]:
        os.remove(f)
    print(f"wrote {out}  ({os.path.getsize(out) / 1e6:.1f} MB at {scale} m, "
          f"values are MNDWI x 10000, nodata -32768)")
    return out


def main():
    col, coeffs, img = report()
    validate(col, ["2012-02-15", "2012-01-30", "2012-04-19"])
    print()
    with open(os.path.join(OUT_CSV, "harmonic_fit_2012_water_area.csv"), "w",
              newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["date", "mndwi_threshold", "open_water_km2", "method"])
        for t in MNDWI_THRESHOLDS:
            a = water_area(img, t)
            print(f"fitted open water on {TARGET} at MNDWI > {t:+.2f}: {a:8,.0f} km2")
            w.writerow([TARGET, f"{t:+.2f}", round(a, 1),
                        "harmonic fit, Landsat 5/7/8, 2010-06 to 2014-06"])
    export(img, f"mndwi_harmonic_{TARGET}")


if __name__ == "__main__":
    main()
