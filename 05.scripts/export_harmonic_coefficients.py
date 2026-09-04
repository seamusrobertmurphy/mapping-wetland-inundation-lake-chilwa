#!/usr/bin/env python3
"""Export the six fitted harmonic coefficients so the optical model can be read
on any date without a live Earth Engine session.

harmonic_fit_2012.py fits, at every 30 m pixel, MNDWI as a level, a linear
drift and annual plus semi-annual sine and cosine pairs, on Landsat 5, 7 and 8
between June 2010 and June 2014. That script evaluates the model on one date
and writes that one raster. The Envisat ASAR series holds sixteen dates between
November 2011 and March 2012, and the plan is to read the optical model on each
radar date exactly rather than resample either product onto the other's
calendar. That needs the coefficients themselves.

Writes
  03.outputs/TIF/harmonic_2012/harmonic_coefficients.tif
      six bands, const, t, cos1, sin1, cos2, sin2, int16 x 10000, nodata -32768
  03.outputs/TIF/harmonic_2012/harmonic_n_obs.tif
      clear observations behind each pixel, int16
  03.outputs/CSV/harmonic_date_terms.csv
      for every ASAR date and the 22 February target, the value of t in years
      since 2012-01-01 exactly as Earth Engine computes it, and the cosine and
      sine terms, so a local evaluation reproduces the server one.

Reading the model locally on date d, with t from the CSV:
  MNDWI = const + t*t_coef + cos1*cos(2 pi t) + sin1*sin(2 pi t)
                + cos2*cos(4 pi t) + sin2*sin(4 pi t)
"""

import csv
import io
import math
import os
import subprocess
import time
import zipfile

import ee
import requests

import harmonic_fit_2012 as hf

ROOT = hf.ROOT
OUT_TIF = hf.OUT_TIF
OUT_CSV = hf.OUT_CSV

ASSET = "projects/murphys-deforisk/assets/harmonic_coefficients_2010_2014"

ASAR_DATES = [
    "2011-11-05", "2011-11-13", "2011-11-16", "2011-12-05", "2011-12-16",
    "2011-12-24", "2012-01-04", "2012-01-12", "2012-01-15", "2012-02-03",
    "2012-02-11", "2012-02-14", "2012-02-22", "2012-03-04", "2012-03-12",
    "2012-03-15",
]


def export(image, name, factor, scale=30, tiles=3, nodata=-32768):
    """Tiled direct download, merged with gdalbuildvrt, as in harmonic_fit_2012."""
    img = image.multiply(factor).toInt16().clip(hf.AOI)
    b = hf.AOI.bounds().getInfo()["coordinates"][0]
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
                continue
            url = img.getDownloadURL(
                {"scale": scale, "region": box.getInfo()["coordinates"],
                 "filePerBand": False, "format": "ZIPPED_GEO_TIFF"})
            r = requests.get(url, timeout=900)
            if not r.ok:
                raise RuntimeError(f"tile {i}{j}: {r.status_code} {r.text[:400]}")
            with zipfile.ZipFile(io.BytesIO(r.content)) as z:
                with open(part, "wb") as fh:
                    fh.write(z.read(z.namelist()[0]))
            parts.append(part)
            print(f"  {name} tile {i}{j}  {os.path.getsize(part) / 1e6:5.1f} MB",
                  flush=True)
    out = os.path.join(OUT_TIF, f"{name}.tif")
    vrt = out.replace(".tif", ".vrt")
    subprocess.run(["gdalbuildvrt", "-srcnodata", str(nodata), "-vrtnodata",
                    str(nodata), vrt] + parts, check=True, capture_output=True)
    subprocess.run(["gdal_translate", "-of", "GTiff", "-co", "COMPRESS=DEFLATE",
                    "-co", "TILED=YES", "-a_nodata", str(nodata), vrt, out],
                   check=True, capture_output=True)
    for f in parts + [vrt]:
        os.remove(f)
    print(f"wrote {out}  ({os.path.getsize(out) / 1e6:.1f} MB)", flush=True)
    return out


def main():
    col, coeffs, n = hf.build()
    names = hf.predictors(trend=True)
    print("coefficient bands:", names, flush=True)

    with open(os.path.join(OUT_CSV, "harmonic_date_terms.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["date", "t_years"] + names[2:])
        for d in ASAR_DATES + [hf.TARGET]:
            t = ee.Date(d).difference(hf.ORIGIN, "year").getInfo()
            row = [d, f"{t:.10f}"]
            for k in hf.HARMONICS:
                row += [f"{math.cos(2 * math.pi * k * t):.10f}",
                        f"{math.sin(2 * math.pi * k * t):.10f}"]
            w.writerow(row)
            print(f"  {d}  t = {t:.6f} yr", flush=True)

    # Set the band order explicitly and check the 22 February prediction from
    # the coefficients matches the raster already on disk, on the server, before
    # anything is downloaded.
    coeffs = coeffs.select(names)
    pred = hf.evaluate(coeffs, hf.TARGET, trend=True)
    area = hf.water_area(pred, 0.0)
    print(f"server check, open water on {hf.TARGET} at MNDWI > 0: {area:,.0f} km2",
          flush=True)

    # The direct download route recomputes the regression over 243 scenes for
    # every tile and exceeds the interactive memory limit. The batch system has
    # no such limit, so the fitted image is stored once as an asset and the
    # tiles are then read from the stored copy, which involves no computation.
    try:
        ee.data.getAsset(ASSET)
        print(f"asset {ASSET} already exists", flush=True)
    except ee.EEException:
        task = ee.batch.Export.image.toAsset(
            image=coeffs.addBands(n.select("n")).toFloat().clip(hf.AOI),
            description="harmonic_coefficients_2010_2014", assetId=ASSET,
            region=hf.AOI.bounds(), scale=30, crs="EPSG:4326", maxPixels=1e9)
        task.start()
        print(f"asset export started, task {task.id}", flush=True)
        while task.active():
            time.sleep(30)
            print(f"  {task.status()['state']}", flush=True)
        state = task.status()
        if state["state"] != "COMPLETED":
            raise SystemExit(f"asset export failed: {state}")
        print("asset export complete", flush=True)

    stored = ee.Image(ASSET)
    export(stored.select(names), "harmonic_coefficients", factor=10000)
    export(stored.select("n"), "harmonic_n_obs", factor=1)


if __name__ == "__main__":
    main()
