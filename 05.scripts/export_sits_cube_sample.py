#!/usr/bin/env python3
"""Download a sample of the Lake Chilwa time series as a local data cube.

Purpose. Everything upstream of this runs in Earth Engine and returns numbers.
This returns pixels, on a regular sixteen-day grid, named so that the sits
package in R can read them as a cube without any conversion step. It is the
handover point between the Earth Engine stage of the analysis and the local one.

Scope, and why 2018 rather than the fieldwork year. The point of this download
is to prove that a cube assembled from instruments with different grains, orbits
and physics converts into something sits can read with metadata it recognises.
2012 cannot test that, because Landsat 7 is the only sensor over this basin that
year and a single-source stack proves nothing about a multisource one. 2018
carries Landsat 7 and 8 at thirty metres, Sentinel-2 at ten, Sentinel-1 at ten,
ALOS PALSAR at twenty-five and the Moderate Resolution Imaging Spectroradiometer
at five hundred, which is every family the final model is meant to fit. It is
also the deepest drawdown in the whole record, falling to 48 square kilometres,
so the test period is the one where the instruments most disagree and the
conversion is hardest.

Everything is resampled onto one thirty-metre grid at sixteen-day steps. That
resampling is the normalisation the study has been calling for, and it is done
here rather than downstream so that what lands on disk is already a cube rather
than a pile of rasters on different grids.

Bands. Optical carries green, near infrared, shortwave infrared and the water
index. Radar carries the two polarisations each instrument transmits, since the
difference between them is what reveals water beneath vegetation. The band names
say which instrument they came from, because a cube that hides its provenance
cannot be audited.

Gaps are preserved, not filled. Where a sixteen-day window holds no clear
observation the pixels are written as nodata rather than interpolated, because
the fusion that fills them lives downstream and a cube that silently
interpolates cannot be told apart from one that measures.

Naming follows the sits local-cube convention, five underscore-delimited fields
read by parse_info = c("X1", "X2", "tile", "band", "date"), so the files are
CHILWA_LANDSAT_lake_<BAND>_<YYYY-MM-DD>.tif

Writes 02.inputs/CUBE/chilwa_2012_30m/
"""

import argparse
import datetime
import io
import json
import os
import subprocess
import time
import zipfile

import ee
import geopandas as gpd
import requests

ee.Initialize(project="murphys-deforisk")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "02.inputs", "CUBE", "chilwa_multisource_30m")
os.makedirs(OUT, exist_ok=True)

_lk = gpd.read_file(os.path.join(ROOT, "02.inputs", "SHP", "lakes_site.shp")).to_crs(4326)
LAKE = ee.Geometry(json.loads(json.dumps(_lk.geometry.union_all().__geo_interface__)))
AOI = LAKE.buffer(10000).bounds()
SCALE = 30
STEP = 16

# Landsat 7 is the only optical sensor over this basin in 2012, which is the
# constraint the whole study has been working around. Landsat 5 stopped imaging
# here in October 2011 and Landsat 8 began in March 2013.
LS = {"LE07": ("SR_B2", "SR_B4", "SR_B5"), "LT05": ("SR_B2", "SR_B4", "SR_B5"),
      "LC08": ("SR_B3", "SR_B5", "SR_B6"), "LC09": ("SR_B3", "SR_B5", "SR_B6")}


def clear(img):
    qa = img.select("QA_PIXEL")
    bad = (qa.bitwiseAnd(1 << 1).Or(qa.bitwiseAnd(1 << 2)).Or(qa.bitwiseAnd(1 << 3))
           .Or(qa.bitwiseAnd(1 << 4)).Or(qa.bitwiseAnd(1 << 5)))
    return img.updateMask(bad.eq(0))


def landsat(t0, t1):
    """Landsat at its native thirty metres, the reference grid."""
    out = None
    for sensor, (g, n, w) in LS.items():
        col = (ee.ImageCollection(f"LANDSAT/{sensor}/C02/T1_L2").filterBounds(AOI)
               .filterDate(t0.isoformat(), t1.isoformat()).map(clear))

        def prep(img, g=g, n=n, w=w):
            sr = img.select("SR_B.").multiply(0.0000275).add(-0.2)
            G = sr.select(g).rename("LSGREEN")
            N = sr.select(n).rename("LSNIR")
            W = sr.select(w).rename("LSSWIR1")
            M = G.subtract(W).divide(G.add(W)).rename("LSMNDWI")
            return G.addBands([N, W, M])

        col = col.map(prep)
        out = col if out is None else out.merge(col)
    return out


def sentinel2(t0, t1):
    """Sentinel-2 at ten metres, resampled to the Landsat grid on export."""
    col = (ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED").filterBounds(AOI)
           .filterDate(t0.isoformat(), t1.isoformat())
           .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", 60)))

    def prep(img):
        scl = img.select("SCL")
        good = scl.eq(4).Or(scl.eq(5)).Or(scl.eq(6)).Or(scl.eq(7)).Or(scl.eq(11))
        sr = img.updateMask(good).divide(10000)
        G = sr.select("B3").rename("S2GREEN")
        N = sr.select("B8").rename("S2NIR")
        W = sr.select("B11").rename("S2SWIR1")
        M = G.subtract(W).divide(G.add(W)).rename("S2MNDWI")
        return G.addBands([N, W, M])
    return col.map(prep)


def sentinel1(t0, t1):
    """Sentinel-1 backscatter, vertical transmit, already terrain corrected."""
    col = (ee.ImageCollection("COPERNICUS/S1_GRD").filterBounds(AOI)
           .filterDate(t0.isoformat(), t1.isoformat())
           .filter(ee.Filter.eq("instrumentMode", "IW"))
           .filter(ee.Filter.listContains("transmitterReceiverPolarisation", "VV")))
    return col.map(lambda im: im.select("VV").rename("S1VV")
                   .addBands(im.select("VH").rename("S1VH")))


def modis(t0, t1):
    """MODIS at five hundred metres, the dense coarse stream."""
    col = (ee.ImageCollection("MODIS/061/MOD09A1").filterBounds(AOI)
           .filterDate(t0.isoformat(), t1.isoformat()))

    def prep(img):
        qa = img.select("StateQA")
        cloud = qa.bitwiseAnd(3).neq(0).Or(qa.bitwiseAnd(1 << 2).neq(0))
        g = img.select("sur_refl_b04").multiply(0.0001)
        w = img.select("sur_refl_b06").multiply(0.0001)
        return (g.subtract(w).divide(g.add(w)).rename("MODMNDWI")
                .updateMask(cloud.Not()))
    return col.map(prep)


def palsar(year):
    """ALOS PALSAR yearly mosaic, horizontal transmit, one image per year."""
    col = (ee.ImageCollection("JAXA/ALOS/PALSAR/YEARLY/SAR_EPOCH")
           .filterBounds(AOI).filterDate(f"{year}-01-01", f"{year + 1}-01-01"))
    if col.size().getInfo() == 0:
        return None
    im = ee.Image(col.mosaic())
    hh = im.select("HH").pow(2).log10().multiply(10).add(-83.0).rename("PSHH")
    hv = im.select("HV").pow(2).log10().multiply(10).add(-83.0).rename("PSHV")
    return hh.addBands(hv)


def download(image, path, tiles=2):
    b = AOI.bounds().getInfo()["coordinates"][0]
    xs = [c[0] for c in b]; ys = [c[1] for c in b]
    x0, x1, y0, y1 = min(xs), max(xs), min(ys), max(ys)
    parts = []
    for i in range(tiles):
        for j in range(tiles):
            box = ee.Geometry.Rectangle([
                x0 + (x1 - x0) * i / tiles, y0 + (y1 - y0) * j / tiles,
                x0 + (x1 - x0) * (i + 1) / tiles, y0 + (y1 - y0) * (j + 1) / tiles])
            p = path.replace(".tif", f"_{i}{j}.tif")
            # Earth Engine returns 429 and 503 under load. These are transient
            # and the whole download must not fail on one of them, which is what
            # an earlier version did after eighteen of twenty-three steps.
            last = None
            for attempt in range(6):
                try:
                    url = image.getDownloadURL({"scale": SCALE,
                                                "region": box.getInfo()["coordinates"],
                                                "filePerBand": False,
                                                "format": "ZIPPED_GEO_TIFF"})
                    r = requests.get(url, timeout=900)
                    if r.ok:
                        with zipfile.ZipFile(io.BytesIO(r.content)) as z:
                            open(p, "wb").write(z.read(z.namelist()[0]))
                        break
                    last = f"HTTP {r.status_code}"
                except Exception as e:
                    last = str(e)[:90]
                time.sleep(min(20 * (attempt + 1), 90))
            else:
                raise RuntimeError(f"{os.path.basename(path)}: gave up after six tries, {last}")
            parts.append(p)
    vrt = path.replace(".tif", ".vrt")
    subprocess.run(["gdalbuildvrt", "-srcnodata", "-32768", "-vrtnodata", "-32768",
                    vrt] + parts, check=True, capture_output=True)
    subprocess.run(["gdal_translate", "-of", "GTiff", "-co", "COMPRESS=DEFLATE",
                    "-co", "TILED=YES", "-a_nodata", "-32768", vrt, path],
                   check=True, capture_output=True)
    for f in parts + [vrt]:
        os.remove(f)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--year", type=int, default=2018)
    a = ap.parse_args()
    start = datetime.date(a.year, 1, 1)
    n = (datetime.date(a.year, 12, 31) - start).days // STEP + 1
    steps = [start + datetime.timedelta(days=STEP * i) for i in range(n)]
    print(f"{a.year}: {len(steps)} sixteen-day steps on one {SCALE} m grid, "
          f"lake plus a ten-kilometre belt", flush=True)

    # PALSAR is annual, so its two bands repeat across the year's steps and the
    # cube records that by writing them once per step with the same values. A
    # reader must not mistake that for sixteen-day radar.
    ps = palsar(a.year)
    if ps is not None:
        print("  ALOS PALSAR: one annual mosaic, written to every step", flush=True)

    written = 0
    for s_ in steps:
        t1 = s_ + datetime.timedelta(days=STEP)
        pieces = []
        for name, col in (("Landsat", landsat(s_, t1)), ("Sentinel-2", sentinel2(s_, t1)),
                          ("Sentinel-1", sentinel1(s_, t1)), ("MODIS", modis(s_, t1))):
            k = col.size().getInfo()
            if k:
                pieces.append((name, k, col.median()))
        if not pieces:
            print(f"  {s_}: nothing observed by any instrument", flush=True)
            continue
        img = pieces[0][2]
        for _, _, m in pieces[1:]:
            img = img.addBands(m)
        if ps is not None:
            img = img.addBands(ps)
        names = img.bandNames().getInfo()
        for band in names:
            path = os.path.join(OUT, f"CHILWA_MULTI_lake_{band}_{s_.isoformat()}.tif")
            if os.path.exists(path) and os.path.getsize(path) > 0:
                continue
            layer = img.select(band)
            # radar is decibels, reflectance and indices are unitless; both are
            # stored as int16 times 10000 so one nodata value serves the cube
            download(layer.multiply(10000).toInt16().unmask(-32768), path)
            written += 1
        print(f"  {s_}: " + ", ".join(f"{nm} {k}" for nm, k, _ in pieces)
              + f" -> {len(names)} bands", flush=True)

    print(f"\n{written} files in {OUT}")
    print("read in R with parse_info = c('X1','X2','tile','band','date'), delim = '_'")


if __name__ == "__main__":
    main()
