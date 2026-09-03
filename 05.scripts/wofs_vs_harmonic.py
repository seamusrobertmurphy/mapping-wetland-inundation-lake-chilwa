#!/usr/bin/env python3
"""Test the 2012 harmonic fit against Water Observations from Space.

WOfS is Digital Earth Africa's Landsat water classification. It is not an
independent observation, since its metadata names the same Landsat 7 Collection
2 Tier 1 scenes the harmonic fit is built from, but it is an independent
decision rule: a regression-tree classifier over the six reflectance bands
(Mueller et al. 2016), against this study's threshold on MNDWI. Where the two
disagree, the disagreement is about method rather than about data.

Both products are brought to a common 120 m grid in EPSG:4326. The coarser grid
is deliberate. It keeps each date to one small download, and the comparison of
interest is agreement over the basin rather than the placement of any single
30 m pixel.

For every date WOfS observed, the comparison is restricted to pixels WOfS marks
clear, so cloud and scan-line gaps are excluded from both sides alike.

Writes 03.outputs/CSV/wofs_vs_harmonic_2012.csv
"""

import csv
import io
import json
import os
import subprocess
import urllib.request
import zipfile

import numpy as np
import rasterio
import requests
import ee

ee.Initialize(project="murphys-deforisk")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRATCH = os.path.join(ROOT, "03.outputs", "TIF", "wofs_2012")
os.makedirs(SCRATCH, exist_ok=True)

SCALE = 120                 # metres, the common grid
MNDWI_THRESHOLD = 0.0
STAC = "https://explorer.digitalearth.africa/stac/search"
S3 = ("s3://deafrica-services/", "https://deafrica-services.s3.af-south-1.amazonaws.com/")
REACHING = ("167070", "167071", "166071")   # path/rows with pixels inside the basin

gj = json.load(open(os.path.join(ROOT, "03.outputs", "JSON", "chilwa_basin.geojson")))
GEOM = gj["features"][0]["geometry"] if gj.get("type") == "FeatureCollection" else gj["geometry"]
AOI = ee.Geometry(GEOM)


def wofs_items(start="2012-01-01", end="2013-01-01"):
    """WOfS water rasters over the basin, grouped by date."""
    body = json.dumps({"collections": ["wofs_ls"], "intersects": GEOM,
                       "datetime": f"{start}/{end}", "limit": 500}).encode()
    req = urllib.request.Request(STAC, data=body, headers={
        "Content-Type": "application/json", "Accept": "application/geo+json"})
    d = json.load(urllib.request.urlopen(req, timeout=180))
    out = {}
    for f in d["features"]:
        p = f["properties"]
        if str(p.get("odc:region_code")) not in REACHING:
            continue
        out.setdefault(p["datetime"][:10], []).append(
            f["assets"]["water"]["href"].replace(*S3))
    return dict(sorted(out.items()))


def fetch_wofs(date, hrefs, grid):
    """Download the date's WOfS tiles and warp them onto the reference grid."""
    parts = []
    for i, h in enumerate(hrefs):
        p = os.path.join(SCRATCH, f"wofs_{date}_{i}.tif")
        if not os.path.exists(p):
            r = requests.get(h, timeout=600)
            r.raise_for_status()
            open(p, "wb").write(r.content)
        parts.append(p)
    out = os.path.join(SCRATCH, f"wofs_{date}_grid.tif")
    subprocess.run(
        ["gdalwarp", "-overwrite", "-t_srs", "EPSG:4326", "-r", "near",
         "-te", *[str(v) for v in grid["bounds"]],
         "-ts", str(grid["width"]), str(grid["height"]),
         "-srcnodata", "1", "-dstnodata", "1", *parts, out],
        check=True, capture_output=True)
    with rasterio.open(out) as src:
        return src.read(1)


def fitted_mndwi(coeffs, date, grid):
    """The harmonic fit on one date, pulled at the reference grid's scale."""
    import harmonic_fit_2012 as h
    img = h.evaluate(coeffs, date, trend=True).multiply(10000).toInt16().clip(AOI)
    url = img.getDownloadURL({"scale": SCALE, "region": grid["region"],
                              "filePerBand": False, "format": "ZIPPED_GEO_TIFF"})
    r = requests.get(url, timeout=900)
    if not r.ok:
        raise RuntimeError(f"{date}: {r.status_code} {r.text[:300]}")
    p = os.path.join(SCRATCH, f"fit_{date}.tif")
    with zipfile.ZipFile(io.BytesIO(r.content)) as z:
        open(p, "wb").write(z.read(z.namelist()[0]))
    with rasterio.open(p) as src:
        a = src.read(1).astype("float32")
        return np.where(a == src.nodata if src.nodata is not None else False,
                        np.nan, a / 10000.0), src


def main():
    import sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import harmonic_fit_2012 as h

    dates = wofs_items()
    print(f"{len(dates)} WOfS dates over the basin in 2012\n")

    _, coeffs, _ = h.build()

    # The reference grid comes from one Earth Engine pull, so both products land
    # on exactly the same pixels.
    probe = ee.Image.constant(0).toInt16().clip(AOI)
    region = AOI.bounds().getInfo()["coordinates"]
    url = probe.getDownloadURL({"scale": SCALE, "region": region,
                                "filePerBand": False, "format": "ZIPPED_GEO_TIFF"})
    p = os.path.join(SCRATCH, "grid.tif")
    with zipfile.ZipFile(io.BytesIO(requests.get(url, timeout=600).content)) as z:
        open(p, "wb").write(z.read(z.namelist()[0]))
    with rasterio.open(p) as src:
        t = src.transform
        grid = {"bounds": (t.c, t.f + src.height * t.e, t.c + src.width * t.a, t.f),
                "width": src.width, "height": src.height, "region": region}
        rows = np.arange(src.height)
        lat = t.f + (rows + 0.5) * t.e
        px_km2 = (abs(t.a) * 111320.0 * np.cos(np.radians(lat)) * abs(t.e) * 110540.0) / 1e6
    print(f"reference grid {grid['width']} x {grid['height']} at {SCALE} m\n")

    out = []
    for date, hrefs in dates.items():
        try:
            w = fetch_wofs(date, hrefs, grid)
            fit, _ = fitted_mndwi(coeffs, date, grid)
        except Exception as e:
            print(f"  {date}  skipped, {type(e).__name__} {str(e)[:80]}")
            continue

        # WOfS bit 7 is water; any other flag set means cloud, shadow, slope or gap.
        wofs_wet = w == 128
        wofs_dry = w == 0
        clear = (wofs_wet | wofs_dry) & np.isfinite(fit)
        if clear.sum() == 0:
            print(f"  {date}  no co-valid pixels")
            continue

        fit_wet = fit > MNDWI_THRESHOLD
        area = px_km2[:, None]
        a_wofs = float(((wofs_wet & clear) * area).sum())
        a_fit = float(((fit_wet & clear) * area).sum())
        a_clear = float((clear * area).sum())
        agree = float(((wofs_wet == fit_wet) & clear).sum()) / clear.sum()
        both = float(((wofs_wet & fit_wet) & clear).sum())
        union = float(((wofs_wet | fit_wet) & clear).sum())
        iou = both / union if union else float("nan")

        print(f"  {date}  clear {a_clear:6,.0f} km2 | WOfS {a_wofs:6,.0f} | "
              f"fit {a_fit:6,.0f} | agree {100*agree:5.1f} % | IoU {iou:.3f}")
        out.append({"date": date, "clear_km2": round(a_clear, 1),
                    "wofs_water_km2": round(a_wofs, 1), "fit_water_km2": round(a_fit, 1),
                    "pixel_agreement": round(agree, 4), "water_iou": round(iou, 4),
                    "n_wofs_tiles": len(hrefs)})

    with open(os.path.join(ROOT, "03.outputs", "CSV", "wofs_vs_harmonic_2012.csv"),
              "w", newline="") as fh:
        wcsv = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
        wcsv.writeheader()
        wcsv.writerows(out)
    print(f"\nwrote 03.outputs/CSV/wofs_vs_harmonic_2012.csv, {len(out)} dates")


if __name__ == "__main__":
    main()
