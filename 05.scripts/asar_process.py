#!/usr/bin/env python3
"""Envisat ASAR Alternating Polarisation Precision (ASA_APP_1P) to calibrated
backscatter on the grid of the 2012 optical harmonic fit.

No SNAP. GDAL reads the Envisat N1 format directly, including the tie-point
geolocation grid as ground control points, and the calibration constant sits in
the product's own Main Processing Parameters record. Each product carries two
detected amplitude images, MDS1 and MDS2, whose polarisations are named in the
Specific Product Header (HH and HV for the products over Lake Chilwa).

Calibration follows the ASAR Product Handbook (ESA 2007, section 2.6) and
Rosich and Meadows (2004): for a precision image, already corrected for the
antenna pattern and range spreading loss,

    sigma0 = DN^2 * sin(alpha) / K

with alpha the incidence angle interpolated from the tie-point grid and K the
external calibration constant of that image. No terrain correction is applied.
The basin floor and the wetland are flat, which is where the join is read;
backscatter on the surrounding hills carries geometric and radiometric terrain
effects and is not used.

Geocoding is a thin-plate-spline warp on the tie-point ground control points,
resampled by averaging from the 12.5 m product spacing onto the 30 m grid of
03.outputs/TIF/harmonic_2012/mndwi_harmonic_2012-02-22.tif, so the radar and
the optical model share one pixel grid. Averaging in linear power is a
multilook, which is the first speckle reduction.

Usage
    python3 asar_process.py FILE.N1 [FILE.N1 ...]

Writes 03.outputs/TIF/asar/<product>_sigma0.tif, float32, three bands,
sigma0 HH (linear), sigma0 HV (linear), incidence angle (degrees), nodata NaN,
and appends one row per product to 03.outputs/CSV/asar_products.csv.
"""

import csv
import os
import struct
import subprocess
import sys
import tempfile

import numpy as np
import rasterio
from rasterio.transform import Affine
from scipy.interpolate import RegularGridInterpolator

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GRID = os.path.join(ROOT, "03.outputs", "TIF", "harmonic_2012",
                    "mndwi_harmonic_2012-02-22.tif")
OUT_TIF = os.path.join(ROOT, "03.outputs", "TIF", "asar")
OUT_CSV = os.path.join(ROOT, "03.outputs", "CSV", "asar_products.csv")
os.makedirs(OUT_TIF, exist_ok=True)

MPH_SIZE = 1247


def read_headers(path):
    """The ASCII main and specific product headers and the data set descriptors."""
    with open(path, "rb") as fh:
        mph = fh.read(MPH_SIZE).decode("latin1")
        fields = dict(line.split("=", 1) for line in mph.splitlines() if "=" in line)
        sph_size = int(fields["SPH_SIZE"].rstrip("<bytes>"))
        sph = fh.read(sph_size).decode("latin1")
    sfields = {}
    for line in sph.splitlines():
        if "=" in line:
            k, v = line.split("=", 1)
            sfields[k] = v.strip().strip('"').strip()
    dsds = []
    for block in sph.split("DS_NAME=")[1:]:
        d = {"DS_NAME": block.split('"')[1]}
        for line in block.splitlines()[1:]:
            if "=" in line:
                k, v = line.split("=", 1)
                d[k] = v.strip().strip('"').rstrip("<bytes>").strip()
        dsds.append(d)
    return fields, sfields, dsds


def find_dsd(dsds, name):
    for d in dsds:
        if d["DS_NAME"].strip() == name:
            return d
    raise KeyError(name)


def geolocation_grid(path, dsds):
    """Tie-point lines, sample numbers, incidence angles, latitudes, longitudes.

    Record layout from the ASAR Products Specification (PO-RS-MDA-GS-2009),
    Geolocation Grid ADSR, big-endian: MJD time (12), attach flag (1), first
    line number (4), number of lines (4), sub-satellite track heading (4), then
    eleven samples each of sample number (unsigned integer), slant range time,
    incidence angle, latitude and longitude (five arrays of 11 x 4 bytes), and
    the same block again for the last line of the granule. The last-line block
    was found by search at byte 254 of the 521-byte record, where its latitudes
    reproduce the SPH_LAST_NEAR_LAT header value; the line number of that block
    is taken as first line plus number of lines minus one."""
    d = find_dsd(dsds, "GEOLOCATION GRID ADS")
    off, n, size = int(d["DS_OFFSET"]), int(d["NUM_DSR"]), int(d["DSR_SIZE"])
    lines, samps, angs, lats, lons = [], [], [], [], []
    with open(path, "rb") as fh:
        fh.seek(off)
        for i in range(n):
            rec = fh.read(size)
            first = struct.unpack(">I", rec[13:17])[0]
            nlines = struct.unpack(">I", rec[17:21])[0]
            # first-line block
            s = np.array(struct.unpack(">11I", rec[25:69]), dtype=float)
            a = np.array(struct.unpack(">11f", rec[113:157]))
            la = np.array(struct.unpack(">11i", rec[157:201])) * 1e-6
            lo = np.array(struct.unpack(">11i", rec[201:245])) * 1e-6
            lines.append(first); samps.append(s); angs.append(a); lats.append(la); lons.append(lo)
            if i == n - 1:
                # last-line block, same layout starting at 267
                b = 254
                last = first + nlines - 1
                s2 = np.array(struct.unpack(">11I", rec[b + 25:b + 69]), dtype=float)
                a2 = np.array(struct.unpack(">11f", rec[b + 113:b + 157]))
                la2 = np.array(struct.unpack(">11i", rec[b + 157:b + 201])) * 1e-6
                lo2 = np.array(struct.unpack(">11i", rec[b + 201:b + 245])) * 1e-6
                lines.append(last); samps.append(s2); angs.append(a2); lats.append(la2); lons.append(lo2)
    return (np.array(lines, dtype=float), np.array(samps), np.array(angs),
            np.array(lats), np.array(lons))


def calibration_constants(path, dsds, sfields, ds):
    """External calibration constant K for MDS1 and MDS2.

    GDAL exposes the Main Processing Parameters record field by field in its
    RECORDS metadata domain; the two constants are read from there."""
    md = ds.tags(ns="RECORDS")
    keys = sorted(k for k in md if k.endswith("EXT_CAL_FACT"))
    return md, keys


def process(path):
    name = os.path.basename(path).replace(".N1", "")
    out = os.path.join(OUT_TIF, f"{name}_sigma0.tif")
    fields, sfields, dsds = read_headers(path)
    pol1 = sfields.get("MDS1_TX_RX_POLAR", "?")
    pol2 = sfields.get("MDS2_TX_RX_POLAR", "?")
    swath = sfields.get("SWATH", "?")
    sensing = fields.get("SENSING_START", "").strip('"').strip()
    print(f"{name}: swath {swath}, MDS1 {pol1}, MDS2 {pol2}, {sensing}", flush=True)

    with rasterio.open(path) as ds:
        md, keys = calibration_constants(path, dsds, sfields, ds)
        ncol, nrow = ds.width, ds.height
        lines, samps, angs, lats, lons = geolocation_grid(path, dsds)
        gcps, _ = ds.gcps
        # tie points must agree with the GCPs GDAL reads from the same record
        g0 = gcps[0]
        print(f"  {ncol} x {nrow} px, {len(lines)} tie lines x 11, "
              f"first tie point lat {lats[0,0]:.5f} lon {lons[0,0]:.5f}; "
              f"GDAL GCP 0 lat {g0.y:.5f} lon {g0.x:.5f} at col {g0.col} row {g0.row}",
              flush=True)
        print(f"  incidence angle range {angs.min():.2f} to {angs.max():.2f} deg",
              flush=True)
        for k in keys:
            print(f"  {k} = {md[k]}", flush=True)
        if len(keys) < 2:
            raise SystemExit("calibration constants not found in RECORDS metadata")
        K = [float(md[k]) for k in keys][:2]

        # incidence angle interpolator over (line, sample), sample numbers in
        # the product are one-based
        assert lines[-1] == nrow, (lines[-1], nrow)
        assert np.all(np.diff(lines) > 0) and angs.min() > 10, "tie-point grid misread"
        interp = RegularGridInterpolator((lines - 1, samps[0] - 1), angs,
                                         bounds_error=False, fill_value=None)

        tmp = os.path.join(tempfile.gettempdir(), f"{name}_tmp.tif")
        profile = {"driver": "GTiff", "width": ncol, "height": nrow, "count": 3,
                   "dtype": "float32", "nodata": np.nan, "tiled": True,
                   "compress": "DEFLATE", "BIGTIFF": "YES"}
        with rasterio.open(tmp, "w", **profile) as dst:
            dst.gcps = (gcps, rasterio.crs.CRS.from_epsg(4326))
            cols = np.arange(ncol, dtype=float)
            for r0 in range(0, nrow, 512):
                r1 = min(r0 + 512, nrow)
                rows = np.arange(r0, r1, dtype=float)
                rr, cc = np.meshgrid(rows, cols, indexing="ij")
                ang = interp(np.stack([rr.ravel(), cc.ravel()], axis=1)).reshape(r1 - r0, ncol)
                sin_a = np.sin(np.deg2rad(ang)).astype(np.float32)
                win = rasterio.windows.Window(0, r0, ncol, r1 - r0)
                for b in (1, 2):
                    dn = ds.read(b, window=win).astype(np.float32)
                    s0 = dn * dn * sin_a / np.float32(K[b - 1])
                    s0[dn == 0] = np.nan
                    dst.write(s0, b, window=win)
                dst.write(ang.astype(np.float32), 3, window=win)

    with rasterio.open(GRID) as g:
        b = g.bounds
        res = g.res[0]
    subprocess.run(["gdalwarp", "-overwrite", "-tps", "-r", "average",
                    "-t_srs", "EPSG:4326", "-te", str(b.left), str(b.bottom),
                    str(b.right), str(b.top), "-tr", str(res), str(res),
                    "-srcnodata", "nan", "-dstnodata", "nan", "-multi", "-wo",
                    "NUM_THREADS=4", "-co", "COMPRESS=DEFLATE", "-co", "TILED=YES",
                    tmp, out], check=True, capture_output=True)
    os.remove(tmp)
    with rasterio.open(out, "r+") as o:
        o.descriptions = (f"sigma0_{pol1.replace('/', '')}", f"sigma0_{pol2.replace('/', '')}",
                          "incidence_deg")
        o.update_tags(SWATH=swath, MDS1_POL=pol1, MDS2_POL=pol2, SENSING_START=sensing,
                      K_MDS1=K[0], K_MDS2=K[1], SOURCE=os.path.basename(path))
        arr = o.read(1)
        valid = np.isfinite(arr).sum()
    print(f"  wrote {out}: {valid * (res * 111.2) ** 2 * np.cos(np.deg2rad(15.3)):,.0f} km2 "
          f"valid inside the harmonic grid", flush=True)

    new = not os.path.exists(OUT_CSV)
    with open(OUT_CSV, "a", newline="") as fh:
        w = csv.writer(fh)
        if new:
            w.writerow(["product", "date", "sensing_start", "swath", "mds1_pol",
                        "mds2_pol", "K_mds1", "K_mds2", "valid_km2_in_grid"])
        w.writerow([name, name[14:18] + "-" + name[18:20] + "-" + name[20:22], sensing,
                    swath, pol1, pol2, K[0], K[1],
                    round(valid * (res * 111.2) ** 2 * np.cos(np.deg2rad(15.3)), 1)])
    return out


if __name__ == "__main__":
    for p in sys.argv[1:]:
        process(p)
