#!/usr/bin/env python3
"""Water standing within vegetation, and open water, from sub-pixel spectral
mixture analysis carried through the same harmonic model as the index series.

Two quantities, one framework. Every clear Landsat observation over the basin
is unmixed against the four committed endmembers, open water, emergent or
flooded vegetation, dry vegetation and bare soil, exactly as
sma_endmember_modelling.R does, with the fractions constrained to be
non-negative and to sum to one. The resulting fraction is then treated as the
signal a harmonic model is fitted to, in place of the water index, so the
seasonal shape, the window, the model order and the observation floor are the
ones the open-water series already settled for that year. Reading the fitted
fraction at twenty-four dates through the target year and summing fraction
times pixel area gives an area in square kilometres without ever thresholding,
which is the sub-pixel measurement Halabisky et al. (2016) make and the reason
that paper is this study's method twin for the inundation core.

Doing it this way means the optical estimate of the vegetated class is built on
the same window and the same observations as the open-water column beside it,
so the two are comparable, and the difference from the radar estimate in the
years that carry both is a measurement of the optical method's error rather
than a difference of design. Hypothesis H2 says that difference should be
large, because optical reflectance sees only the top of the reed stand while
the radar sees the water beneath it.

Endmembers are the committed image-derived spectra in
03.outputs/CSV/sma_endmember_spectra.csv, the Landsat 5 set for the Thematic
Mapper and Enhanced Thematic Mapper Plus and the Landsat 8 set for the
Operational Land Imager, because band positions differ between them.

Writes 03.outputs/CSV/sma_harmonic_fractions.csv
"""

import argparse
import csv
import os
import time

import ee

import harmonic_open_water_series as H

ee.Initialize(project="murphys-deforisk")

ROOT = H.ROOT
AOI = H.AOI
AREA = H.AREA_KM2
OUT = os.path.join(ROOT, "03.outputs", "CSV", "sma_harmonic_fractions.csv")
SPECTRA = os.path.join(ROOT, "03.outputs", "CSV", "sma_endmember_spectra.csv")

CLASSES = ["water", "flooded_veg", "dry_veg", "bare_soil"]
SIX = ["blue", "green", "red", "nir", "swir1", "swir2"]
BANDS6 = {"LT04": ["SR_B1", "SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B7"],
          "LT05": ["SR_B1", "SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B7"],
          "LE07": ["SR_B1", "SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B7"],
          "LC08": ["SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B6", "SR_B7"],
          "LC09": ["SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B6", "SR_B7"]}
ERA = {"LT04": "L5", "LT05": "L5", "LE07": "L5", "LC08": "L8", "LC09": "L8"}


def endmembers():
    em = {}
    for r in csv.DictReader(open(SPECTRA)):
        era = r["era"].strip('"')
        em.setdefault(era, {})[r["class"].strip('"')] = [float(r[b]) for b in SIX]
    return em


EM = endmembers()


def observations(start, end, sensors, harmonics, target):
    """Clear scenes, unmixed, with the design row attached.

    `target` names which endmember fraction becomes the signal."""
    out = None
    for s in sensors:
        col = (ee.ImageCollection(f"LANDSAT/{s}/C02/T1_L2")
               .filterBounds(AOI).filterDate(start, end).map(H.clear))
        spec = [EM[ERA[s]][c] for c in CLASSES]
        bands = BANDS6[s]

        def prep(img, bands=bands, spec=spec):
            d = ee.Date(img.get("system:time_start"))
            ref = (img.select(bands).multiply(0.0000275).add(-0.2).rename(SIX))
            frac = ref.unmix(spec, True, True).rename(CLASSES)
            sig = frac.select(target).rename("y")
            return (H.design(d, harmonics).addBands(sig).updateMask(sig.mask())
                    .set("system:time_start", img.get("system:time_start")))

        col = col.map(prep)
        out = col if out is None else out.merge(col)
    return out


def window_for(year):
    path = os.path.join(ROOT, "03.outputs", "CSV", "harmonic_open_water_annual.csv")
    for r in csv.DictReader(open(path)):
        if int(r["year"]) == year:
            harm = (H.FULL_HARMONICS if r["model_order"] == "annual+semiannual"
                    else H.SIMPLE_HARMONICS)
            return int(r["half_window"]), harm
    raise SystemExit(f"{year} is not in the open-water series yet; run that first")


def year_row(year, retries=3):
    hw, harmonics = window_for(year)
    y0, y1 = max(year - hw, 1984), min(year + hw + 1, 2026)
    sens = H.sensors_for(y0, y1)
    px = H.predictors(harmonics)
    floor = H.OBS_PER_PARAM * len(px)

    out = {"year": year, "window_start": y0, "window_end": y1 - 1,
           "model_order": "annual+semiannual" if len(harmonics) == 2 else "annual",
           "sensors": "+".join(sens)}

    for target, label in (("flooded_veg", "flooded"), ("water", "open")):
        col = observations(f"{y0}-01-01", f"{y1}-01-01", sens, harmonics, target)
        n = col.select("y").count().rename("n")
        coeffs = (col.select(px + ["y"])
                  .reduce(ee.Reducer.linearRegression(numX=len(px), numY=1))
                  .select("coefficients").arrayProject([0]).arrayFlatten([px])
                  .updateMask(n.gte(floor)))
        km2 = []
        for i in range(H.N_DATES):
            d = ee.Date(f"{year}-01-01").advance((i + 0.5) / H.N_DATES * 365.25, "day")
            f = coeffs.multiply(H.design(d, harmonics)).reduce(ee.Reducer.sum())
            km2.append(f.clamp(0, 1).multiply(AREA).rename(f"d{i:02d}"))
        area = ee.Image.cat(km2).reduceRegion(
            ee.Reducer.sum(), AOI, H.SCALE, maxPixels=1e10, bestEffort=True)

        def resid(img):
            pred = coeffs.multiply(img.select(px)).reduce(ee.Reducer.sum())
            return img.select("y").subtract(pred).pow(2).rename("sq")

        stats = ee.Dictionary({
            "area": area, "n_scenes": col.size(),
            "rmse": col.map(resid).mean().sqrt().reduceRegion(
                ee.Reducer.median(), AOI, 300, maxPixels=1e10, bestEffort=True).get("sq"),
            "cov": coeffs.select("const").mask().multiply(AREA).rename("km2").reduceRegion(
                ee.Reducer.sum(), AOI, H.SCALE, maxPixels=1e10, bestEffort=True).get("km2"),
        })
        last = None
        for a in range(retries):
            try:
                g = stats.getInfo(); break
            except Exception as e:
                last = e; time.sleep(20 * (a + 1))
        else:
            raise RuntimeError(f"{year} {target}: {last}")

        vals = [g["area"][f"d{i:02d}"] for i in range(H.N_DATES)]
        mean = sum(vals) / len(vals)
        # the fraction's own fit residual, in fraction units, times the basin
        # gives the area equivalent of the model error
        sd = (g["rmse"] or 0.0) * H.BASIN_KM2_NOMINAL / (len(vals) ** 0.5)
        out[f"{label}_km2"] = round(mean, 1)
        out[f"{label}_sd_km2"] = round(sd, 1)
        out[f"{label}_max_km2"] = round(max(vals), 1)
        out[f"{label}_min_km2"] = round(min(vals), 1)
        out[f"{label}_fit_residual"] = round(g["rmse"] or 0.0, 4)
        out[f"{label}_n_scenes"] = g["n_scenes"]
        out[f"{label}_coverage"] = round((g["cov"] or 0.0) / H.BASIN_KM2_NOMINAL, 4)
    out["method"] = (f"four-endmember spectral mixture analysis per scene, non-negative and "
                     f"summing to one, then a harmonic fit on the fraction over "
                     f"{y0} to {y1 - 1}, read at {H.N_DATES} dates and integrated sub-pixel")
    return out


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
    for year in range(args.start, args.end + 1):
        if year in done:
            continue
        r = year_row(year)
        done[year] = r
        print(f"  {year}  vegetated {r['flooded_km2']:7.0f} +/- {r['flooded_sd_km2']:5.0f} km2, "
              f"open {r['open_km2']:7.0f} +/- {r['open_sd_km2']:5.0f} km2, "
              f"window {r['window_start']}-{r['window_end']}, "
              f"residuals {r['flooded_fit_residual']:.3f} and {r['open_fit_residual']:.3f}",
              flush=True)
        with open(OUT, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(done[max(done)].keys()))
            w.writeheader()
            for y in sorted(done):
                w.writerow(done[y])
    print(f"wrote {OUT}, {len(done)} years")


if __name__ == "__main__":
    main()
