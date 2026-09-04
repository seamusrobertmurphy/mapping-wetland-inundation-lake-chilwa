#!/usr/bin/env python3
"""Does HH minus HV separate flooded reeds from open water and from dry land
on one Envisat ASAR date, read against the optical harmonic fit for the same
date?

The physics being tested. A horizontally polarised wave passes between the
vertical reed stems, reflects off the water surface, bounces into the stem
base and returns doubled, so HH is bright over flooded vegetation. HV is
produced by volume scattering within the stand and does not share that gain.
Open water reflects the wave away from the sensor in both channels. Dry
vegetation scatters moderately in both. The difference HH minus HV in decibels
should therefore be highest over flooded reeds, and the incidence-angle term of
the calibration cancels in it exactly because both channels share one geometry.

Reference strata come from the optical model on the same date and from the
lake outline, not from field points, which carry no coordinates. The optical
fit places open water; it cannot place flooded reeds, so the test is whether
the radar difference is bimodal within the optically non-water belt around the
open water and whether the bright mode lies against the open water, where the
reed belt is, rather than scattered across the basin.

Usage
    python3 asar_hhhv_test.py 2012-02-22

Reads 03.outputs/TIF/asar/*<date>*_sigma0.tif, the harmonic coefficients, the
date terms and the lake and basin outlines. Writes
    03.outputs/CSV/asar_hhhv_test_<date>.csv        stratum statistics
    03.outputs/CSV/asar_hhhv_mixture_<date>.csv     the two-component mixture
    03.outputs/TIF/asar/asar_join_<date>.tif        HH dB, HV dB, HH-HV dB,
                                                    fitted MNDWI, stratum, class
    03.outputs/PNG/asar_hhhv_test_<date>.png
"""

import csv
import glob
import os
import sys

import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from rasterio import features
from scipy import ndimage

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TIF = os.path.join(ROOT, "03.outputs", "TIF")
CSV = os.path.join(ROOT, "03.outputs", "CSV")
PNG = os.path.join(ROOT, "03.outputs", "PNG")
COEF = os.path.join(TIF, "harmonic_2012", "harmonic_coefficients.tif")
TERMS = os.path.join(CSV, "harmonic_date_terms.csv")
LAKE = os.path.join(ROOT, "02.inputs", "SHP", "lakes_site.shp")
BASIN = os.path.join(ROOT, "03.outputs", "JSON", "chilwa_basin.geojson")

MEDIAN_WINDOW = 5      # pixels at 30 m, applied to the decibel images
WATER = 0.10           # fitted MNDWI above which a pixel is confidently open water
DRY = -0.10            # fitted MNDWI below which a pixel is confidently not open water
BELT_KM = 5.0          # how far outside the lake outline the reed belt may extend
FAR_KM = 10.0          # beyond this from the lake outline is the dry-land reference
MIN_COVERAGE = 0.25    # radar must cover this share of the lake vicinity for a class map


def fitted_mndwi(date):
    """Read the optical model on one date from the exported coefficients, or
    from the raster harmonic_fit_2012.py wrote for that date if it exists."""
    direct = os.path.join(TIF, "harmonic_2012", f"mndwi_harmonic_{date}.tif")
    if os.path.exists(direct):
        with rasterio.open(direct) as d:
            m = d.read(1).astype(np.float32)
            m[d.read(1) == d.nodata] = np.nan
            return m / 10000.0, d.profile
    with open(TERMS) as fh:
        rows = {r["date"]: r for r in csv.DictReader(fh)}
    r = rows[date]
    with rasterio.open(COEF) as c:
        coef = c.read().astype(np.float32) / 10000.0
        nodata = c.read(1) == c.nodata
        profile = c.profile
    t = float(r["t_years"])
    terms = np.array([1.0, t, float(r["cos1"]), float(r["sin1"]),
                      float(r["cos2"]), float(r["sin2"])], dtype=np.float32)
    m = np.tensordot(terms, coef, axes=1)
    m[nodata] = np.nan
    return m, profile


def db(x):
    with np.errstate(divide="ignore", invalid="ignore"):
        return 10.0 * np.log10(x)


def smooth(a, size):
    """Median filter that ignores NaN by filling with the median first."""
    fill = np.nanmedian(a)
    filled = np.where(np.isfinite(a), a, fill)
    out = ndimage.median_filter(filled, size=size)
    out[~np.isfinite(a)] = np.nan
    return out


def pixel_area_km2(profile):
    """Area of each row's pixels on a geographic grid."""
    tr = profile["transform"]
    lat = tr.f + tr.e * (np.arange(profile["height"]) + 0.5)
    dx = tr.a * 111320.0 * np.cos(np.deg2rad(lat))
    dy = -tr.e * 111320.0
    return (dx * dy / 1e6).astype(np.float32)[:, None]


def gmm(x, k, iters=200, seed=1):
    """One-dimensional Gaussian mixture by expectation maximisation.
    Returns means, standard deviations, weights and the Bayesian information
    criterion, so a one-component and a two-component fit can be compared."""
    from scipy.stats import norm
    x = np.asarray(x, dtype=np.float64).ravel()
    q = np.percentile(x, np.linspace(20, 80, k))
    mu, sd, wt = q.copy(), np.full(k, x.std()), np.full(k, 1.0 / k)
    for _ in range(iters):
        dens = np.stack([wt[j] * norm.pdf(x, mu[j], sd[j]) for j in range(k)])
        tot = dens.sum(axis=0) + 1e-300
        resp = dens / tot
        nk = resp.sum(axis=1)
        mu = (resp * x).sum(axis=1) / nk
        sd = np.sqrt((resp * (x - mu[:, None]) ** 2).sum(axis=1) / nk) + 1e-6
        wt = nk / len(x)
    ll = np.log(np.stack([wt[j] * norm.pdf(x, mu[j], sd[j]) for j in range(k)]).sum(axis=0)
                + 1e-300).sum()
    bic = (3 * k - 1) * np.log(len(x)) - 2 * ll
    return mu, sd, wt, bic


def m_stat(a, b):
    """Separability index M of Kaufman and Remer (1994), as in the study's
    own spectral screen: distance between means over the sum of the standard
    deviations, with 1 the usual threshold for good separation."""
    return abs(np.nanmean(a) - np.nanmean(b)) / (np.nanstd(a) + np.nanstd(b))


def main(date, figure=True):
    tag = date.replace("-", "")
    files = sorted(glob.glob(os.path.join(TIF, "asar", f"ASA_APP_1PNESA{tag}_*_sigma0.tif")))
    if not files:
        raise SystemExit(f"no processed ASAR frames for {date}")
    print(f"{date}: {len(files)} frame(s)")
    hh = hv = inc = None
    for f in files:
        with rasterio.open(f) as ds:
            # the product order of the two polarisations varies, so bands are
            # taken by the label asar_process.py wrote, never by position
            desc = list(ds.descriptions)
            a = np.stack([ds.read(desc.index("sigma0_HH") + 1),
                          ds.read(desc.index("sigma0_HV") + 1),
                          ds.read(desc.index("incidence_deg") + 1)]).astype(np.float32)
        if hh is None:
            hh, hv, inc = a[0], a[1], a[2]
        else:
            # frames along one track overlap by a few lines; take the mean there
            hh = np.where(np.isfinite(hh) & np.isfinite(a[0]), (hh + a[0]) / 2,
                          np.where(np.isfinite(hh), hh, a[0]))
            hv = np.where(np.isfinite(hv) & np.isfinite(a[1]), (hv + a[1]) / 2,
                          np.where(np.isfinite(hv), hv, a[1]))
            inc = np.where(np.isfinite(inc), inc, a[2])

    mndwi, profile = fitted_mndwi(date)
    assert mndwi.shape == hh.shape, (mndwi.shape, hh.shape)

    hh_db = smooth(db(hh), MEDIAN_WINDOW)
    hv_db = smooth(db(hv), MEDIAN_WINDOW)
    diff = hh_db - hv_db

    lake = gpd.read_file(LAKE).to_crs(4326)
    basin = gpd.read_file(BASIN).to_crs(4326)
    lake_m = lake.to_crs(32736)
    shp = (profile["height"], profile["width"])
    tr = profile["transform"]
    in_basin = features.rasterize(basin.geometry, out_shape=shp, transform=tr, fill=0,
                                  default_value=1).astype(bool)
    belt = features.rasterize(lake_m.buffer(BELT_KM * 1000).to_crs(4326).geometry,
                              out_shape=shp, transform=tr, fill=0, default_value=1).astype(bool)
    near = features.rasterize(lake_m.buffer(FAR_KM * 1000).to_crs(4326).geometry,
                              out_shape=shp, transform=tr, fill=0, default_value=1).astype(bool)

    valid = np.isfinite(diff) & np.isfinite(mndwi) & in_basin
    stratum = np.zeros(shp, dtype=np.uint8)
    stratum[valid & (mndwi > WATER)] = 1                      # open water, optical
    stratum[valid & belt & (mndwi < DRY)] = 2                 # reed belt candidates
    stratum[valid & ~near & (mndwi < DRY)] = 3                # dry land, far from lake
    stratum[valid & (mndwi >= DRY) & (mndwi <= WATER)] = 4    # optical transition
    names = {1: "open water (fitted MNDWI > 0.10)",
             2: f"belt within {BELT_KM:.0f} km of the lake outline, fitted MNDWI < -0.10",
             3: f"dry land beyond {FAR_KM:.0f} km of the lake outline, fitted MNDWI < -0.10",
             4: "optical transition, -0.10 to 0.10"}
    area = pixel_area_km2(profile)
    area_img = np.broadcast_to(area, shp)

    # Coverage of the lake vicinity by this date's frames. The four tracks
    # place the swath differently and two of them miss the lake, so an area
    # is only meaningful where the radar looked. Open water from the optical
    # model alone is kept for every date as the context the radar joins to.
    lake_px = features.rasterize(lake.geometry, out_shape=shp, transform=tr, fill=0,
                                 default_value=1).astype(bool)
    dist_km = ndimage.distance_transform_edt(~lake_px) * 30.0 / 1000.0
    vicinity = in_basin & (lake_px | (dist_km <= 10))
    covered = float(area_img[vicinity & np.isfinite(diff)].sum() / area_img[vicinity].sum())
    optical_all = float(area_img[in_basin & np.isfinite(mndwi) & (mndwi > WATER)].sum())
    optical_all0 = float(area_img[in_basin & np.isfinite(mndwi) & (mndwi > 0)].sum())
    print(f"  lake vicinity covered by radar: {100 * covered:.1f} per cent; optical open water "
          f"on this date {optical_all:,.0f} km2 at MNDWI > 0.10, {optical_all0:,.0f} km2 at > 0")
    summary = {"date": date, "frames": len(files), "vicinity_covered": round(covered, 4),
               "optical_open_water_km2_gt_0.10": round(optical_all, 1),
               "optical_open_water_km2_gt_0": round(optical_all0, 1)}
    if covered < MIN_COVERAGE or (stratum == 2).sum() < 5000:
        print(f"  radar covers under {100 * MIN_COVERAGE:.0f} per cent of the lake vicinity; "
              "no class map for this date")
        return summary

    rows = []
    for s, name in names.items():
        m = stratum == s
        if m.sum() == 0:
            continue
        rows.append({"stratum": name, "n_px": int(m.sum()),
                     "area_km2": round(float(area_img[m].sum()), 1),
                     "HH_dB_median": round(float(np.median(hh_db[m])), 2),
                     "HH_dB_q25": round(float(np.percentile(hh_db[m], 25)), 2),
                     "HH_dB_q75": round(float(np.percentile(hh_db[m], 75)), 2),
                     "HV_dB_median": round(float(np.median(hv_db[m])), 2),
                     "HV_dB_q25": round(float(np.percentile(hv_db[m], 25)), 2),
                     "HV_dB_q75": round(float(np.percentile(hv_db[m], 75)), 2),
                     "HHminusHV_dB_median": round(float(np.median(diff[m])), 2),
                     "HHminusHV_dB_q25": round(float(np.percentile(diff[m], 25)), 2),
                     "HHminusHV_dB_q75": round(float(np.percentile(diff[m], 75)), 2),
                     "incidence_deg_median": round(float(np.nanmedian(inc[m])), 1)})
    with open(os.path.join(CSV, f"asar_hhhv_test_{date}.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print(f"  {r['stratum']}: {r['area_km2']:,.0f} km2, HH {r['HH_dB_median']} dB, "
              f"HV {r['HV_dB_median']} dB, HH-HV {r['HHminusHV_dB_median']} dB")

    w_, b_, d_ = (stratum == 1), (stratum == 2), (stratum == 3)
    print("  separability M, water against dry land: "
          f"HH {m_stat(hh_db[w_], hh_db[d_]):.2f}, HV {m_stat(hv_db[w_], hv_db[d_]):.2f}, "
          f"HH-HV {m_stat(diff[w_], diff[d_]):.2f}")

    # Is the belt bimodal in HH minus HV, and where does its bright mode sit?
    x = diff[b_].reshape(-1, 1).astype(np.float64)
    rng = np.random.default_rng(1)
    xs = x if len(x) <= 200000 else x[rng.choice(len(x), 200000, replace=False)]
    _, _, _, bic1 = gmm(xs, 1)
    mu, sd, wt, bic2 = gmm(xs, 2)
    order = np.argsort(mu)
    mu, sd, wt = mu[order], sd[order], wt[order]
    ashman_d = abs(mu[1] - mu[0]) / np.sqrt((sd[0] ** 2 + sd[1] ** 2) / 2)
    # the threshold is where the two weighted densities cross
    grid = np.linspace(mu[0], mu[1], 2001)
    from scipy.stats import norm
    d0 = wt[0] * norm.pdf(grid, mu[0], sd[0])
    d1 = wt[1] * norm.pdf(grid, mu[1], sd[1])
    cross = grid[np.argmin(np.abs(d0 - d1))]
    print(f"  belt mixture: modes {mu[0]:.2f} and {mu[1]:.2f} dB, sd {sd[0]:.2f} and {sd[1]:.2f}, "
          f"weights {wt[0]:.2f} and {wt[1]:.2f}, Ashman D {ashman_d:.2f}, "
          f"BIC one {bic1:,.0f} two {bic2:,.0f}, crossing {cross:.2f} dB")
    with open(os.path.join(CSV, f"asar_hhhv_mixture_{date}.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["date", "component", "mean_dB", "sd_dB", "weight", "ashman_D",
                    "bic_one", "bic_two", "crossing_dB", "n_fit"])
        for i in range(2):
            w.writerow([date, i + 1, round(mu[i], 3), round(sd[i], 3), round(wt[i], 3),
                        round(ashman_d, 3), round(bic1, 1), round(bic2, 1),
                        round(cross, 3), len(xs)])

    # Class map. Flooded vegetation is brighter than dry land in HH, because
    # of the double bounce, and brighter in HH minus HV, because HV does not
    # share that gain. Both thresholds are the 95th percentile of dry land far
    # from the lake, so the rule admits at most 5 per cent of dry land on
    # either channel, and the joint rate on dry land is reported as the noise
    # floor. The rule is applied to every optically non-water pixel in the
    # basin, not only to the belt, so that the fall-off with distance from the
    # lake is itself a test of whether the class is the wetland fringe.
    t_diff = float(np.percentile(diff[d_], 95))
    t_hh = float(np.percentile(hh_db[d_], 95))
    nonwater = valid & (mndwi <= WATER)
    bright = nonwater & (diff > t_diff) & (hh_db > t_hh)
    floor = float(area_img[bright & d_].sum() / area_img[d_].sum())
    print(f"  thresholds, dry-land 95th percentiles: HH-HV > {t_diff:.2f} dB and HH > {t_hh:.2f} dB; "
          f"joint rate on dry land {100 * floor:.1f} per cent")
    bands = [(0, 1), (1, 2), (2, 3), (3, 5), (5, 7), (7, 10), (10, 15), (15, 25), (25, 80)]
    dist_rows = []
    m0 = nonwater & lake_px
    dist_rows.append(("inside the lake outline", float(area_img[m0].sum()),
                      float(area_img[m0 & bright].sum())))
    for lo, hi in bands:
        m = nonwater & ~lake_px & (dist_km >= lo) & (dist_km < hi)
        dist_rows.append((f"{lo} to {hi} km outside", float(area_img[m].sum()),
                          float(area_img[m & bright].sum())))
    for lab, a_all, a_br in dist_rows:
        print(f"    {lab:26s} {a_all:8,.0f} km2 non-water, {a_br:7,.0f} km2 flooded "
              f"({100 * a_br / max(a_all, 1e-9):4.1f} per cent)")
    cls = np.zeros(shp, dtype=np.uint8)
    cls[valid & (mndwi > WATER)] = 1
    cls[bright] = 2
    a_water = float(area_img[cls == 1].sum())
    a_flood = float(area_img[cls == 2].sum())
    a_flood_near = float(area_img[(cls == 2) & (dist_km <= 10)].sum())
    print(f"  open water {a_water:,.0f} km2; flooded vegetation {a_flood:,.0f} km2 in the basin, "
          f"{a_flood_near:,.0f} km2 of it within 10 km of the lake outline; "
          f"combined within 10 km {a_water + a_flood_near:,.0f} km2")
    with open(os.path.join(CSV, f"asar_hhhv_test_{date}.csv"), "a", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([])
        w.writerow(["quantity", "km2"])
        w.writerow(["open water, optical fit, MNDWI > 0.10", round(a_water, 1)])
        w.writerow([f"flooded vegetation, HH-HV > {t_diff:.2f} dB and HH > {t_hh:.2f} dB, whole basin",
                    round(a_flood, 1)])
        w.writerow(["flooded vegetation within 10 km of the lake outline", round(a_flood_near, 1)])
        w.writerow(["combined, open water plus flooded within 10 km", round(a_water + a_flood_near, 1)])
        w.writerow(["joint rule rate on dry land beyond 10 km (noise floor)", round(floor, 4)])
        w.writerow([])
        w.writerow(["distance from lake outline", "non-water km2", "flooded km2", "fraction"])
        for lab, a_all, a_br in dist_rows:
            w.writerow([lab, round(a_all, 1), round(a_br, 1), round(a_br / max(a_all, 1e-9), 4)])
    cross = t_diff  # drawn on the histogram panel

    belt_nonwater = nonwater & (dist_km <= 10) & np.isfinite(diff)
    summary.update({"open_water_km2_mndwi_gt_0.10": round(a_water, 1),
               "open_water_km2_mndwi_gt_0": round(float(area_img[valid & (mndwi > 0)].sum()), 1),
               "flooded_km2_within_10km": round(a_flood_near, 1),
               "flooded_share_of_nonwater_within_10km": round(
                   float(area_img[bright & (dist_km <= 10)].sum() / max(area_img[belt_nonwater].sum(), 1e-9)), 4),
               "combined_km2_within_10km": round(a_water + a_flood_near, 1),
               "t_hhhv_dB": round(t_diff, 2), "t_hh_dB": round(t_hh, 2),
               "dry_land_false_positive": round(floor, 4),
               "HH_water_median_dB": round(float(np.median(hh_db[w_])), 2),
               "HH_dry_median_dB": round(float(np.median(hh_db[d_])), 2),
               "HV_water_median_dB": round(float(np.median(hv_db[w_])), 2),
               "HV_dry_median_dB": round(float(np.median(hv_db[d_])), 2),
               "incidence_lake_median_deg": round(float(np.nanmedian(inc[w_])), 1)})

    out = os.path.join(TIF, "asar", f"asar_join_{date}.tif")
    prof = dict(profile, count=6, dtype="float32", nodata=np.nan, compress="DEFLATE",
                tiled=True)
    with rasterio.open(out, "w", **prof) as dst:
        for i, arr in enumerate([hh_db, hv_db, diff, mndwi, stratum.astype(np.float32),
                                 cls.astype(np.float32)], 1):
            dst.write(arr.astype(np.float32), i)
        dst.descriptions = ("HH_dB", "HV_dB", "HHminusHV_dB", "fitted_MNDWI", "stratum",
                            "class_1water_2flooded")
    print(f"  wrote {out}")
    if not figure:
        return summary

    # Figure: the difference image over the lake, the class map, the histograms.
    rows_, cols_ = np.where(belt)
    r0, r1 = max(rows_.min() - 60, 0), min(rows_.max() + 60, shp[0])
    c0, c1 = max(cols_.min() - 60, 0), min(cols_.max() + 60, shp[1])
    ext = [tr.c + tr.a * c0, tr.c + tr.a * c1, tr.f + tr.e * r1, tr.f + tr.e * r0]
    fig, ax = plt.subplots(1, 3, figsize=(15, 5.2))
    im = ax[0].imshow(diff[r0:r1, c0:c1], cmap="viridis", vmin=-2, vmax=14, extent=ext)
    ax[0].contour(np.isfinite(mndwi[r0:r1, c0:c1]) & (mndwi[r0:r1, c0:c1] > 0), levels=[0.5],
                  colors="white", linewidths=0.6, extent=ext, origin="upper")
    plt.colorbar(im, ax=ax[0], fraction=0.04, label="HH minus HV (dB)")
    ax[0].set_title(f"ASAR HH minus HV, {date}\nwhite line: fitted MNDWI = 0")
    cm = matplotlib.colors.ListedColormap(["#e8e4d8", "#2166ac", "#1a9641"])
    ax[1].imshow(cls[r0:r1, c0:c1], cmap=cm, vmin=0, vmax=2, extent=ext, interpolation="nearest")
    ax[1].set_title(f"open water {a_water:,.0f} km2 (blue)\nflooded vegetation {a_flood_near:,.0f} km2 within 10 km (green)")
    bins = np.linspace(-6, 18, 97)
    for s, lab, col in ((1, "open water", "#2166ac"), (2, "belt", "#1a9641"),
                        (3, "dry land", "#a6611a")):
        ax[2].hist(diff[stratum == s], bins=bins, density=True, histtype="step", lw=1.4,
                   color=col, label=lab)
    ax[2].axvline(cross, color="k", lw=0.8, ls="--", label=f"dry-land 95th, {cross:.1f} dB")
    ax[2].set_xlabel("HH minus HV (dB)")
    ax[2].set_ylabel("density")
    ax[2].legend(frameon=False)
    ax[2].set_title("distribution by optical stratum")
    for a in ax[:2]:
        a.set_xlabel("longitude")
        a.set_ylabel("latitude")
    fig.tight_layout()
    fig.savefig(os.path.join(PNG, f"asar_hhhv_test_{date}.png"), dpi=160)
    print(f"  wrote {os.path.join(PNG, f'asar_hhhv_test_{date}.png')}")
    return summary


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "2012-02-22")
