#!/usr/bin/env python3
"""Join the Envisat ASAR series to the optical harmonic fit on every radar
date, 5 November 2011 to 15 March 2012.

Each date runs the same test as asar_hhhv_test.py: the optical model read on
that date supplies open water, the radar supplies flooded vegetation by the
double-bounce rule with thresholds set from that date's own dry land, and the
lake vicinity's radar coverage is recorded because the four tracks place the
swath differently. The radar is never fitted through time; sixteen dates over
150 days cannot constrain an annual model.

Writes 03.outputs/CSV/asar_series_join.csv and
       03.outputs/PNG/asar_series_join.png
"""

import csv
import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np

import asar_hhhv_test as t

DATES = sorted({os.path.basename(f)[14:22] for f in
                glob.glob(os.path.join(t.TIF, "asar", "ASA_APP_1PNESA*_sigma0.tif"))})
DATES = [f"{d[:4]}-{d[4:6]}-{d[6:]}" for d in DATES]


def main():
    rows = []
    for d in DATES:
        try:
            rows.append(t.main(d, figure=(d == "2012-02-22")))
        except SystemExit as e:
            print(d, e)
    out = os.path.join(t.CSV, "asar_series_join.csv")
    keys = max((list(r.keys()) for r in rows), key=len)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out}")

    dates = np.array([np.datetime64(r["date"]) for r in rows])
    ow_all = np.array([r["optical_open_water_km2_gt_0.10"] for r in rows])
    cov = np.array([r["vicinity_covered"] for r in rows])
    has = np.array(["flooded_km2_within_10km" in r for r in rows])
    ow = np.array([r.get("open_water_km2_mndwi_gt_0.10", np.nan) for r in rows])
    fl = np.array([r.get("flooded_km2_within_10km", np.nan) for r in rows])
    share = np.array([r.get("flooded_share_of_nonwater_within_10km", np.nan) for r in rows])
    fpr = [r["dry_land_false_positive"] for r in rows if "dry_land_false_positive" in r]
    full = has & (cov >= 0.90)
    part = has & (cov < 0.90)
    fig, ax = plt.subplots(2, 1, figsize=(8, 6.5), sharex=True)
    ax[0].plot(dates, ow_all, "o-", color="#2166ac", label="open water, optical fit (MNDWI > 0.10)")
    ax[0].plot(dates[full], (ow + fl)[full], "s", color="#1a9641",
               label="open water plus flooded vegetation within 10 km, radar covers 90 per cent or more")
    ax[0].plot(dates[part], (ow + fl)[part], "s", mfc="none", color="#1a9641",
               label="same, radar covers 25 to 90 per cent (area only where it looked)")
    ax[0].set_ylabel("km$^2$")
    ax[0].legend(frameon=False, fontsize=7.5)
    ax[1].bar(dates[has], 100 * share[has], width=2.5,
              color=np.where(cov[has] >= 0.90, "#1a9641", "#a6dba0"))
    ax[1].axhline(100 * np.mean(fpr), color="k",
                  lw=0.8, ls="--", label="mean false-positive rate on dry land")
    ax[1].set_ylabel("flooded share of non-water\nwithin 10 km (%)")
    ax[1].legend(frameon=False, fontsize=8)
    ax[1].xaxis.set_major_formatter(mdates.DateFormatter("%d %b\n%Y"))
    fig.tight_layout()
    png = os.path.join(t.PNG, "asar_series_join.png")
    fig.savefig(png, dpi=160)
    print(f"wrote {png}")


if __name__ == "__main__":
    main()
