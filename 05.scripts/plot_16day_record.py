#!/usr/bin/env python3
"""The sixteen-day inundation record, drawn.

Two panels. The upper one is the record itself, open water with its ninety-five
per cent credible band and the vegetated fraction stacked on top, across all 936
steps. The band is the point of the figure as much as the line: it is narrow
where MODIS observes every step and wide through the 1980s and 1990s where
Landsat alone supplies about four usable looks a year, and a reader should be
able to see that difference without being told.

The lower panel is why the band behaves that way. It shows how many instruments
actually reported at each step, so the widening of the interval above can be
read directly off the thinning of the observations below.

Writes 03.outputs/PNG/fig08_inundation_16day.png and the PDF twin.
"""

import csv
import datetime
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "03.outputs", "CSV", "chilwa_inundation_16day.csv")

rows = list(csv.DictReader(open(SRC)))
d = np.array([np.datetime64(r["date"]) for r in rows])
ow = np.array([float(r["open_water_km2"]) for r in rows])
lo = np.array([float(r["open_water_lo95"]) for r in rows])
hi = np.array([float(r["open_water_hi95"]) for r in rows])
veg = np.array([float(r["flooded_veg_km2"]) for r in rows])
obs = np.array([r["dominant_source"] == "observed" for r in rows])
nobs = np.array([int(r["n_obs"]) for r in rows])

fig, ax = plt.subplots(2, 1, figsize=(12, 7), sharex=True,
                       gridspec_kw={"height_ratios": [3, 1], "hspace": 0.08})

for x0, x1 in (("1984-01-01", "2000-01-01"), ("2015-01-01", "2025-01-01")):
    ax[0].axvspan(np.datetime64(x0), np.datetime64(x1), color="#1F6E8C", alpha=0.04)

ax[0].fill_between(d, lo, hi, color="#1F6E8C", alpha=0.22, lw=0,
                   label="open water, 95 per cent credible interval")
ax[0].fill_between(d, ow, ow + veg, color="#5F8836", alpha=0.55, lw=0,
                   label="water standing within vegetation")
ax[0].plot(d, ow, color="#12455A", lw=1.0, label="open water")
ax[0].set_ylabel("inundated area (km$^2$)")
ax[0].set_ylim(0, max(hi.max(), (ow + veg).max()) * 1.05)
ax[0].legend(frameon=False, fontsize=9, loc="upper left", ncol=3)
ax[0].set_title("Lake Chilwa inundation at sixteen-day steps, 1984 to 2024, "
                "fused from Landsat, MODIS, Envisat ASAR and ALOS PALSAR", fontsize=11)

for y, lab in ((1995, "1995\nrecession"), (2012, "2012\nfieldwork"),
               (2018, "2018\nnear-desiccation"), (2023, "2023-24\nrefill")):
    i = int(np.argmin(np.abs(d - np.datetime64(f"{y}-06-01"))))
    ax[0].annotate(lab, xy=(d[i], ow[i]), xytext=(d[i], ow[i] + 620),
                   ha="center", fontsize=8, color="#3A4A52",
                   arrowprops=dict(arrowstyle="-", color="#8A9AA2", lw=0.7))

ax[1].fill_between(d, 0, nobs, step="mid", color="#B08447", alpha=0.65, lw=0)
ax[1].plot(d, np.where(obs, np.nan, 0), ".", ms=2.5, color="#A6483A")
ax[1].set_ylabel("scenes\nbehind\neach step", fontsize=9)
ax[1].set_xlabel("")
ax[1].text(np.datetime64("1985-06-01"), nobs.max() * 0.72,
           "red marks along the axis are steps with no clear observation,\n"
           "where the record is carried by the model", fontsize=8, color="#A6483A")
ax[1].xaxis.set_major_locator(mdates.YearLocator(5))
ax[1].xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
for a in ax:
    a.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
for ext in ("PNG", "PDF"):
    out = os.path.join(ROOT, "03.outputs", ext, f"fig08_inundation_16day.{ext.lower()}")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=170 if ext == "PNG" else None)
    print(f"wrote {out}")
