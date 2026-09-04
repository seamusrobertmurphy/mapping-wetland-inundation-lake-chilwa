#!/usr/bin/env python3
"""Exploratory statistics for the multi-source inundation record.

Three tables, written both as CSV for the manuscript to read at render time and
as GitHub-flavoured markdown for the README.

Table 1 is the inventory. What each instrument is, when it flew, what it
resolves, which bands the water retrieval uses, and how much of the finished
record it actually carries. The last column is the one that matters for bias:
a stream contributing four per cent of the steps cannot be blamed or credited
for much, and a stream contributing sixty per cent sets the record's character.

Table 2 is the distributional check that has to happen before any model is
fitted. For each stream it gives the centre and spread, the shape through
skewness and excess kurtosis, a formal normality test, and then the regression
of that stream against Landsat, which is the reference the others are calibrated
to. The intercept is the additive bias in square kilometres and the slope is the
multiplicative one; a slope away from one means the stream stretches or
compresses the range rather than merely sitting high or low, and that distorts
every downstream fit differently at the wet and dry ends. Reporting the two
separately is the point, because a single correlation hides both.

Table 3 is the record by era, which is where the honest limits live. The share
of steps carrying a real observation and the mean width of the credible interval
say plainly how much of each era is measured and how much is carried.

Writes 03.outputs/CSV/eda_table1_inventory.csv, _table2_distributions.csv,
       _table3_by_era.csv and the matching .md files in 03.outputs/TABLES/
"""

import csv
import datetime
import math
import os
import statistics as st

import numpy as np
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CSVD = os.path.join(ROOT, "03.outputs", "CSV")
TBL = os.path.join(ROOT, "03.outputs", "TABLES")
os.makedirs(TBL, exist_ok=True)

# Native grain and the bands each water retrieval actually reads. Counts of
# images over the basin come from the archive survey, which queried the
# catalogue against the basin geometry rather than a bounding box.
INVENTORY = [
    ("Landsat 4 TM", "1988-1992", "30 m", "16 d", "green 0.52-0.60, SWIR1 1.55-1.75", "MNDWI"),
    ("Landsat 5 TM", "1984-2011", "30 m", "16 d", "green 0.52-0.60, SWIR1 1.55-1.75", "MNDWI"),
    ("Landsat 7 ETM+", "1999-2024", "30 m", "16 d", "green 0.52-0.60, SWIR1 1.55-1.75", "MNDWI"),
    ("Landsat 8 OLI", "2013-2024", "30 m", "16 d", "green 0.53-0.59, SWIR1 1.57-1.65", "MNDWI"),
    ("Landsat 9 OLI-2", "2021-2024", "30 m", "16 d", "green 0.53-0.59, SWIR1 1.57-1.65", "MNDWI"),
    ("MODIS Terra MOD09A1", "2000-2024", "500 m", "8 d", "green 0.545-0.565, SWIR1 1.628-1.652", "MNDWI"),
    ("MODIS Aqua MYD09A1", "2002-2024", "500 m", "8 d", "green 0.545-0.565, SWIR1 1.628-1.652", "MNDWI"),
    ("Envisat ASAR APP", "2011-2012", "30 m", "35 d", "C-band 5.3 GHz, HH and HV", "HH minus HV"),
    ("ALOS PALSAR mosaic", "2007-2024", "25 m", "annual", "L-band 1.27 GHz, HH and HV", "HH minus HV"),
]


def load_streams():
    """Per-step water area by instrument family, and the fused record."""
    obs = {}
    for r in csv.DictReader(open(os.path.join(CSVD, "optical_water_16day.csv"))):
        cf = float(r["clear_fraction"] or 0)
        if cf >= 0.30:
            obs.setdefault(r["sensor_family"], {})[r["date"]] = float(r["water_km2"]) / cf
    fused = {r["date"]: float(r["open_water_km2"])
             for r in csv.DictReader(open(os.path.join(CSVD, "chilwa_inundation_16day.csv")))}
    return obs, fused


def counts():
    """Images over the basin, from the archive survey."""
    n = {}
    p = os.path.join(CSVD, "optical_archive_survey.csv")
    if os.path.exists(p):
        for r in csv.DictReader(open(p)):
            n[r["sensor"]] = n.get(r["sensor"], 0) + int(r["n_cloud_lt_101"])
    return n


def md_table(rows, cols, path, title):
    with open(path, "w") as fh:
        fh.write(f"**{title}**\n\n")
        fh.write("| " + " | ".join(cols) + " |\n")
        fh.write("|" + "|".join([":---"] * len(cols)) + "|\n")
        for r in rows:
            fh.write("| " + " | ".join(str(r.get(c, "")) for c in cols) + " |\n")


def table1(obs, fused):
    surv = counts()
    key = {"Landsat 4 TM": "L4_TM", "Landsat 5 TM": "L5_TM", "Landsat 7 ETM+": "L7_ETM",
           "Landsat 8 OLI": "L8_OLI", "Landsat 9 OLI-2": "L9_OLI"}
    n_ls = len(obs.get("landsat", {})); n_md = len(obs.get("modis", {}))
    rows = []
    for name, span, res, rev, bands, idx in INVENTORY:
        if name in key:
            n_img = surv.get(key[name], 0); contrib = "part of the Landsat stream"
        elif "MODIS" in name:
            n_img = ""; contrib = f"{n_md} steps ({100 * n_md / 936:.0f}% of the record)"
        elif "ASAR" in name:
            n_img = 40; contrib = "7 usable dates, vegetated class only"
        else:
            n_img = 14; contrib = "14 annual mosaics, vegetated class only"
        rows.append({"Instrument": name, "Span over the basin": span, "Grain": res,
                     "Repeat": rev, "Bands read": bands, "Water rule": idx,
                     "Scenes over basin": n_img if n_img != "" else "-",
                     "Contribution to the record": contrib})
    rows.append({"Instrument": "**Landsat, all platforms**", "Span over the basin": "1984-2024",
                 "Grain": "30 m", "Repeat": "16 d", "Bands read": "green, SWIR1",
                 "Water rule": "MNDWI > 0", "Scenes over basin": sum(surv.get(v, 0) for v in key.values()),
                 "Contribution to the record": f"{n_ls} steps ({100 * n_ls / 936:.0f}% of the record)"})
    cols = ["Instrument", "Span over the basin", "Grain", "Repeat", "Bands read",
            "Water rule", "Scenes over basin", "Contribution to the record"]
    with open(os.path.join(CSVD, "eda_table1_inventory.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols); w.writeheader(); w.writerows(rows)
    md_table(rows, cols, os.path.join(TBL, "table1_inventory.md"),
             "Table 1. The instruments behind the sixteen-day inundation record, 1984 to 2024.")
    return rows, cols


def table2(obs, fused):
    ref = obs.get("landsat", {})
    rows = []
    series = [("Landsat, 30 m", obs.get("landsat", {})),
              ("MODIS, 500 m", obs.get("modis", {})),
              ("Fused record", fused)]
    for name, d in series:
        v = np.array(list(d.values()), dtype=float)
        if len(v) < 10:
            continue
        sk = float(stats.skew(v)); ku = float(stats.kurtosis(v))
        samp = v if len(v) <= 5000 else np.random.default_rng(1).choice(v, 5000, replace=False)
        W, pW = stats.shapiro(samp)
        row = {"Stream": name, "n": len(v),
               "Mean (km2)": f"{v.mean():,.0f}", "SD (km2)": f"{v.std(ddof=1):,.0f}",
               "Median (km2)": f"{np.median(v):,.0f}",
               "IQR (km2)": f"{np.percentile(v, 75) - np.percentile(v, 25):,.0f}",
               "Skewness": f"{sk:+.2f}", "Excess kurtosis": f"{ku:+.2f}",
               "Shapiro-Wilk W": f"{W:.3f}",
               "p (normality)": ("< 0.001" if pW < 0.001 else f"{pW:.3f}")}
        common = sorted(set(d) & set(ref))
        if name != "Landsat, 30 m" and len(common) >= 30:
            x = np.array([ref[k] for k in common]); yv = np.array([d[k] for k in common])
            slope, inter, r, p, se = stats.linregress(yv, x)
            resid = x - (slope * yv + inter)
            row.update({"n paired with Landsat": len(common),
                        "Intercept (km2)": f"{inter:+,.0f}",
                        "Slope": f"{slope:.3f}",
                        "Slope 95% CI": f"{slope - 1.96 * se:.3f} to {slope + 1.96 * se:.3f}",
                        "r": f"{r:.3f}",
                        "Residual SD (km2)": f"{resid.std(ddof=1):,.0f}",
                        "Slope differs from 1": ("yes" if abs(slope - 1) > 1.96 * se else "no")})
        else:
            row.update({"n paired with Landsat": "reference", "Intercept (km2)": "0 by definition",
                        "Slope": "1 by definition", "Slope 95% CI": "-", "r": "-",
                        "Residual SD (km2)": "-", "Slope differs from 1": "-"})
        rows.append(row)
    cols = ["Stream", "n", "Mean (km2)", "SD (km2)", "Median (km2)", "IQR (km2)",
            "Skewness", "Excess kurtosis", "Shapiro-Wilk W", "p (normality)",
            "n paired with Landsat", "Intercept (km2)", "Slope", "Slope 95% CI", "r",
            "Residual SD (km2)", "Slope differs from 1"]
    with open(os.path.join(CSVD, "eda_table2_distributions.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols); w.writeheader(); w.writerows(rows)
    md_table(rows, cols, os.path.join(TBL, "table2_distributions.md"),
             "Table 2. Distribution of each observation stream and its bias against the "
             "Landsat reference. The intercept is additive bias, the slope multiplicative.")
    return rows, cols


def table3():
    rec = list(csv.DictReader(open(os.path.join(CSVD, "chilwa_inundation_16day.csv"))))
    rows = []
    for era in ("1984-1999", "2000-2014", "2015-2024"):
        sub = [r for r in rec if r["era"] == era]
        ow = np.array([float(r["open_water_km2"]) for r in sub])
        w = np.array([float(r["open_water_hi95"]) - float(r["open_water_lo95"]) for r in sub])
        seen = sum(1 for r in sub if r["dominant_source"] == "observed")
        rows.append({"Era": era, "Steps": len(sub),
                     "Observed": f"{seen} ({100 * seen / len(sub):.0f}%)",
                     "Carried by model": f"{len(sub) - seen} ({100 * (len(sub) - seen) / len(sub):.0f}%)",
                     "Mean open water (km2)": f"{ow.mean():,.0f}",
                     "SD (km2)": f"{ow.std(ddof=1):,.0f}",
                     "Minimum (km2)": f"{ow.min():,.0f}", "Maximum (km2)": f"{ow.max():,.0f}",
                     "Mean 95% interval width (km2)": f"{w.mean():,.0f}"})
    cols = list(rows[0].keys())
    with open(os.path.join(CSVD, "eda_table3_by_era.csv"), "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=cols); wr.writeheader(); wr.writerows(rows)
    md_table(rows, cols, os.path.join(TBL, "table3_by_era.md"),
             "Table 3. The finished record by era, showing how much of each is measured "
             "and how much is carried by the model.")
    return rows, cols


def main():
    obs, fused = load_streams()
    r1, c1 = table1(obs, fused)
    r2, c2 = table2(obs, fused)
    r3, c3 = table3()
    for name, rows, cols in (("Table 1 inventory", r1, c1),
                             ("Table 2 distributions", r2, c2),
                             ("Table 3 by era", r3, c3)):
        print(f"\n=== {name}")
        for r in rows:
            print("  " + " | ".join(f"{c}={r.get(c,'')}" for c in cols[:6]))
    print(f"\nwrote three CSVs to 03.outputs/CSV and three markdown tables to 03.outputs/TABLES")


if __name__ == "__main__":
    main()
