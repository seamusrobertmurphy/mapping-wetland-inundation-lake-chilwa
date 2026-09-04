#!/usr/bin/env python3
"""Assemble the Lake Chilwa inundation time series, 1984 to 2024.

One row a year, carrying open water, water standing within vegetation, their
sum, an uncertainty on each, the sensors behind them and the method that
produced them. Nothing here is computed: every number is read from a file that
a named script wrote, and the provenance record names each of those files with
its checksum, so any value in the output can be traced to the code that made
it.

Sources
  harmonic_open_water_annual.csv    open water, per-pixel harmonic fit of MNDWI
                                    on an adaptive window, read at 24 dates
  sma_harmonic_fractions.csv        the four-endmember mixture model carried
                                    through the same harmonic, giving a
                                    sub-pixel vegetated fraction and a second,
                                    independent open-water estimate
  palsar_flooded_vegetation.csv     L-band double bounce, 2007 to 2010 and
                                    2015 to 2024
  asar_series_join.csv              Envisat, the 2011 and 2012 fieldwork window

Which number goes in which column, and why. Open water is the index threshold
estimate, because that is the quantity the surface-water literature reports and
the one the 2012 work validated against Water Observations from Space. The
vegetated column is the adopted estimate: the radar measurement wherever a
horizontally polarised radar record exists, since hypothesis H2 holds that
optical reflectance cannot see water beneath a reed stand, and the optical
mixture estimate elsewhere. The raw optical estimate is kept in its own column
for every year, so the difference between the two methods is visible in the
fourteen years that carry both rather than hidden by the choice.

Run with --check to rebuild into a temporary file and compare, which is how the
acceptance test proves the file is reproducible.

Writes 03.outputs/CSV/chilwa_inundation_timeseries.csv
       03.outputs/CSV/chilwa_inundation_validation.csv
       03.outputs/JSON/chilwa_inundation_provenance.json
"""

import argparse
import csv
import datetime
import hashlib
import json
import math
import os
import statistics as st
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CSVD = os.path.join(ROOT, "03.outputs", "CSV")
JSOND = os.path.join(ROOT, "03.outputs", "JSON")
TARGET = os.path.join(CSVD, "chilwa_inundation_timeseries.csv")
VALID = os.path.join(CSVD, "chilwa_inundation_validation.csv")
PROV = os.path.join(JSOND, "chilwa_inundation_provenance.json")

SRC = {
    "open_water": os.path.join(CSVD, "harmonic_open_water_annual.csv"),
    "mixture": os.path.join(CSVD, "sma_harmonic_fractions.csv"),
    "palsar": os.path.join(CSVD, "palsar_flooded_vegetation.csv"),
    "asar": os.path.join(CSVD, "asar_series_join.csv"),
    "modis": os.path.join(CSVD, "modis_water_fraction_8day.csv"),
    "survey": os.path.join(CSVD, "optical_archive_survey.csv"),
}
YEARS = list(range(1984, 2025))
BASIN_KM2 = 8752.0

# The Envisat dates whose track actually covers the lake. The other track sees
# 57 per cent of the lake vicinity at 21 degrees, where open water is not dark
# and the double-bounce rule is weak, so it is not used for an annual figure.
ASAR_USABLE = {2011: ["2011-12-24"], 2012: ["2012-02-22"]}

COLUMNS = ["year", "open_water_km2", "open_water_sd_km2", "flooded_veg_km2",
           "flooded_veg_sd_km2", "total_inundation_km2", "total_sd_km2",
           "flooded_veg_optical_km2", "flooded_veg_source", "open_water_mixture_km2",
           "n_clear_obs", "optical_sensors", "radar_sensor", "radar_flooded_km2",
           "radar_date", "optical_flooded_same_date_km2",
           "window_start", "window_end", "model_order", "basin_coverage",
           "method", "quality", "note"]


def rows(path, key="year"):
    if not os.path.exists(path):
        return {}
    out = {}
    for r in csv.DictReader(open(path)):
        v = (r.get(key) or "").strip().strip('"')
        if v[:4].isdigit():
            out.setdefault(int(v[:4]), []).append(r)
    return out


def f(r, k):
    v = (r or {}).get(k)
    if v is None or str(v).strip() == "":
        return None
    try:
        x = float(v)
    except ValueError:
        return None
    return None if math.isnan(x) else x


def quality(cov, obs, hw):
    if cov is None or obs is None:
        return "low"
    if cov >= 0.99 and obs >= 25 and hw <= 2:
        return "high"
    if cov >= 0.95 and obs >= 15:
        return "medium"
    return "low"


def build():
    ow = {y: v[0] for y, v in rows(SRC["open_water"]).items()}
    mx = {y: v[0] for y, v in rows(SRC["mixture"]).items()}
    pal = {y: v[0] for y, v in rows(SRC["palsar"]).items()}
    asar = rows(SRC["asar"], "date")
    survey = rows(SRC["survey"])

    missing = [y for y in YEARS if y not in ow]
    if missing:
        print(f"open-water series is short by {len(missing)} years: {missing}", file=sys.stderr)

    out = []
    for y in YEARS:
        o, m, p = ow.get(y), mx.get(y), pal.get(y)
        if o is None:
            continue
        open_km2, open_sd = f(o, "open_water_mean_km2"), f(o, "open_water_sd_km2")
        cov, obs = f(o, "basin_coverage"), f(o, "median_obs_per_pixel")
        hw = f(o, "half_window")

        opt_veg = f(m, "flooded_km2") if m else None
        opt_veg_sd = f(m, "flooded_sd_km2") if m else None
        mix_open = f(m, "open_km2") if m else None

        radar_sensor, radar_km2, radar_date, same_date = "", None, "", None
        if p is not None:
            radar_sensor = "ALOS PALSAR yearly mosaic"
            radar_km2 = f(p, "flooded_veg_belt_km2")
            radar_date = p.get("mosaic_date", "")
            same_date = f(p, "optical_flooded_same_date_km2")
        elif y in ASAR_USABLE:
            vals = []
            for d in ASAR_USABLE[y]:
                for r in asar.get(y, []):
                    if r.get("date") == d:
                        v = f(r, "flooded_km2_within_10km")
                        if v is not None:
                            vals.append(v)
            if vals:
                radar_sensor = "Envisat ASAR"
                radar_km2 = st.mean(vals)
                radar_date = ", ".join(ASAR_USABLE[y])

        if radar_km2 is not None:
            veg, veg_src = radar_km2, radar_sensor
            # the radar rule admits a few per cent of dry land, and C-band and
            # L-band are both attenuated by a dense stand, so the interval is
            # taken as a quarter of the value rather than from a fit residual
            veg_sd = 0.25 * radar_km2
        elif opt_veg is not None:
            veg, veg_src, veg_sd = opt_veg, "optical mixture model", opt_veg_sd
        else:
            veg, veg_src, veg_sd = None, "", None

        note = ""
        n_year = f(o, "n_scenes_in_year")
        if n_year is not None and n_year <= 1:
            note = ("no clear scene of its own; the level is set by the window and the "
                    "drift term, and the interval is correspondingly wide")
        elif hw is not None and hw >= 5:
            note = (f"carried by a {int(2 * hw + 1)}-year window because the archive is thin "
                    f"here; adjacent years are not independent")
        if y in ASAR_USABLE and radar_sensor == "Envisat ASAR":
            note = ((note + "; ") if note else "") + (
                f"the radar figure is a single Envisat date, "
                f"{', '.join(ASAR_USABLE[y])}, not an annual mean")

        total = (open_km2 + veg) if (open_km2 is not None and veg is not None) else None
        total_sd = (math.hypot(open_sd or 0, veg_sd or 0)
                    if total is not None else None)

        out.append({
            "year": y,
            "open_water_km2": round(open_km2, 1) if open_km2 is not None else "",
            "open_water_sd_km2": round(open_sd, 1) if open_sd else "",
            "flooded_veg_km2": round(veg, 1) if veg is not None else "",
            "flooded_veg_sd_km2": round(veg_sd, 1) if veg_sd else "",
            "total_inundation_km2": round(total, 1) if total is not None else "",
            "total_sd_km2": round(total_sd, 1) if total_sd else "",
            "flooded_veg_optical_km2": round(opt_veg, 1) if opt_veg is not None else "",
            "flooded_veg_source": veg_src,
            "open_water_mixture_km2": round(mix_open, 1) if mix_open is not None else "",
            "n_clear_obs": int(obs) if obs is not None else "",
            "optical_sensors": o.get("sensors", ""),
            "radar_sensor": radar_sensor,
            "radar_flooded_km2": round(radar_km2, 1) if radar_km2 is not None else "",
            "radar_date": radar_date,
            "optical_flooded_same_date_km2": round(same_date, 1) if same_date is not None else "",
            "window_start": o.get("window_start", ""),
            "window_end": o.get("window_end", ""),
            "model_order": o.get("model_order", ""),
            "basin_coverage": cov if cov is not None else "",
            "method": o.get("method", ""),
            "quality": quality(cov, obs, hw),
            "note": note,
        })
    return out


def validation(built):
    """Errors that can be computed from what is on disk, in km2 where possible."""
    v = []
    pairs = [(r["open_water_km2"], r["open_water_mixture_km2"]) for r in built
             if r["open_water_km2"] != "" and r["open_water_mixture_km2"] != ""]
    if pairs:
        d = [abs(a - b) for a, b in pairs]
        v.append({"test": "open water, index threshold against sub-pixel mixture model",
                  "held_out": f"{len(pairs)} years, two independent retrievals",
                  "metric": "mean absolute difference", "value": round(st.mean(d), 1),
                  "units": "km2"})
        v.append({"test": "open water, index threshold against sub-pixel mixture model",
                  "held_out": f"{len(pairs)} years, two independent retrievals",
                  "metric": "root mean square difference",
                  "value": round(math.sqrt(sum(x * x for x in d) / len(d)), 1),
                  "units": "km2"})

    per = {}
    if os.path.exists(SRC["modis"]):
        for r in csv.DictReader(open(SRC["modis"])):
            try:
                yy = int(r["year"]); wf = float(r["water_fraction"])
                nc = float(r["n_clear_cells"] or 0)
            except (ValueError, TypeError, KeyError):
                continue
            if nc > 0 and not math.isnan(wf):
                per.setdefault(yy, []).append(wf)
    mod = {y: st.mean(x) * BASIN_KM2 for y, x in per.items() if len(x) >= 20}
    ours = {r["year"]: r["open_water_km2"] for r in built if r["open_water_km2"] != ""}
    common = sorted(y for y in mod if y in ours)
    if len(common) >= 10:
        a = [ours[y] for y in common]; b = [mod[y] for y in common]
        ma, mb = st.mean(a), st.mean(b)
        cov_ = sum((x - ma) * (y - mb) for x, y in zip(a, b))
        den = math.sqrt(sum((x - ma) ** 2 for x in a) * sum((y - mb) ** 2 for y in b)) or 1e-9
        v.append({"test": "open water against the MODIS eight-day water fraction",
                  "held_out": f"{len(common)} years, a different sensor and grain",
                  "metric": "Pearson r", "value": round(cov_ / den, 3), "units": "correlation"})
        v.append({"test": "open water against the MODIS eight-day water fraction",
                  "held_out": f"{len(common)} years, a different sensor and grain",
                  "metric": "root mean square error",
                  "value": round(math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)) / len(a)), 1),
                  "units": "km2"})

    # like for like: the optical model read on the day the radar looked
    sd_pairs = [(r["optical_flooded_same_date_km2"], r["radar_flooded_km2"]) for r in built
                if r["optical_flooded_same_date_km2"] != "" and r["radar_flooded_km2"] != ""]
    if sd_pairs:
        d = [a - b for a, b in sd_pairs]
        v.append({"test": "vegetated water, optical mixture model against radar double bounce, "
                          "both read on the radar acquisition date",
                  "held_out": f"{len(sd_pairs)} years carrying both",
                  "metric": "mean optical minus radar", "value": round(st.mean(d), 1),
                  "units": "km2"})
        v.append({"test": "vegetated water, optical mixture model against radar double bounce, "
                          "both read on the radar acquisition date",
                  "held_out": f"{len(sd_pairs)} years carrying both",
                  "metric": "standard deviation of the difference",
                  "value": round(st.pstdev(d) if len(d) > 1 else 0.0, 1), "units": "km2"})

    vp = [(r["flooded_veg_optical_km2"], r["radar_flooded_km2"]) for r in built
          if r["flooded_veg_optical_km2"] != "" and r["radar_flooded_km2"] != ""]
    if vp:
        d = [a - b for a, b in vp]
        v.append({"test": "vegetated water, optical annual mean against radar; the radar date "
                          "sits in the dry season so this pair spans different times of year",
                  "held_out": f"{len(vp)} years carrying both",
                  "metric": "mean optical minus radar", "value": round(st.mean(d), 1),
                  "units": "km2"})
    return v


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write(built, check=False):
    tgt = TARGET + (".check" if check else "")
    with open(tgt, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS)
        w.writeheader()
        w.writerows(built)
    if check:
        same = open(tgt, "rb").read() == open(TARGET, "rb").read()
        os.remove(tgt)
        if not same:
            print("rebuild differs from the committed file", file=sys.stderr)
            return 1
        print("rebuild is byte identical")
        return 0

    v = validation(built)
    with open(VALID, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["test", "held_out", "metric", "value", "units"])
        w.writeheader()
        w.writerows(v)

    os.makedirs(JSOND, exist_ok=True)
    prov = {
        "product": "Lake Chilwa inundation time series, 1984 to 2024",
        "builder": "05.scripts/build_inundation_timeseries.py",
        "built": datetime.date.today().isoformat(),
        "rows": len(built),
        "basin_km2": BASIN_KM2,
        "inputs": [{"role": k, "path": os.path.relpath(p, ROOT), "sha256": sha(p),
                    "bytes": os.path.getsize(p)}
                   for k, p in sorted(SRC.items()) if os.path.exists(p)],
        "producers": {
            "harmonic_open_water_annual.csv": "05.scripts/harmonic_open_water_series.py",
            "sma_harmonic_fractions.csv": "05.scripts/sma_harmonic_flooded_veg.py",
            "palsar_flooded_vegetation.csv": "05.scripts/palsar_flooded_vegetation.py",
            "asar_series_join.csv": "05.scripts/asar_series_join.py",
            "optical_archive_survey.csv": "05.scripts/survey_optical_archive.py",
            "modis_water_fraction_8day.csv": "05.scripts/modis_gap_fill_2012.py",
        },
    }
    with open(PROV, "w") as fh:
        json.dump(prov, fh, indent=2)
    print(f"wrote {TARGET}  ({len(built)} rows)")
    print(f"wrote {VALID}  ({len(v)} tests)")
    print(f"wrote {PROV}   ({len(prov['inputs'])} inputs with checksums)")
    for t in v:
        print(f"  {t['metric']}: {t['value']} {t['units']}  ({t['test']})")
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true",
                    help="rebuild and compare, without writing over the file")
    a = ap.parse_args()
    built = build()
    if not built:
        print("no rows built; the source series are not ready", file=sys.stderr)
        return 1
    return write(built, check=a.check)


if __name__ == "__main__":
    sys.exit(main())
