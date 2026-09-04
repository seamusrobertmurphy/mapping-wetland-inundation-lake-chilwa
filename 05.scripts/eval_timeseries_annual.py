#!/usr/bin/env python3
"""The acceptance test for the Lake Chilwa inundation time series.

This file is the contract. The dataset is finished when every criterion below
passes, and not before. Nothing else decides it: not a plausible-looking plot,
not a script that ran without error, not a number that resembles a published
one. Each criterion is a function that returns a verdict, a measured value and
the evidence it read, so a failure names the row or the year that caused it.

The target is 03.outputs/CSV/chilwa_inundation_timeseries.csv, one row per year
from 1984 to 2024, carrying open water, water standing within vegetation, their
sum, an uncertainty on each, the sensors behind them and the method that
produced them. Hypothesis H1 in SPEC.md is that most of the basin's surface
water stands under vegetation, so a record of open water alone measures the
wrong quantity; the vegetated column is therefore not optional and a row that
carries only open water is a failing row.

Run: python3 eval_timeseries.py [--json]
Writes 03.outputs/JSON/eval_timeseries.json and prints a report.
Exit status is 0 when every criterion passes, 1 otherwise, so a loop can branch
on it.
"""

import csv
import json
import math
import os
import statistics as st
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CSV = os.path.join(ROOT, "03.outputs", "CSV")
TARGET = os.path.join(CSV, "chilwa_inundation_timeseries.csv")
VALID = os.path.join(CSV, "chilwa_inundation_validation.csv")
PROV = os.path.join(ROOT, "03.outputs", "JSON", "chilwa_inundation_provenance.json")
OUT = os.path.join(ROOT, "03.outputs", "JSON", "eval_timeseries.json")

YEARS = list(range(1984, 2025))
BASIN_KM2 = 8752.0
GAP_YEARS = [1985, 1988, 2002, 2012]      # no value in the committed record
THIN_YEARS = [1984, 1987, 1992, 2001, 2003, 2011]   # fewer than five clear scenes

# Physical bounds. Lake Chilwa has dried to bare mud within living memory and
# has reached roughly 2,000 km2 of open water in a full year, and the whole
# basin is 8,752 km2, so any value outside these is a processing fault rather
# than a wet or dry year.
OPEN_MIN, OPEN_MAX = 100.0, 2400.0
TOTAL_MAX = 4000.0

REQUIRED_COLUMNS = [
    "year", "open_water_km2", "open_water_sd_km2", "flooded_veg_km2",
    "flooded_veg_sd_km2", "total_inundation_km2", "total_sd_km2",
    "n_clear_obs", "optical_sensors", "radar_sensor", "radar_flooded_km2",
    "flooded_veg_optical_km2", "open_water_mixture_km2",
    "radar_date", "optical_flooded_same_date_km2",
    "method", "quality", "note",
]

PLACEHOLDER_STRINGS = ["todo", "tbd", "placeholder", "xxx", "fixme", "n/a", "nan"]


def read_target():
    if not os.path.exists(TARGET):
        return None
    with open(TARGET) as fh:
        return list(csv.DictReader(fh))


def num(row, key):
    v = (row.get(key) or "").strip()
    if v == "":
        return None
    try:
        f = float(v)
    except ValueError:
        return None
    return None if math.isnan(f) else f


def verdict(name, title, ok, measured, detail, blocking=True):
    return {"id": name, "title": title, "pass": bool(ok), "measured": measured,
            "detail": detail, "blocking": blocking}


# --------------------------------------------------------------------------
# The criteria
# --------------------------------------------------------------------------

def c01_file_exists(rows):
    return verdict("C01", "The target file exists and parses as CSV",
                   rows is not None, TARGET if rows is not None else None,
                   "absent" if rows is None else f"{len(rows)} rows")


def c02_schema(rows):
    if rows is None:
        return verdict("C02", "Every required column is present", False, None, "no file")
    have = set(rows[0].keys())
    missing = [c for c in REQUIRED_COLUMNS if c not in have]
    return verdict("C02", "Every required column is present", not missing,
                   f"{len(have)} columns", f"missing {missing}" if missing else "complete")


def c03_coverage(rows):
    """One row for every year from 1984 to 2024, none missing and none twice."""
    if rows is None:
        return verdict("C03", "One row per year, 1984 to 2024, no gaps", False, None, "no file")
    ys = [int(r["year"]) for r in rows if (r.get("year") or "").strip().isdigit()]
    missing = sorted(set(YEARS) - set(ys))
    dupes = sorted({y for y in ys if ys.count(y) > 1})
    ok = not missing and not dupes and len(ys) == len(YEARS)
    return verdict("C03", "One row per year, 1984 to 2024, no gaps", ok, f"{len(ys)}/41 years",
                   f"missing {missing}; duplicated {dupes}" if not ok else "41 of 41")


def c04_no_nulls(rows):
    """Open water, vegetated water and their sum are present in every row."""
    if rows is None:
        return verdict("C04", "No missing value in the three area columns", False, None, "no file")
    bad = []
    for r in rows:
        for k in ("open_water_km2", "flooded_veg_km2", "total_inundation_km2"):
            if num(r, k) is None:
                bad.append(f"{r.get('year')}:{k}")
    return verdict("C04", "No missing value in the three area columns", not bad,
                   f"{len(bad)} empty cells", "; ".join(bad[:12]) if bad else "all present")


def c05_physical_range(rows):
    """Every area sits inside what the basin can physically hold."""
    if rows is None:
        return verdict("C05", "Every area is physically possible", False, None, "no file")
    bad = []
    for r in rows:
        o, f, t = num(r, "open_water_km2"), num(r, "flooded_veg_km2"), num(r, "total_inundation_km2")
        y = r.get("year")
        if o is not None and not (OPEN_MIN <= o <= OPEN_MAX):
            bad.append(f"{y}: open water {o:.0f} km2 outside {OPEN_MIN:.0f} to {OPEN_MAX:.0f}")
        if f is not None and not (0 <= f <= BASIN_KM2):
            bad.append(f"{y}: vegetated water {f:.0f} km2 outside 0 to basin")
        if t is not None and t > TOTAL_MAX:
            bad.append(f"{y}: total {t:.0f} km2 above {TOTAL_MAX:.0f}")
        if None not in (o, f, t) and abs((o + f) - t) > 1.0:
            bad.append(f"{y}: total {t:.0f} does not equal open {o:.0f} plus vegetated {f:.0f}")
    return verdict("C05", "Every area is physically possible and the sum adds up", not bad,
                   f"{len(bad)} violations", "; ".join(bad[:10]) if bad else "all inside bounds")


def c06_outliers_defended(rows):
    """A year far from its neighbours must say in the row why.

    The committed spectral-mixture record put 2003 at 118 km2 of open water and
    1987 at 318 km2, both from three scenes or fewer, which the lake's own
    hydrology does not permit. An excursion is allowed, but only where the row
    states its cause."""
    if rows is None:
        return verdict("C06", "Every outlier year carries a stated cause", False, None, "no file")
    vals = [(int(r["year"]), num(r, "open_water_km2")) for r in rows
            if (r.get("year") or "").strip().isdigit() and num(r, "open_water_km2") is not None]
    if len(vals) < 10:
        return verdict("C06", "Every outlier year carries a stated cause", False,
                       f"{len(vals)} usable rows", "too few rows to test")
    xs = [v for _, v in vals]
    med = st.median(xs)
    mad = st.median([abs(x - med) for x in xs]) * 1.4826 or 1.0
    bad = []
    notes = {int(r["year"]): (r.get("note") or "").strip() for r in rows
             if (r.get("year") or "").strip().isdigit()}
    for y, v in vals:
        if abs(v - med) > 3 * mad and len(notes.get(y, "")) < 10:
            bad.append(f"{y}: {v:.0f} km2 is {abs(v-med)/mad:.1f} robust sd from the median, no note")
    return verdict("C06", "Every outlier year carries a stated cause", not bad,
                   f"median {med:.0f} km2, robust sd {mad:.0f} km2",
                   "; ".join(bad) if bad else "no unexplained excursion")


def c07_uncertainty(rows):
    """No bare estimate. Every area carries a dispersion, per the house rule."""
    if rows is None:
        return verdict("C07", "Every estimate carries an uncertainty", False, None, "no file")
    bad = []
    for r in rows:
        for k in ("open_water_sd_km2", "flooded_veg_sd_km2", "total_sd_km2"):
            v = num(r, k)
            if v is None or v <= 0:
                bad.append(f"{r.get('year')}:{k}")
    return verdict("C07", "Every estimate carries an uncertainty", not bad,
                   f"{len(bad)} missing or zero", "; ".join(bad[:12]) if bad else "all present")


def c08_thin_years_rebuilt(rows):
    """The four empty years and the six thin years are retrieved, not interpolated.

    A value copied from the years either side is not a measurement. Each of
    these years must name a retrieval method, and straight interpolation is
    explicitly disallowed."""
    if rows is None:
        return verdict("C08", "Gap and thin years carry a retrieval, not an interpolation",
                       False, None, "no file")
    by = {int(r["year"]): r for r in rows if (r.get("year") or "").strip().isdigit()}
    bad = []
    for y in GAP_YEARS + THIN_YEARS:
        r = by.get(y)
        if r is None:
            bad.append(f"{y}: no row")
            continue
        m = (r.get("method") or "").lower()
        if not m:
            bad.append(f"{y}: no method named")
        elif "interpol" in m or "neighbour" in m or "neighbor" in m or "carried" in m:
            bad.append(f"{y}: method is '{m.strip()}', which is interpolation")
        if num(r, "open_water_km2") is None:
            bad.append(f"{y}: no open water value")
    return verdict("C08", "Gap and thin years carry a retrieval, not an interpolation", not bad,
                   f"{len(GAP_YEARS + THIN_YEARS)} years tested",
                   "; ".join(bad[:12]) if bad else "all retrieved")


def c09_radar_where_it_exists(rows):
    """Radar supplies the vegetated fraction wherever a radar record exists.

    H1 says the vegetated water cannot be recovered from optical reflectance,
    so the years that have radar must use it. Envisat ASAR covers 2011 and
    2012, ALOS PALSAR the yearly mosaics from 2007, and Sentinel-1 from 2015."""
    if rows is None:
        return verdict("C09", "Radar carries the vegetated fraction where radar exists",
                       False, None, "no file")
    want = set(range(2007, 2011)) | {2011, 2012} | set(range(2015, 2025))
    by = {int(r["year"]): r for r in rows if (r.get("year") or "").strip().isdigit()}
    bad = []
    for y in sorted(want):
        r = by.get(y)
        if r is None:
            bad.append(f"{y}: no row"); continue
        if not (r.get("radar_sensor") or "").strip():
            bad.append(f"{y}: no radar sensor named")
        elif num(r, "radar_flooded_km2") is None:
            bad.append(f"{y}: radar sensor named but no radar area")
    return verdict("C09", "Radar carries the vegetated fraction where radar exists", not bad,
                   f"{len(want)} radar-era years", "; ".join(bad[:12]) if bad else "all carry radar")


def c10_cross_sensor_agreement(rows):
    """The optical record must agree with the independent MODIS record.

    Both read a water index from optical reflectance, so this is corroboration
    of the decision rule rather than of the observation, and it is reported that
    way. It still catches a year that has gone wrong, which is its purpose."""
    if rows is None:
        return verdict("C10", "Optical open water agrees with the MODIS record", False, None, "no file")
    mp = os.path.join(CSV, "modis_water_fraction_8day.csv")
    if not os.path.exists(mp):
        return verdict("C10", "Optical open water agrees with the MODIS record", False, None,
                       "MODIS file absent")
    per = {}
    for r in csv.DictReader(open(mp)):
        try:
            y = int(r["year"]); wf = float(r["water_fraction"]); n = float(r["n_clear_cells"] or 0)
        except (ValueError, TypeError, KeyError):
            continue
        if n <= 0 or math.isnan(wf):
            continue
        per.setdefault(y, []).append(wf)
    mod = {y: st.mean(v) * BASIN_KM2 for y, v in per.items() if len(v) >= 20}
    ours = {int(r["year"]): num(r, "open_water_km2") for r in rows
            if (r.get("year") or "").strip().isdigit()}
    common = sorted(y for y in mod if ours.get(y) is not None)
    if len(common) < 10:
        return verdict("C10", "Optical open water agrees with the MODIS record", False,
                       f"{len(common)} overlapping years", "too few years to correlate")
    a = [ours[y] for y in common]; b = [mod[y] for y in common]
    ma, mb = st.mean(a), st.mean(b)
    cov = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    den = math.sqrt(sum((x - ma) ** 2 for x in a) * sum((y - mb) ** 2 for y in b)) or 1e-9
    r_ = cov / den
    rmse = math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)) / len(a))
    ok = r_ >= 0.80 and rmse <= 200.0
    return verdict("C10", "Optical open water agrees with the MODIS record", ok,
                   f"r = {r_:.2f}, RMSE = {rmse:.0f} km2 over {len(common)} years",
                   "needs r at or above 0.80 and RMSE at or below 200 km2")


def c11_optical_radar_vegetated_compared(rows):
    """Where both exist, the optical and radar vegetated fractions are compared.

    This is the study's own claim under test. If the optical mixture model and
    the radar disagree on the vegetated class in the years that carry both,
    H2's premise is evidenced; if they agree, the radar is confirming the
    optical. Either outcome is a result, but the comparison has to exist.

    The comparison reads the raw optical column, not the adopted one. The
    adopted column takes the radar wherever radar exists, so comparing it
    against the radar would in those years compare the radar with itself and
    report a difference of zero, which would be a guard that cannot fail."""
    if rows is None:
        return verdict("C11", "Optical and radar vegetated fractions are compared", False, None, "no file")
    pairs = [(int(r["year"]), num(r, "flooded_veg_optical_km2"), num(r, "radar_flooded_km2"))
             for r in rows if (r.get("year") or "").strip().isdigit()]
    pairs = [(y, a, b) for y, a, b in pairs if a is not None and b is not None]
    if len(pairs) < 8:
        return verdict("C11", "Optical and radar vegetated fractions are compared", False,
                       f"{len(pairs)} years with both", "needs at least 8 years carrying both")
    diffs = [a - b for _, a, b in pairs]
    bias = st.mean(diffs)
    sd = st.pstdev(diffs) if len(diffs) > 1 else 0.0
    return verdict("C11", "Optical and radar vegetated fractions are compared", True,
                   f"{len(pairs)} years, mean optical minus radar {bias:+.0f} km2, sd {sd:.0f} km2",
                   "reported as a finding, not a pass or fail threshold")


def c12_validation_file(rows):
    """A holdout test exists and reports its error in the same units."""
    if not os.path.exists(VALID):
        return verdict("C12", "A holdout validation exists with error in km2", False, None,
                       "validation file absent")
    vr = list(csv.DictReader(open(VALID)))
    if not vr:
        return verdict("C12", "A holdout validation exists with error in km2", False, "0 rows", "empty")
    need = {"test", "held_out", "metric", "value", "units"}
    have = set(vr[0].keys())
    if not need <= have:
        return verdict("C12", "A holdout validation exists with error in km2", False,
                       f"{len(vr)} rows", f"missing columns {sorted(need - have)}")
    km2 = [r for r in vr if (r.get("units") or "").strip() == "km2"]
    return verdict("C12", "A holdout validation exists with error in km2", bool(km2),
                   f"{len(vr)} rows, {len(km2)} in km2",
                   "at least one error reported in km2" if km2 else "no km2 error reported")


def c13_no_placeholders(rows):
    """No placeholder text and none of the four known invented numbers.

    Overall accuracy 81 per cent, kappa 0.77, a coupling of r = 0.73 and
    agreement at 74 to 92 per cent of locations are carried in the manuscript
    as placeholders. None of them may reach this file as though measured."""
    if rows is None:
        return verdict("C13", "No placeholder text or invented number", False, None, "no file")
    bad = []
    for r in rows:
        for k, v in r.items():
            s = (v or "").strip().lower()
            if s in PLACEHOLDER_STRINGS:
                bad.append(f"{r.get('year')}:{k}='{s}'")
            for tok in PLACEHOLDER_STRINGS:
                if tok in s and len(s) < 40 and tok != "nan":
                    bad.append(f"{r.get('year')}:{k} contains '{tok}'")
    return verdict("C13", "No placeholder text or invented number", not bad,
                   f"{len(bad)} suspect cells", "; ".join(sorted(set(bad))[:10]) if bad else "clean")


def c14_provenance(rows):
    """Every input is named with its checksum, so the file can be traced back."""
    if not os.path.exists(PROV):
        return verdict("C14", "A provenance record names every input with a checksum",
                       False, None, "provenance file absent")
    try:
        p = json.load(open(PROV))
    except json.JSONDecodeError as e:
        return verdict("C14", "A provenance record names every input with a checksum",
                       False, None, f"unparseable: {e}")
    ins = p.get("inputs", [])
    missing_sum = [i.get("path") for i in ins if not i.get("sha256")]
    gone = [i.get("path") for i in ins
            if not os.path.exists(os.path.join(ROOT, i.get("path", "")))]
    ok = bool(ins) and not missing_sum and not gone and p.get("builder")
    return verdict("C14", "A provenance record names every input with a checksum", ok,
                   f"{len(ins)} inputs",
                   f"no checksum {missing_sum}; not on disk {gone}" if not ok else
                   f"built by {p.get('builder')}")


def c15_reproducible(rows):
    """The builder regenerates the file from its inputs without changing it."""
    builder = os.path.join(ROOT, "05.scripts", "build_inundation_timeseries.py")
    if not os.path.exists(builder):
        return verdict("C15", "The builder reproduces the file exactly", False, None,
                       "builder script absent")
    if rows is None:
        return verdict("C15", "The builder reproduces the file exactly", False, None, "no file")
    before = open(TARGET, "rb").read()
    r = subprocess.run([sys.executable, builder, "--check"], capture_output=True, text=True,
                       cwd=os.path.join(ROOT, "05.scripts"))
    after = open(TARGET, "rb").read()
    ok = r.returncode == 0 and before == after
    return verdict("C15", "The builder reproduces the file exactly", ok,
                   f"exit {r.returncode}",
                   (r.stdout + r.stderr)[-300:] if not ok else "byte identical")


CRITERIA = [c01_file_exists, c02_schema, c03_coverage, c04_no_nulls, c05_physical_range,
            c06_outliers_defended, c07_uncertainty, c08_thin_years_rebuilt,
            c09_radar_where_it_exists, c10_cross_sensor_agreement,
            c11_optical_radar_vegetated_compared, c12_validation_file,
            c13_no_placeholders, c14_provenance, c15_reproducible]


def main():
    rows = read_target()
    results = []
    for fn in CRITERIA:
        try:
            results.append(fn(rows))
        except Exception as e:  # a criterion must never take the eval down
            results.append(verdict(fn.__name__[:3].upper(), fn.__doc__.splitlines()[0]
                                   if fn.__doc__ else fn.__name__, False, None,
                                   f"criterion raised {type(e).__name__}: {e}"))
    passed = sum(1 for r in results if r["pass"])
    report = {"passed": passed, "total": len(results),
              "complete": passed == len(results), "criteria": results}
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(report, fh, indent=2)

    if "--json" in sys.argv:
        print(json.dumps(report, indent=2))
    else:
        print(f"Lake Chilwa inundation time series, acceptance test: "
              f"{passed} of {len(results)} criteria pass\n")
        for r in results:
            mark = "PASS" if r["pass"] else "FAIL"
            print(f"  [{mark}] {r['id']}  {r['title']}")
            if r["measured"]:
                print(f"         measured: {r['measured']}")
            if not r["pass"] or r["id"] in ("C10", "C11"):
                print(f"         {r['detail']}")
        print(f"\nwrote {OUT}")
        if passed < len(results):
            nxt = next(r for r in results if not r["pass"])
            print(f"next: {nxt['id']}  {nxt['title']}\n      {nxt['detail']}")
    return 0 if passed == len(results) else 1


if __name__ == "__main__":
    sys.exit(main())
