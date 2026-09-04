#!/usr/bin/env python3
"""The acceptance test for the Lake Chilwa inundation time series.

This file is the contract. The dataset is finished when every criterion passes
and not before, and nothing else decides it: not a plausible-looking plot, not
a script that ran without error, not a number that resembles a published one.

The target is a sixteen-day record of inundated area over the basin from 1984
to 2024, carrying open water and water standing within vegetation separately,
each with a credible interval. Sixteen days is the Landsat repeat, and the
record has to run at that step across the whole span because the study is about
two things at once: a recession and refilling cycle that turns over across
roughly fifteen years, and a shoreline that moves kilometres within a season.
An annual series answers neither, and measuring how far the lake moves inside a
year turned out to need more care than it looks. Three figures were tried. A
committed MODIS eight-day file gives a median within-year range of 737 square
kilometres, and that number is inflated, almost certainly by cloud passing as
water. Direct observation on 25 years whose clear looks span both seasons gives
291. But a maximum minus a minimum is the statistic most exposed to noise, and a
perfectly still lake observed 31 times with the measured 59 square kilometre
retrieval error would show a spread of about 227 square kilometres from noise
alone, so most of that 291 is not movement either.

The criterion therefore compares a robust spread, the tenth to ninetieth
percentile, between the finished record and its own verified observations on the
same years. That is like for like, it is self-calibrating, and it cannot be
moved by a bad file elsewhere. Setting the earlier threshold from an
unverified external number was a mistake and this replaces it.

No single sensor spans the record. Landsat carries the detail from 1984 but is
cloud-limited to a handful of clear looks a year in the wet season; AVHRR is
daily from 1981 at 5 km; MODIS is daily from 2000 at 500 m; Sentinel-1, 2 and 3
and PROBA-V densify from 2014. The record is therefore a fusion, and the honest
consequence is that precision degrades backwards. Criterion C07 requires the
intervals to say so rather than reporting the same confidence in 1987 as in
2020.

Run: python3 eval_timeseries.py [--json]
Exit status is 0 when every criterion passes, so a loop can branch on it.
"""

import csv
import datetime
import json
import math
import os
import statistics as st
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CSVD = os.path.join(ROOT, "03.outputs", "CSV")
TARGET = os.path.join(CSVD, "chilwa_inundation_16day.csv")
VALID = os.path.join(CSVD, "chilwa_inundation_16day_validation.csv")
PROV = os.path.join(ROOT, "03.outputs", "JSON", "chilwa_inundation_16day_provenance.json")
OUT = os.path.join(ROOT, "03.outputs", "JSON", "eval_timeseries.json")
BUILDER = os.path.join(ROOT, "05.scripts", "build_inundation_16day.py")

START = datetime.date(1984, 1, 1)
END = datetime.date(2024, 12, 31)
STEP = 16
N_STEPS = (END - START).days // STEP + 1

BASIN_KM2 = 8752.0
# Lake Chilwa genuinely empties. Direct Landsat observation on 14 November 2018,
# six scenes with 99 per cent of the basin clear, gives 7 km2 of open water, and
# the MODIS record puts that year's seasonal minimum at 3 km2. A floor above zero
# would reject the most important event in the record, so there is none.
OPEN_MIN, OPEN_MAX = 0.0, 2600.0
TOTAL_MAX = 4500.0

# The record must retain most of the seasonal movement its own observations show.
# Some loss is correct: a maximum minus a minimum over raw observations is moved
# by a single noisy step, and removing those is what a smoother is for.
OBSERVED = os.path.join(CSVD, "optical_water_16day.csv")
MIN_RANGE_RETAINED = 0.60
PLO, PHI = 0.10, 0.90        # the robust spread compared on both sides

# first light of each sensor, so a step cannot claim one that was not in orbit
FIRST_LIGHT = {
    "AVHRR": 1981, "Landsat 4": 1982, "Landsat 5": 1984, "Landsat 7": 1999,
    "MODIS Terra": 2000, "MODIS Aqua": 2002, "Envisat ASAR": 2002,
    "ALOS PALSAR": 2006, "Landsat 8": 2013, "PROBA-V": 2013,
    "Sentinel-1": 2014, "Sentinel-2": 2015, "Sentinel-3": 2016, "Landsat 9": 2021,
}

# Criterion C16 was added on 3 September 2026 after the record had already
# passed fifteen of fifteen. Nothing in the test suite had asked whether the
# fusion model was correctly specified, so a record could satisfy every
# criterion while its state equation was unable to follow what the data did.
# The one-step-ahead innovations of a correct filter are white with unit
# variance once standardised; autocorrelation in them means the model lags the
# lake, and a standardised variance away from one means the credible intervals
# are the wrong size.
CONSISTENCY = os.path.join(CSVD, "cube_consistency_tests.csv")

REQUIRED_COLUMNS = [
    "date", "year", "open_water_km2", "open_water_lo95", "open_water_hi95",
    "flooded_veg_km2", "flooded_veg_lo95", "flooded_veg_hi95",
    "total_inundation_km2", "total_lo95", "total_hi95",
    "n_obs", "sensors", "dominant_source", "era", "quality", "note",
]
PLACEHOLDERS = ["todo", "tbd", "placeholder", "xxx", "fixme"]


def grid():
    return [START + datetime.timedelta(days=STEP * i) for i in range(N_STEPS)]


def read():
    if not os.path.exists(TARGET):
        return None
    with open(TARGET) as fh:
        return list(csv.DictReader(fh))


def num(r, k):
    v = (r or {}).get(k)
    if v is None or str(v).strip() == "":
        return None
    try:
        f = float(v)
    except ValueError:
        return None
    return None if math.isnan(f) else f


def d(r):
    try:
        return datetime.date.fromisoformat((r.get("date") or "").strip())
    except ValueError:
        return None


def V(i, t, ok, measured, detail):
    return {"id": i, "title": t, "pass": bool(ok), "measured": measured, "detail": detail}


def c01(rows):
    return V("C01", "The target file exists and parses as CSV", rows is not None,
             f"{len(rows)} rows" if rows else None, "absent" if rows is None else "")


def c02(rows):
    if rows is None:
        return V("C02", "Every required column is present", False, None, "no file")
    miss = [c for c in REQUIRED_COLUMNS if c not in rows[0]]
    return V("C02", "Every required column is present", not miss, f"{len(rows[0])} columns",
             f"missing {miss}" if miss else "")


def c03(rows):
    """A complete, regular sixteen-day grid with nothing missing or repeated."""
    if rows is None:
        return V("C03", "The sixteen-day grid is complete and regular", False, None, "no file")
    got = [d(r) for r in rows]
    if any(x is None for x in got):
        return V("C03", "The sixteen-day grid is complete and regular", False,
                 f"{len(rows)} rows", "some dates do not parse")
    want = grid()
    missing = sorted(set(want) - set(got))
    dupes = sorted({x for x in got if got.count(x) > 1})
    spacing = {(b - a).days for a, b in zip(got, got[1:])}
    ok = not missing and not dupes and spacing <= {STEP} and len(got) == N_STEPS
    return V("C03", "The sixteen-day grid is complete and regular", ok,
             f"{len(got)} of {N_STEPS} steps, {START} to {END}",
             f"missing {len(missing)} (first {missing[:3]}); duplicated {dupes[:3]}; "
             f"spacings {sorted(spacing)[:5]}" if not ok else "")


def c04(rows):
    if rows is None:
        return V("C04", "No missing value in the three area columns", False, None, "no file")
    bad = [f"{r.get('date')}:{k}" for r in rows
           for k in ("open_water_km2", "flooded_veg_km2", "total_inundation_km2")
           if num(r, k) is None]
    return V("C04", "No missing value in the three area columns", not bad,
             f"{len(bad)} empty cells", "; ".join(bad[:8]))


def c05(rows):
    """Every area physically possible, and the parts add to the whole."""
    if rows is None:
        return V("C05", "Every area is physically possible and the sum adds up", False, None, "no file")
    bad = []
    for r in rows:
        o, f, t = (num(r, "open_water_km2"), num(r, "flooded_veg_km2"),
                   num(r, "total_inundation_km2"))
        if o is not None and not (OPEN_MIN <= o <= OPEN_MAX):
            bad.append(f"{r.get('date')}: open water {o:.0f} outside {OPEN_MIN:.0f}-{OPEN_MAX:.0f}")
        if f is not None and not (0 <= f <= BASIN_KM2):
            bad.append(f"{r.get('date')}: vegetated {f:.0f} out of range")
        if t is not None and t > TOTAL_MAX:
            bad.append(f"{r.get('date')}: total {t:.0f} above {TOTAL_MAX:.0f}")
        if None not in (o, f, t) and abs(o + f - t) > 1.0:
            bad.append(f"{r.get('date')}: total does not equal the parts")
    return V("C05", "Every area is physically possible and the sum adds up", not bad,
             f"{len(bad)} violations", "; ".join(bad[:6]))


def c06(rows):
    """Every step carries a credible interval and the bounds are the right way round."""
    if rows is None:
        return V("C06", "Every step carries an ordered credible interval", False, None, "no file")
    bad = []
    for r in rows:
        for base in ("open_water", "flooded_veg", "total"):
            k = f"{base}_km2" if base != "total" else "total_inundation_km2"
            v, lo, hi = num(r, k), num(r, f"{base}_lo95"), num(r, f"{base}_hi95")
            if lo is None or hi is None:
                bad.append(f"{r.get('date')}:{base} no interval")
            elif v is not None and not (lo <= v <= hi):
                bad.append(f"{r.get('date')}:{base} {lo:.0f}/{v:.0f}/{hi:.0f} out of order")
            elif hi - lo <= 0:
                bad.append(f"{r.get('date')}:{base} zero-width interval")
    return V("C06", "Every step carries an ordered credible interval", not bad,
             f"{len(bad)} faults", "; ".join(bad[:6]))


def c07(rows):
    """The interval must be wider in the thin early record than in the dense recent one.

    Reporting the same confidence in 1987, when the only dense sensor is AVHRR
    at five kilometres, as in 2020, when Sentinel-1, 2 and 3 and two Landsats
    are all in orbit, is a false claim about how much is known."""
    if rows is None:
        return V("C07", "Uncertainty widens where the record is thin", False, None, "no file")
    early = [num(r, "open_water_hi95") - num(r, "open_water_lo95") for r in rows
             if d(r) and d(r).year <= 1999
             and num(r, "open_water_hi95") is not None and num(r, "open_water_lo95") is not None]
    late = [num(r, "open_water_hi95") - num(r, "open_water_lo95") for r in rows
            if d(r) and d(r).year >= 2015
            and num(r, "open_water_hi95") is not None and num(r, "open_water_lo95") is not None]
    if len(early) < 50 or len(late) < 50:
        return V("C07", "Uncertainty widens where the record is thin", False,
                 f"{len(early)} early, {len(late)} recent steps", "too few steps to compare")
    me, ml = st.mean(early), st.mean(late)
    return V("C07", "Uncertainty widens where the record is thin", me > ml * 1.3,
             f"mean 95 per cent width {me:.0f} km2 before 2000 against {ml:.0f} km2 from 2015",
             "the early interval must be at least 30 per cent wider")


def c09(rows):
    """The record must show the lake actually moving within the year.

    Measured against its own observations rather than an external figure, and
    on a robust spread rather than a maximum minus a minimum, because that
    difference is what separates movement from noise. The record must retain at
    least sixty per cent of the tenth-to-ninetieth-percentile spread its own
    clear observations show on the same years. Some loss is correct, since
    smoothing is what removes the retrieval error, but a record keeping well
    under that has flattened the flood pulse."""
    if rows is None:
        return V("C09", "Within-year range matches a lake that actually moves", False, None, "no file")
    if not os.path.exists(OBSERVED):
        return V("C09", "Within-year range matches a lake that actually moves", False, None,
                 "the observation file is absent, so there is nothing to compare against")
    obs = {}
    for r in csv.DictReader(open(OBSERVED)):
        try:
            cf = float(r["clear_fraction"])
        except (ValueError, KeyError):
            continue
        if cf >= 0.8:
            obs.setdefault(int(r["year"]), []).append(
                (int(r["date"][5:7]), float(r["water_km2"]) / cf))
    seasonal = {}
    for y, v in obs.items():
        wet = [a for m, a in v if m in (1, 2, 3, 4)]
        dry = [a for m, a in v if m in (8, 9, 10, 11)]
        if len(wet) >= 2 and len(dry) >= 2 and len(v) >= 10:
            seasonal[y] = [a for _, a in v]
    got = {}
    for r in rows:
        dt, v = d(r), num(r, "open_water_km2")
        if dt and v is not None:
            got.setdefault(dt.year, []).append(v)
    common = sorted(y for y in seasonal if len(got.get(y, [])) >= 15)
    if len(common) < 15:
        return V("C09", "Within-year range matches a lake that actually moves", False,
                 f"{len(common)} comparable years", "too few years whose looks span both seasons")
    def spread(v):
        z = sorted(v); n = len(z)
        return z[min(n - 1, int(PHI * n))] - z[int(PLO * n)]

    o = st.median([spread(seasonal[y]) for y in common])
    f = st.median([spread(got[y]) for y in common])
    return V("C09", "Within-year range matches a lake that actually moves",
             f >= MIN_RANGE_RETAINED * o,
             f"record keeps {f:.0f} km2 of the {o:.0f} km2 robust spread its own observations show "
             f"({100 * f / max(o, 1e-9):.0f} per cent) over {len(common)} years",
             f"must retain at least {100 * MIN_RANGE_RETAINED:.0f} per cent")


def c08(rows):
    """Every step names its sources, and none claims a sensor not yet in orbit."""
    if rows is None:
        return V("C08", "Every step names sources that existed at the time", False, None, "no file")
    bad = []
    for r in rows:
        dt = d(r)
        s = (r.get("sensors") or "").strip()
        if not s:
            bad.append(f"{r.get('date')}: no sensor named"); continue
        for name in [x.strip() for x in s.split("+") if x.strip()]:
            yr = FIRST_LIGHT.get(name)
            if yr is None:
                bad.append(f"{r.get('date')}: unknown sensor '{name}'")
            elif dt and dt.year < yr:
                bad.append(f"{r.get('date')}: claims {name}, first light {yr}")
    return V("C08", "Every step names sources that existed at the time", not bad,
             f"{len(bad)} faults", "; ".join(sorted(set(bad))[:6]))


def c10(rows):
    """Radar carries the vegetated class wherever a horizontally polarised record exists.

    Optical reflectance sees the top of a reed stand, not the water beneath it,
    which is hypothesis H2. Envisat covers the 2011 to 2012 fieldwork window and
    ALOS PALSAR the yearly mosaics from 2007, and those are the only
    horizontally polarised records over this basin: of 1,966 Sentinel-1 scenes
    exactly one transmits horizontally."""
    if rows is None:
        return V("C10", "Radar carries the vegetated class where it exists", False, None, "no file")
    want = set(range(2007, 2013)) | set(range(2015, 2025))
    seen = {dt.year for r in rows for dt in [d(r)]
            if dt and "ASAR" in (r.get("sensors") or "") or
            dt and "PALSAR" in (r.get("sensors") or "")}
    missing = sorted(want - seen)
    return V("C10", "Radar carries the vegetated class where it exists", not missing,
             f"{len(seen & want)} of {len(want)} radar-era years carry a radar step",
             f"no radar step in {missing}" if missing else "")


def c11(rows):
    """A holdout validation exists and reports its error in km2."""
    if not os.path.exists(VALID):
        return V("C11", "A holdout validation exists with error in km2", False, None, "absent")
    vr = list(csv.DictReader(open(VALID)))
    if not vr:
        return V("C11", "A holdout validation exists with error in km2", False, "0 rows", "empty")
    need = {"test", "held_out", "metric", "value", "units"}
    if not need <= set(vr[0]):
        return V("C11", "A holdout validation exists with error in km2", False, f"{len(vr)} rows",
                 f"missing {sorted(need - set(vr[0]))}")
    km2 = [r for r in vr if (r.get("units") or "").strip() == "km2"]
    return V("C11", "A holdout validation exists with error in km2", bool(km2),
             f"{len(vr)} tests, {len(km2)} in km2", "" if km2 else "no km2 error reported")


def c12(rows):
    """The record reproduces the recessions the basin is known to have had.

    1995 is the deepest recession in the instrumental record and 2012 is the
    year of the fieldwork, both drier than their surrounding decade; 2023 and
    2024 are a refill. A series that does not rank these correctly is not
    measuring this lake."""
    if rows is None:
        return V("C12", "Known recessions and refills appear in the right order", False, None, "no file")
    per = {}
    for r in rows:
        dt, v = d(r), num(r, "open_water_km2")
        if dt and v is not None:
            per.setdefault(dt.year, []).append(v)
    ann = {y: st.mean(v) for y, v in per.items() if len(v) >= 15}
    if len(ann) < 30:
        return V("C12", "Known recessions and refills appear in the right order", False,
                 f"{len(ann)} complete years", "too few years")
    med = st.median(ann.values())
    checks, bad = [], []
    for y, why in ((1995, "the deepest recession in the record"), (2012, "the fieldwork recession")):
        if y in ann:
            checks.append((y, ann[y] < med, why))
    for y, why in ((2023, "the refill"), (2024, "the refill")):
        if y in ann:
            checks.append((y, ann[y] > med, why))
    for y, ok, why in checks:
        if not ok:
            bad.append(f"{y} ({why}) at {ann[y]:.0f} km2 against a median of {med:.0f}")
    return V("C12", "Known recessions and refills appear in the right order", not bad,
             f"{len(checks)} events tested against a median of {med:.0f} km2",
             "; ".join(bad))


def c13(rows):
    if rows is None:
        return V("C13", "No placeholder text", False, None, "no file")
    bad = [f"{r.get('date')}:{k}" for r in rows for k, v in r.items()
           if any(t in (v or "").strip().lower() for t in PLACEHOLDERS)]
    return V("C13", "No placeholder text", not bad, f"{len(bad)} suspect cells",
             "; ".join(sorted(set(bad))[:6]))


def c14(rows):
    if not os.path.exists(PROV):
        return V("C14", "A provenance record names every input with a checksum", False, None, "absent")
    try:
        p = json.load(open(PROV))
    except json.JSONDecodeError as e:
        return V("C14", "A provenance record names every input with a checksum", False, None, str(e))
    ins = p.get("inputs", [])
    nosum = [i.get("path") for i in ins if not i.get("sha256")]
    gone = [i.get("path") for i in ins if not os.path.exists(os.path.join(ROOT, i.get("path", "")))]
    ok = bool(ins) and not nosum and not gone and p.get("builder")
    return V("C14", "A provenance record names every input with a checksum", ok,
             f"{len(ins)} inputs", f"no checksum {nosum}; absent {gone}" if not ok else "")


def c15(rows):
    if not os.path.exists(BUILDER):
        return V("C15", "The builder reproduces the file exactly", False, None, "builder absent")
    if rows is None:
        return V("C15", "The builder reproduces the file exactly", False, None, "no file")
    before = open(TARGET, "rb").read()
    r = subprocess.run([sys.executable, BUILDER, "--check"], capture_output=True, text=True,
                       cwd=os.path.join(ROOT, "05.scripts"))
    after = open(TARGET, "rb").read()
    ok = r.returncode == 0 and before == after
    return V("C15", "The builder reproduces the file exactly", ok, f"exit {r.returncode}",
             (r.stdout + r.stderr)[-300:] if not ok else "")


def c16(rows):
    """The fusion model must be correctly specified, not merely plausible.

    Read from the diagnostics that 05.scripts/diagnose_cube_consistency.py
    writes. Three things have to hold: the standardised innovations are white,
    their variance is near one so the stated intervals are the right size, and
    no single instrument pulls the record when it speaks alone."""
    if not os.path.exists(CONSISTENCY):
        return V("C16", "The fusion model is correctly specified", False, None,
                 "run diagnose_cube_consistency.py; its results are missing")
    t = {r["test"]: r for r in csv.DictReader(open(CONSISTENCY))}
    bad = []
    lb = t.get("Ljung-Box on standardised innovations")
    if lb and lb["verdict"] != "white":
        bad.append(f"innovations are not white (Q = {lb['value']}, p = {lb['p']})")
    cal = t.get("calibration of the stated uncertainty")
    if cal and cal["verdict"] != "well calibrated":
        bad.append(f"stated intervals are mis-sized (standardised variance {cal['value']})")
    pulls = [r for k, r in t.items() if k.startswith("innovation bias,")
             and "pulls" in r["verdict"]]
    for r in pulls:
        bad.append(f"{r['test'].split(', ')[1]} pulls the record (mean {r['value']})")
    return V("C16", "The fusion model is correctly specified", not bad,
             f"{len(t)} diagnostics read", "; ".join(bad))


CRITERIA = [c01, c02, c03, c04, c05, c06, c07, c08, c09, c10, c11, c12, c13, c14, c15, c16]


def main():
    rows = read()
    res = []
    for fn in CRITERIA:
        try:
            res.append(fn(rows))
        except Exception as e:
            res.append(V(fn.__name__.upper(), fn.__doc__.splitlines()[0] if fn.__doc__ else
                         fn.__name__, False, None, f"{type(e).__name__}: {e}"))
    res.sort(key=lambda r: r["id"])
    passed = sum(1 for r in res if r["pass"])
    rep = {"passed": passed, "total": len(res), "complete": passed == len(res), "criteria": res}
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(rep, open(OUT, "w"), indent=2)
    if "--json" in sys.argv:
        print(json.dumps(rep, indent=2))
    else:
        print(f"Lake Chilwa sixteen-day inundation record: {passed} of {len(res)} criteria pass\n")
        for r in res:
            print(f"  [{'PASS' if r['pass'] else 'FAIL'}] {r['id']}  {r['title']}")
            if r["measured"]:
                print(f"         measured: {r['measured']}")
            if not r["pass"] and r["detail"]:
                print(f"         {r['detail']}")
        if passed < len(res):
            n = next(r for r in res if not r["pass"])
            print(f"\nnext: {n['id']}  {n['title']}\n      {n['detail']}")
    return 0 if passed == len(res) else 1


if __name__ == "__main__":
    sys.exit(main())
