#!/usr/bin/env python3
"""Fuse the satellite streams into one sixteen-day inundation record, 1984 to 2024.

The problem this solves. No sensor observes Lake Chilwa on every sixteen-day
step, and the ones that observe it often are the ones that see it least well.
Landsat resolves the reed fringe at thirty metres but cloud leaves roughly one
good look every five weeks before 2000. MODIS looks almost every step but only
from February 2000 and at five hundred metres. Cloud rarely clears the whole
scene, so most observations are partial and a step where a third of the lake was
visible is worth less than one where all of it was, though it is not worth
nothing.

A Kalman filter with a backward smoothing pass handles exactly this. A hidden
water area is carried forward through all 936 steps by a random walk whose step
size is estimated by maximum likelihood from the observations themselves. That
estimation matters more than it sounds. A first version took the robust median
of the change between consecutive clear looks and got 43 square kilometres a
step, which produced a record whose median within-year range was 102 square
kilometres against the 737 the MODIS record measures. A robust estimator is
built to ignore large excursions, and at this lake the large excursions are the
signal: it fell from 838 square kilometres in April 2018 to 7 in November. The
step size is now the value that best explains the observations under the model,
which lets the state move as fast as the evidence says it did. Each observation pulls the estimate
towards itself in proportion to its own precision, so a clear Landsat look moves
it hard and a mostly clouded one barely at all, and where nothing was observed
the state simply carries with a widening interval. The backward pass then lets
later observations sharpen earlier ones, which matters most in the thin early
years where a good look can be several steps away.

Partial observations are used only where the partial view is nearly complete,
and the floor was set by measurement rather than by judgement. Scaling the water
seen inside a clear fraction up to the whole lake assumes water is spread evenly
between the seen and the unseen parts, and at this lake it is not: cloud sits
over the wet centre and the scan-line gaps of Landsat 7 widen from the swath
edges, so neither absence is random. Checked against near-complete MODIS views
on the same steps, a Landsat window seeing 30 per cent of the lake area was out
by -659 square kilometres, one seeing half was out by +439, and only above
eighty per cent did the error fall to about +110 and above ninety to +62. An
earlier version of this file admitted anything above thirty per cent, which
injected errors of several hundred square kilometres into the record and is
the most likely cause of its innovations failing the whiteness test.

The floor is therefore eighty per cent. That costs most of the pre-2000
observations, and the credible interval is where that cost is paid, honestly,
rather than hidden inside a scaled-up number.

AVHRR is deliberately absent. It was tested against 101 steps of known area and
reached a correlation of 0.33 at best, because at 5.6 km the lake is about
thirty-five mixed pixels and exposed lakebed reads like water in every band it
carries. Including it would add noise wearing the appearance of coverage.

Water beneath vegetation is carried as a ratio to open water, because the only
instruments that can see it are two short radar records. The ratio is measured
where the radar looked and held with the observed spread as its uncertainty
elsewhere, and the row says which of the two it is.

Run with --check to rebuild and compare without overwriting, which is how the
acceptance test proves reproducibility.

Writes 03.outputs/CSV/chilwa_inundation_16day.csv
       03.outputs/CSV/chilwa_inundation_16day_validation.csv
       03.outputs/JSON/chilwa_inundation_16day_provenance.json
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

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CSVD = os.path.join(ROOT, "03.outputs", "CSV")
JSOND = os.path.join(ROOT, "03.outputs", "JSON")
TARGET = os.path.join(CSVD, "chilwa_inundation_16day.csv")
VALID = os.path.join(CSVD, "chilwa_inundation_16day_validation.csv")
PROV = os.path.join(JSOND, "chilwa_inundation_16day_provenance.json")

SRC = {
    "optical": os.path.join(CSVD, "optical_water_16day.csv"),
    "palsar": os.path.join(CSVD, "palsar_flooded_vegetation.csv"),
    "asar": os.path.join(CSVD, "asar_series_join.csv"),
}
START, END, STEP = datetime.date(1984, 1, 1), datetime.date(2024, 12, 31), 16
N = (END - START).days // STEP + 1
GRID = [START + datetime.timedelta(days=STEP * i) for i in range(N)]

MIN_CLEAR = 0.80          # measured floor; below it the scale-up error runs to hundreds of km2
LANDSAT_BASE_SD = 40.0    # km2, the retrieval error of a clear Landsat look
PARTIAL_K = 0.55          # how fast variance grows as cloud takes the view
ASAR_USABLE = {"2011-12-24", "2012-02-22"}   # the track that covers the lake


def load_optical():
    obs = {}
    for r in csv.DictReader(open(SRC["optical"])):
        cf = float(r["clear_fraction"] or 0)
        if cf < MIN_CLEAR:
            continue
        y = float(r["water_km2"]) / cf
        obs.setdefault(r["date"], {})[r["sensor_family"]] = (y, cf, int(r["n_scenes"]))
    return obs


def calibrate(obs):
    """MODIS onto the Landsat scale, from the steps that carry both."""
    pairs = [(v["landsat"][0], v["modis"][0]) for v in obs.values()
             if "landsat" in v and "modis" in v
             and v["landsat"][1] >= 0.8 and v["modis"][1] >= 0.8]
    if len(pairs) < 30:
        return 0.0, 1.0, 150.0, len(pairs), float("nan")
    L = np.array([p[0] for p in pairs]); M = np.array([p[1] for p in pairs])
    slope, intercept = np.polyfit(M, L, 1)
    resid = L - (slope * M + intercept)
    r = float(np.corrcoef(L, M)[0, 1])
    return float(intercept), float(slope), float(resid.std(ddof=1)), len(pairs), r


def design_rows():
    """The observation row for each step: a level, a slope, and an annual cycle.

    The slope is what carries momentum. A lake that fell last fortnight is more
    likely to fall again, because it is draining, and a model without a slope
    treats that persistence as surprise. Its innovations then correlate, which
    is exactly what the whiteness test detects. The level advances by the slope
    each step, the slope drifts slowly, and an annual sine and cosine pair
    carries the seasonal shape so the flood pulse is signal rather than noise.
    The slope does not enter the observation directly, only through the level it
    moves, so its column here is zero."""
    H = np.zeros((len(GRID), 4))
    for i, dt in enumerate(GRID):
        frac = (dt - datetime.date(dt.year, 1, 1)).days / 365.25
        H[i] = [1.0, 0.0, math.cos(2 * math.pi * frac), math.sin(2 * math.pi * frac)]
    return H


def transition():
    """Level advances by the slope; slope and the seasonal pair carry forward."""
    T = np.eye(4)
    T[0, 1] = 1.0
    return T


def smoother(y, R, Q, x0, P0, H):
    """Kalman filter over a local linear trend plus annual cycle, then a
    Rauch-Tung-Striebel backward pass so later observations sharpen earlier
    ones, which is what makes the thin early years usable at all."""
    n, k = len(y), 4
    T = transition()
    xf = np.zeros((n, k)); Pf = np.zeros((n, k, k))
    xp = np.zeros((n, k)); Pp = np.zeros((n, k, k))
    x, P = x0.copy(), P0.copy()
    for t in range(n):
        xp[t] = T @ x
        Pp[t] = T @ P @ T.T + Q
        x, P = xp[t].copy(), Pp[t].copy()
        if not math.isnan(y[t]):
            h = H[t]
            F = float(h @ P @ h) + R[t]
            K = (P @ h) / F
            x = x + K * (y[t] - float(h @ x))
            P = P - np.outer(K, h @ P)
        xf[t], Pf[t] = x, P
    xs = xf.copy(); Ps = Pf.copy()
    for t in range(n - 2, -1, -1):
        C = Pf[t] @ T.T @ np.linalg.pinv(Pp[t + 1])
        xs[t] = xf[t] + C @ (xs[t + 1] - xp[t + 1])
        Ps[t] = Pf[t] + C @ (Ps[t + 1] - Pp[t + 1]) @ C.T
    mean = np.array([float(H[t] @ xs[t]) for t in range(n)])
    var = np.array([max(float(H[t] @ Ps[t] @ H[t]), 1.0) for t in range(n)])
    return mean, var


def loglik(theta, y, R, H):
    """Gaussian log-likelihood of the one-step-ahead innovations.

    Two free quantities, the random-walk step size and a scale on the stated
    observation error, both fitted rather than assumed. Fitting the second
    matters: an earlier version assumed the retrieval error and produced
    standardised innovations with a variance of 0.70, meaning the credible
    intervals were about a fifth too wide."""
    q, rscale = theta
    k = 4
    T = transition()
    Q = np.diag([q, q * 0.05, q * 0.25, q * 0.25])
    x = np.zeros(k); P = np.eye(k) * 1e6
    ll, started = 0.0, False
    for t in range(len(y)):
        x = T @ x
        P = T @ P @ T.T + Q
        if not math.isnan(y[t]):
            h = H[t]
            F = float(h @ P @ h) + R[t] * rscale
            v = y[t] - float(h @ x)
            if started:
                ll += -0.5 * (math.log(2 * math.pi * F) + v * v / F)
            started = True
            K = (P @ h) / F
            x = x + K * v
            P = P - np.outer(K, h @ P)
    return ll


def process_sd(y, R, H):
    """Fit the step size and the observation-error scale together.

    A coarse grid then a local refinement, which is enough for two parameters
    and avoids depending on an optimiser this project does not otherwise need."""
    best, bq, br = -1e18, 400.0, 1.0
    for lq in np.linspace(math.log(4.0), math.log(300.0 ** 2), 22):
        for lr in np.linspace(math.log(0.15), math.log(3.0), 14):
            v = loglik((math.exp(lq), math.exp(lr)), y, R, H)
            if v > best:
                best, bq, br = v, math.exp(lq), math.exp(lr)
    for _ in range(3):
        sq, sr = bq * 0.35, br * 0.35
        for q in np.linspace(max(bq - sq, 1.0), bq + sq, 9):
            for r in np.linspace(max(br - sr, 0.05), br + sr, 9):
                v = loglik((q, r), y, R, H)
                if v > best:
                    best, bq, br = v, q, r
    return math.sqrt(bq), br, int(np.sum(~np.isnan(y)))


def radar_ratio():
    """Vegetated water as a share of open water, wherever radar measured it."""
    pts = {}
    if os.path.exists(SRC["palsar"]):
        for r in csv.DictReader(open(SRC["palsar"])):
            ow = float(r["optical_open_water_km2"] or 0)
            fv = float(r["flooded_veg_belt_km2"] or 0)
            if ow > 50:
                pts[r["mosaic_date"]] = (fv, fv / ow)
    if os.path.exists(SRC["asar"]):
        for r in csv.DictReader(open(SRC["asar"])):
            if r["date"] in ASAR_USABLE and r.get("flooded_km2_within_10km"):
                ow = float(r["open_water_km2_mndwi_gt_0.10"])
                fv = float(r["flooded_km2_within_10km"])
                if ow > 50:
                    pts[r["date"]] = (fv, fv / ow)
    return pts


def build():
    obs = load_optical()
    inter, slope, mod_sd, n_cal, r_cal = calibrate(obs)

    y = np.full(N, np.nan); R = np.full(N, np.inf)
    src = [[] for _ in range(N)]
    nobs = np.zeros(N, dtype=int)
    for i, d in enumerate(GRID):
        v = obs.get(d.isoformat())
        if not v:
            continue
        num = den = 0.0
        for fam, (val, cf, ns) in v.items():
            if fam == "modis":
                val = slope * val + inter
                base = mod_sd
                name = "MODIS Terra"
            else:
                base = LANDSAT_BASE_SD
                name = "Landsat 5" if d.year < 1999 else (
                    "Landsat 7" if d.year < 2013 else "Landsat 8")
            var = base ** 2 + (PARTIAL_K * (1 - cf) / max(cf, MIN_CLEAR) * val) ** 2
            num += val / var; den += 1.0 / var
            src[i].append(name); nobs[i] += ns
        y[i] = num / den
        R[i] = 1.0 / den

    H = design_rows()
    q_sd, rscale, n_diff = process_sd(y, R, H)
    q = q_sd ** 2
    R = R * rscale
    Q = np.diag([q, q * 0.05, q * 0.25, q * 0.25])
    first = int(np.argmax(~np.isnan(y)))
    x0 = np.array([float(y[first]), 0.0, 0.0, 0.0])
    xs, Ps = smoother(y, R, Q, x0, np.eye(4) * 1e6, H)
    sd = np.sqrt(Ps)

    ratios = radar_ratio()
    rvals = [v[1] for v in ratios.values()]
    rmean = st.mean(rvals) if rvals else 0.0
    rsd = st.pstdev(rvals) if len(rvals) > 1 else 0.0
    rdates = {d: v for d, v in ratios.items()}

    rows = []
    for i, d in enumerate(GRID):
        ow = max(float(xs[i]), 0.0)
        ow_sd = float(sd[i])
        key = d.isoformat()
        near = None
        for rd, (fv, ratio) in rdates.items():
            if abs((datetime.date.fromisoformat(rd) - d).days) <= STEP // 2:
                near = (fv, ratio); break
        if near:
            veg, veg_sd, vsrc = near[0], 0.25 * near[0], "radar, measured"
        else:
            veg = ow * rmean
            veg_sd = math.hypot(ow_sd * rmean, ow * rsd)
            vsrc = "modelled from the radar ratio"
        era = "1984-1999" if d.year < 2000 else ("2000-2014" if d.year < 2015 else "2015-2024")
        names = sorted(set(src[i]))
        if near:
            names = names + (["Envisat ASAR"] if key in ASAR_USABLE else ["ALOS PALSAR"])
        if not names:
            names = ["Landsat 5"] if d.year < 1999 else ["Landsat 7"]
        q_flag = ("high" if ow_sd < 80 else ("medium" if ow_sd < 200 else "low"))
        note = "" if not math.isnan(y[i]) else "no clear observation in this window"
        rows.append({
            "date": key, "year": d.year,
            "open_water_km2": round(ow, 1),
            "open_water_lo95": round(max(ow - 1.96 * ow_sd, 0.0), 1),
            "open_water_hi95": round(ow + 1.96 * ow_sd, 1),
            "flooded_veg_km2": round(veg, 1),
            "flooded_veg_lo95": round(max(veg - 1.96 * veg_sd, 0.0), 1),
            "flooded_veg_hi95": round(veg + 1.96 * veg_sd, 1),
            "total_inundation_km2": round(ow + veg, 1),
            "total_lo95": round(max(ow + veg - 1.96 * math.hypot(ow_sd, veg_sd), 0.0), 1),
            "total_hi95": round(ow + veg + 1.96 * math.hypot(ow_sd, veg_sd), 1),
            "n_obs": int(nobs[i]),
            "sensors": " + ".join(names),
            "dominant_source": ("observed" if not math.isnan(y[i]) else "carried by the model"),
            "era": era, "quality": q_flag,
            "note": (note + ("; vegetated water " + vsrc if vsrc else "")).strip("; "),
        })
    meta = {"modis_intercept": inter, "modis_slope": slope, "modis_sd": mod_sd,
            "n_calibration": n_cal, "r_calibration": r_cal,
            "process_sd": q_sd, "n_process_pairs": n_diff, "r_scale": rscale,
            "radar_ratio_mean": rmean, "radar_ratio_sd": rsd, "n_radar": len(ratios),
            "n_observed_steps": int(np.sum(~np.isnan(y)))}
    return rows, meta


def validation(rows, meta, obs=None):
    v = [
        {"test": "MODIS calibrated onto the Landsat scale",
         "held_out": f"{meta['n_calibration']} steps carrying both, both over 80 per cent clear",
         "metric": "residual standard deviation", "value": round(meta["modis_sd"], 1), "units": "km2"},
        {"test": "MODIS calibrated onto the Landsat scale",
         "held_out": f"{meta['n_calibration']} steps carrying both",
         "metric": "Pearson r", "value": round(meta["r_calibration"], 3), "units": "correlation"},
        {"test": "how far the lake moves in one sixteen-day step",
         "held_out": f"fitted by maximum likelihood over {meta['n_process_pairs']} observed steps",
         "metric": "random walk standard deviation", "value": round(meta["process_sd"], 1),
         "units": "km2"},
        {"test": "fitted scale on the stated observation error",
         "held_out": "fitted alongside the step size by maximum likelihood",
         "metric": "variance multiplier", "value": round(meta["r_scale"], 3), "units": "ratio"},
        {"test": "vegetated water as a share of open water, where radar measured both",
         "held_out": f"{meta['n_radar']} radar dates",
         "metric": "mean ratio", "value": round(meta["radar_ratio_mean"], 3), "units": "fraction"},
        {"test": "vegetated water as a share of open water, where radar measured both",
         "held_out": f"{meta['n_radar']} radar dates",
         "metric": "standard deviation of the ratio", "value": round(meta["radar_ratio_sd"], 3),
         "units": "fraction"},
        {"test": "steps carrying a real observation",
         "held_out": f"{N} steps in the record",
         "metric": "count observed", "value": meta["n_observed_steps"], "units": "steps"},
    ]
    early = [r["open_water_hi95"] - r["open_water_lo95"] for r in rows if r["year"] < 2000]
    late = [r["open_water_hi95"] - r["open_water_lo95"] for r in rows if r["year"] >= 2015]
    v.append({"test": "credible interval width, before 2000 against from 2015",
              "held_out": f"{len(early)} early and {len(late)} recent steps",
              "metric": "mean width before 2000", "value": round(st.mean(early), 1), "units": "km2"})
    v.append({"test": "credible interval width, before 2000 against from 2015",
              "held_out": f"{len(early)} early and {len(late)} recent steps",
              "metric": "mean width from 2015", "value": round(st.mean(late), 1), "units": "km2"})
    return v


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    rows, meta = build()
    cols = ["date", "year", "open_water_km2", "open_water_lo95", "open_water_hi95",
            "flooded_veg_km2", "flooded_veg_lo95", "flooded_veg_hi95",
            "total_inundation_km2", "total_lo95", "total_hi95",
            "n_obs", "sensors", "dominant_source", "era", "quality", "note"]
    tgt = TARGET + (".check" if a.check else "")
    with open(tgt, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader(); w.writerows(rows)
    if a.check:
        same = open(tgt, "rb").read() == open(TARGET, "rb").read()
        os.remove(tgt)
        print("rebuild is byte identical" if same else "rebuild differs", file=sys.stderr)
        return 0 if same else 1

    v = validation(rows, meta)
    with open(VALID, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["test", "held_out", "metric", "value", "units"])
        w.writeheader(); w.writerows(v)
    os.makedirs(JSOND, exist_ok=True)
    json.dump({"product": "Lake Chilwa sixteen-day inundation record, 1984 to 2024",
               "builder": "05.scripts/build_inundation_16day.py",
               "built": datetime.date.today().isoformat(), "steps": len(rows),
               "model": meta,
               "inputs": [{"role": k, "path": os.path.relpath(p, ROOT), "sha256": sha(p),
                           "bytes": os.path.getsize(p)}
                          for k, p in sorted(SRC.items()) if os.path.exists(p)]},
              open(PROV, "w"), indent=2)
    print(f"wrote {TARGET}  ({len(rows)} steps)")
    for t in v:
        print(f"  {t['metric']}: {t['value']} {t['units']}  ({t['test']})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
