#!/usr/bin/env python3
"""Is the fused record internally consistent, or does it carry the seams of the
sensors that built it?

A record assembled from five instruments across forty-one years can look
perfectly smooth and still be wrong in the way that matters most: a step where
one sensor hands over to another, a period where one instrument quietly pulls
the estimate up, an uncertainty that is stated but not earned. None of that
shows in a plot. These are the tests that find it.

Four questions, each with a test that can fail.

1. Is the model correctly specified? In a properly specified Kalman filter the
   one-step-ahead innovations, the surprise in each observation, are white noise
   with unit variance once standardised by the model's own predicted variance.
   Autocorrelation in them means the state cannot follow what the data are doing.
   A standardised variance far from one means the stated credible intervals are
   the wrong size, too confident below one and too timid above it.

2. Does any sensor pull the record? The innovations are grouped by which
   instrument produced the observation and tested for a non-zero mean. A sensor
   whose innovations are systematically positive is dragging the estimate up
   whenever it speaks, which is a seam in the temporal axis even where no step
   change is visible.

3. Is there a discontinuity where the instruments change hands? The level either
   side of each handover is compared, on windows short enough that real
   hydrology cannot explain a jump. Landsat 7 arrives in 1999, MODIS in 2000,
   Landsat 8 in 2013.

4. Does the record depend on any one sensor? The whole thing is rebuilt with
   each instrument withheld in turn. If dropping one moves the series materially,
   the fusion is not a consensus but a single sensor with company.

Writes 03.outputs/CSV/cube_consistency_tests.csv
"""

import csv
import datetime
import math
import os
import statistics as st

import numpy as np

import build_inundation_16day as B

OUT = os.path.join(B.CSVD, "cube_consistency_tests.csv")


def assemble(families=("landsat", "modis")):
    """Observation vector and variance, restricted to the named instruments."""
    obs = B.load_optical()
    inter, slope, mod_sd, n_cal, r_cal = B.calibrate(obs)
    y = np.full(B.N, np.nan); R = np.full(B.N, np.inf)
    who = [None] * B.N
    for i, d in enumerate(B.GRID):
        v = obs.get(d.isoformat())
        if not v:
            continue
        num = den = 0.0; names = []
        for fam, (val, cf, ns) in v.items():
            if fam not in families:
                continue
            if fam == "modis":
                val = slope * val + inter; base = mod_sd
            else:
                base = B.LANDSAT_BASE_SD
            var = base ** 2 + (B.PARTIAL_K * (1 - cf) / max(cf, B.MIN_CLEAR) * val) ** 2
            num += val / var; den += 1.0 / var; names.append(fam)
        if den > 0:
            y[i] = num / den; R[i] = 1.0 / den; who[i] = "+".join(sorted(names))
    return y, R, who


def innovations(y, R, q, H):
    """One-step-ahead surprise and its predicted variance, standardised."""
    k = 4
    T = B.transition()
    Q = np.diag([q, q * 0.05, q * 0.25, q * 0.25])
    x = np.zeros(k); P = np.eye(k) * 1e6
    out = []
    started = False
    for t in range(len(y)):
        x = T @ x
        P = T @ P @ T.T + Q
        if not math.isnan(y[t]):
            h = H[t]
            F = float(h @ P @ h) + R[t]
            v = y[t] - float(h @ x)
            if started:
                out.append((t, v, F, v / math.sqrt(F)))
            started = True
            K = (P @ h) / F
            x = x + K * v
            P = P - np.outer(K, h @ P)
    return out


def ljung_box(z, lags=20):
    """Whiteness of the standardised innovations. Large Q means structure left."""
    z = np.asarray(z, dtype=float); n = len(z)
    z = z - z.mean()
    denom = float(np.sum(z * z)) or 1e-12
    Q = 0.0
    for k in range(1, lags + 1):
        rk = float(np.sum(z[k:] * z[:-k])) / denom
        Q += rk * rk / (n - k)
    Q *= n * (n + 2)
    # chi-square survival with `lags` degrees of freedom, Wilson-Hilferty
    d = lags
    x = Q / d
    zscore = (x ** (1 / 3) - (1 - 2 / (9 * d))) / math.sqrt(2 / (9 * d))
    p = 0.5 * math.erfc(zscore / math.sqrt(2))
    return Q, p, lags


def welch(a, b):
    """Two-sample test with unequal variances, returning t and a two-sided p."""
    if len(a) < 5 or len(b) < 5:
        return float("nan"), float("nan"), float("nan")
    ma, mb = st.mean(a), st.mean(b)
    va, vb = st.variance(a), st.variance(b)
    se = math.sqrt(va / len(a) + vb / len(b)) or 1e-12
    t = (ma - mb) / se
    df = (va / len(a) + vb / len(b)) ** 2 / (
        (va / len(a)) ** 2 / (len(a) - 1) + (vb / len(b)) ** 2 / (len(b) - 1))
    p = 2 * (1 - _t_cdf(abs(t), df))
    return t, p, ma - mb


def _t_cdf(t, df):
    x = df / (df + t * t)
    return 1 - 0.5 * _betainc(df / 2, 0.5, x)


def _betainc(a, b, x):
    """Regularised incomplete beta, continued fraction, enough for a p-value."""
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    lbeta = math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)
    front = math.exp(math.log(x) * a + math.log(1 - x) * b - lbeta) / a
    f, c, d = 1.0, 1.0, 0.0
    for i in range(200):
        m = i // 2
        if i == 0:
            num = 1.0
        elif i % 2 == 0:
            num = (m * (b - m) * x) / ((a + 2 * m - 1) * (a + 2 * m))
        else:
            num = -((a + m) * (a + b + m) * x) / ((a + 2 * m) * (a + 2 * m + 1))
        d = 1.0 + num * d
        if abs(d) < 1e-30:
            d = 1e-30
        d = 1.0 / d
        c = 1.0 + num / c
        if abs(c) < 1e-30:
            c = 1e-30
        f *= c * d
        if abs(1 - c * d) < 1e-10:
            break
    return front * (f - 1)


def one_sided_p_mean_zero(v):
    """Two-sided p for the mean of standardised innovations being zero."""
    n = len(v)
    if n < 5:
        return float("nan"), float("nan")
    m, s = st.mean(v), st.stdev(v)
    t = m / (s / math.sqrt(n) or 1e-12)
    return t, 2 * (1 - _t_cdf(abs(t), n - 1))


def main():
    rows = []
    y, R, who = assemble()
    H = B.design_rows()
    q_sd, rscale, _ = B.process_sd(y, R, H)
    R = R * rscale
    inn = innovations(y, R, q_sd ** 2, H)
    z = [a[3] for a in inn]

    print(f"{len(inn)} innovations from {int(np.sum(~np.isnan(y)))} observed steps\n")
    print("1. Is the model correctly specified?")
    Q, p, lags = ljung_box(z)
    print(f"   Ljung-Box on the standardised innovations, {lags} lags: "
          f"Q = {Q:,.1f}, p = {p:.4f}"
          f"   {'WHITE, no structure left' if p > 0.05 else 'STRUCTURE REMAINS'}")
    v = st.pvariance(z)
    print(f"   variance of the standardised innovations = {v:.3f} "
          f"(one means the stated intervals are the right size; "
          f"{'too wide' if v < 0.8 else ('too narrow' if v > 1.25 else 'well calibrated')})")
    tm, pm = one_sided_p_mean_zero(z)
    print(f"   mean = {st.mean(z):+.3f}, t = {tm:+.2f}, p = {pm:.3f}"
          f"   {'no overall bias' if pm > 0.05 else 'BIASED'}")
    rows += [
        {"axis": "temporal", "test": "Ljung-Box on standardised innovations",
         "statistic": "Q", "value": round(Q, 1), "p": round(p, 4),
         "verdict": "white" if p > 0.05 else "structure remains"},
        {"axis": "temporal", "test": "calibration of the stated uncertainty",
         "statistic": "variance of standardised innovations", "value": round(v, 3),
         "p": "", "verdict": "well calibrated" if 0.8 <= v <= 1.25 else "mis-sized"},
        {"axis": "temporal", "test": "overall bias of the fused record",
         "statistic": "mean standardised innovation", "value": round(st.mean(z), 3),
         "p": round(pm, 4), "verdict": "unbiased" if pm > 0.05 else "biased"},
    ]

    print("\n2. Does any single instrument pull the record?")
    by = {}
    for t, vv, F, zz in inn:
        by.setdefault(who[t] or "?", []).append(zz)
    for fam in sorted(by):
        g = by[fam]
        tt, pp = one_sided_p_mean_zero(g)
        flag = "no pull" if (math.isnan(pp) or pp > 0.05) else "PULLS THE RECORD"
        print(f"   {fam:18s} n = {len(g):4d}   mean {st.mean(g):+.3f}   "
              f"t = {tt:+.2f}   p = {pp:.3f}   {flag}")
        rows.append({"axis": "sensor", "test": f"innovation bias, {fam}",
                     "statistic": "mean standardised innovation",
                     "value": round(st.mean(g), 3), "p": round(pp, 4) if pp == pp else "",
                     "verdict": flag.lower()})

    print("\n3. Is there a step where the instruments change hands?")
    for name, year in (("Landsat 7 arrives", 1999), ("MODIS arrives", 2000),
                       ("Landsat 8 arrives", 2013)):
        before = [a[3] for a in inn
                  if B.GRID[a[0]].year in (year - 2, year - 1)]
        after = [a[3] for a in inn if B.GRID[a[0]].year in (year, year + 1)]
        tt, pp, diff = welch(after, before)
        flag = ("no step" if (math.isnan(pp) or pp > 0.05) else "STEP DETECTED")
        print(f"   {name:20s} {year}   n = {len(before)}/{len(after)}   "
              f"shift {diff:+.3f} sd   t = {tt:+.2f}   p = {pp:.3f}   {flag}")
        rows.append({"axis": "temporal", "test": f"handover step, {name} {year}",
                     "statistic": "shift in standardised innovation",
                     "value": round(diff, 3) if diff == diff else "",
                     "p": round(pp, 4) if pp == pp else "", "verdict": flag.lower()})

    print("\n4. Does the record depend on any one instrument?")
    full, _ = _run(y, R, H, q_sd)
    for drop in ("modis", "landsat"):
        keep = tuple(f for f in ("landsat", "modis") if f != drop)
        y2, R2, _ = assemble(keep)
        n2 = int(np.sum(~np.isnan(y2)))
        q2, r2, _ = B.process_sd(y2, R2, H)
        alt, _ = _run(y2, R2 * r2, H, q2)
        both = ~np.isnan(full) & ~np.isnan(alt)
        d = full[both] - alt[both]
        print(f"   without {drop:9s} {n2:3d} steps remain   "
              f"mean shift {np.mean(d):+7.1f} km2   "
              f"largest {np.max(np.abs(d)):7.1f} km2   "
              f"correlation {np.corrcoef(full[both], alt[both])[0,1]:.4f}")
        rows.append({"axis": "sensor", "test": f"record rebuilt without {drop}",
                     "statistic": "mean shift", "value": round(float(np.mean(d)), 1),
                     "p": "", "verdict": f"r = {np.corrcoef(full[both], alt[both])[0,1]:.4f}, "
                                         f"largest {np.max(np.abs(d)):.0f} km2"})

    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["axis", "test", "statistic", "value", "p", "verdict"])
        w.writeheader(); w.writerows(rows)
    print(f"\nwrote {OUT}")


def _run(y, R, H, q_sd):
    q = q_sd ** 2
    first = int(np.argmax(~np.isnan(y)))
    x0 = np.array([float(y[first]), 0.0, 0.0, 0.0])
    return B.smoother(y, R, np.diag([q, q * .05, q * .25, q * .25]),
                      x0, np.eye(4) * 1e6, H)


if __name__ == "__main__":
    main()
