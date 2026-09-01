# Harmonic fit for 22 February 2012

2026-09-01. The 2012 water map is produced by a per-pixel harmonic fit, not by LandTrendr.
`05.scripts/harmonic_fit_2012.py` fits MNDWI(t) = a0 + a1*t + annual and semi-annual sine and cosine
pairs, six coefficients, on 243 Landsat 5, 7 and 8 scenes over 2010-06-01 to 2014-06-01, then
evaluates the model on 22 February 2012. Matters because LandTrendr reduces each year to one value
and joins those with straight lines, which would erase the flood pulse this study measures; a
harmonic model keeps the repeating shape of the year. See [[lemma-landtrendr-gap]].

2026-09-01. A fit inside 2012 alone was tried first and fails. The median pixel carries 8 clear
Landsat 7 observations that year against 5 free parameters, so the regression returned MNDWI of 9.96
on a scale bounded at 1 and an in-fit RMSE of 0.024, both artefacts of fitting through as many
observations as parameters. Requiring 12 observations left 20.7 per cent of the basin. Matters
because the observation floor and the widened window are load-bearing, not tuning.

2026-09-01. On the widened window every pixel in the basin clears the 12-observation floor, the
median pixel carries 48 observations, and the fitted MNDWI for 22 February 2012 runs from -0.51 at
the 5th percentile to 0.76 at the 95th, all inside the possible range. In-fit RMSE is 0.058 MNDWI.

2026-09-01. Holdout, whole dates withheld and then predicted, gives mean absolute error of 0.068
MNDWI on 2012-02-15, 0.073 on 2012-01-30 and 0.056 on 2012-04-19, medians 0.049, 0.051 and 0.029,
90th percentiles about 0.15. Matters because the holdout error is barely worse than the in-fit
residual, so the model is not overfitting, but an error of 0.05 to 0.15 flips the class of any pixel
sitting near the MNDWI = 0 water boundary, which is exactly the reed fringe this study is about.

2026-09-01. Fitted open water on 22 February 2012 is 1,158 km2 at MNDWI > 0, 1,205 km2 at > -0.10 and
1,130 km2 at > +0.10. The screened MODIS record (70 per cent clear-cell floor) gives a February to
March 2012 mean of 1,153 km2 over seven eight-day steps, a difference of 5 km2 against a MODIS
calibration RMSE of 128 km2. Matters as corroboration rather than independent validation, since both
routes read a water index from optical reflectance; the sensors, grain and dates differ, the physics
does not. Per-step MODIS values scatter from 1,050 to 1,390 km2, so the agreement of the means is
tighter than the underlying data warrant.

2026-09-01. The docstring of `05.scripts/export_rgb_photo_dates.py` states "2011-12 to 2012-12 no
Landsat at all", which the coverage measurement disproves. Matters because that claim is committed and
will mislead a later reader. See [[landsat7-2012-coverage]].
