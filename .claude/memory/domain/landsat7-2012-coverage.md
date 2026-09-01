# Landsat 7 coverage over the basin in 2012, measured

2026-09-01. Stacking Landsat 7 fills the SLC-off gaps over this basin almost completely. Measured
with `05.scripts/l7_2012_coverage_test.py` against `03.outputs/JSON/chilwa_basin.geojson`, cloud and
shadow masked from QA_PIXEL: over February to June 2012, 15 scenes leave only 18 km2 of the 8,752 km2
basin never observed, 0.2 per cent, and 95.2 per cent of the basin gets two or more clear looks. Over
the whole of 2012, 41 scenes leave nothing unobserved and 96.7 per cent of the basin gets five or more
clear looks. Matters because it removes the technical basis for treating 2012 as an optical void.

2026-09-01. The SLC-off gap wedges move between acquisitions of the same path and row, so the
CLAUDE.md statement that "the gaps repeat in the same geometry on every pass of a given path" is
wrong. For path 167 row 71 in 2012, one scene observes 6,740 km2 of the basin, the union of all 13
scenes observes 7,789 km2, and only 1,814 km2 is observed in every scene. Matters because the
directional-bias argument for excluding Landsat 7 rested on gaps being static, and they are not.

2026-09-01. Path 166 covers only 627 km2 of the basin, 7.2 per cent; path 167 covers effectively all
of it. The 43 to 48 per cent within-basin pixel loss recorded on 29 July 2026 was measured on path
166, which barely touches the basin. Matters because the closed Landsat 7 exclusion in `DECISIONS.md`
line 15 cites that figure as its evidence.

2026-09-01. Near a single date Landsat 7 remains thin, which is the real constraint. Around
22 February 2012, the date of the 43-photograph flooded-vegetation session, a plus or minus 8 day
window gives 2 scenes covering 22.3 per cent of the basin, and plus or minus 24 days gives 5 scenes
covering 64.5 per cent. Matters because a single-date 2012 map still needs a temporal model or a
coarser sensor, not raw Landsat alone. See [[lemma-landtrendr-gap]].
