# Lake Chilwa 2012 layers for QGIS

Open `chilwa-2012.qgs` from this folder. Paths inside it are relative, so the
project only works while it sits here alongside `ground_truth_aid/` and
`harmonic_2012/`. Only the 2012-02-20 true colour layer is switched on; tick the
rest in the layer panel.

Fifteen layers, all EPSG:4326, all 30 m Landsat 7 ETM+ surface reflectance
Collection 2 Level 2 unless stated.

| Layer | File | Bands |
|---|---|---|
| True colour, six 2012 dates | `ground_truth_aid/2012-*_rgb.tif` | 3, byte, stretched for viewing |
| MNDWI, the same six dates | `ground_truth_aid/2012-*_mndwi.tif` | 1, scaled by 10,000 |
| MNDWI harmonic fit | `harmonic_2012/mndwi_harmonic_2012-02-22.tif` | 1 |
| Clear observation count | `harmonic_2012/harmonic_n_obs.tif` | 1 |
| Water occurrence 1984 to 2021 | `ground_truth_aid/water_occurrence_1984_2021.tif` | 1, per cent |

## The stripes

They are Scan Line Corrector gaps. The corrector on Landsat 7 ETM+ failed on
31 May 2003 and every acquisition since is SLC-off, losing about 22 per cent of
a scene in wedges that are zero at nadir and widen toward the swath edges.

The gaps are missing data, not wrong data. Measured on
`ground_truth_aid/2012-02-20_rgb.tif`, 31.24 per cent of the frame is NoData and
exactly 0.00 per cent carries the value zero, so nothing in the gaps enters an
average, an index or an area. The black you see is the viewer painting NoData,
not a reflectance of zero.

## Why Landsat 7 for 2012

Because there was nothing else at 30 m. From `03.outputs/CSV/optical_archive_survey.csv`,
the last year Landsat 5 TM appears over this basin is 2011 and the first year
Landsat 8 OLI appears is 2013. For 2012 the survey returns one sensor,
`L7_ETM`, with 41 scenes, 26 of them under 30 per cent cloud.

Per-scene loss for the 2012 renders is recorded in `02.inputs/RGB/README.md` and
runs from 10.1 to 48.3 per cent, independent of cloud cover. Path 166 row 071
loses over 40 per cent on every date; path 167 loses 10 to 16 per cent. That is
geometry, not weather, and it is why the composite draws on four scenes.

## What to check before trusting a number

Load `harmonic_2012/harmonic_n_obs.tif` under any result. Where the clear
observation count is low the harmonic fit is extrapolating, and a shoreline
drawn there is a model output rather than an observation.
