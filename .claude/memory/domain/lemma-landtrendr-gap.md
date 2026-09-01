# LEMMA, LandTrendr and the 2012 imagery gap

2026-09-01. LEMMA (Oregon State) does not treat 2011 to 2013 as a gap, because it uses Landsat 7
throughout. Their project page states the imagery is "Landsat satellite imagery that was 'temporally
normalized' using the LandTrendr algorithms (Kennedy et al. 2010) to produce a yearly time-series of
imagery mosaics from 1984 to 2012", and the CA-Biomass page names "multiple Landsat missions
(Landsat 5, 7, and 8)". Matters because our Landsat 7 exclusion is closed and measured, so LEMMA's
route is unavailable to us on the same terms.

2026-09-01. The gap-covering mechanism is LandTrendr's temporal segmentation, not a fill step.
Kennedy et al. (2010, Remote Sensing of Environment 114, 2897-2910) build each pixel's trajectory
from "the series of single best values for each pixel for each year", so that "clouds, cloud shadows,
and gaps caused by the Landsat 7 scan-line corrector failure can thus be avoided on a pixel-by-pixel
basis", and where the target date is clouded the value is taken "iteratively from the next closest in
date". Straight-line segments are then fitted between vertices, so a weak year reads off the fitted
line rather than off its own observation. Matters because the published value for a thin year is
modelled, and the same idea underlies any defensible 2012 estimate in this project.

2026-09-01. LandTrendr's own paper warns that the first and last years are the weakest, since "the
spectral deviations of the first and last years of the trajectory are by definition more difficult to
judge than deviations in all other years". LEMMA's 20-year Northwest Forest Plan maps are 1993 and
2012, so their 2012 map sits at exactly that endpoint. Matters because it caps how much authority the
LEMMA precedent can carry for our own 2012 estimate.

2026-09-01. Neither Ohmann et al. (2014, RSE 151, 3-15) nor Bell et al. (2015) discusses the
2011 to 2013 window; Ohmann's series stops at 2011 and Bell's runs to 2012, and neither text contains
"SLC" or "gap". Matters because the documentation the author was looking for on the LEMMA site does
not exist there; it is in Kennedy et al. (2010).
