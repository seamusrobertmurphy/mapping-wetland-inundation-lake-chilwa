# Ground truth locating layers

Imagery to relocate the field photographs, on the days they were taken. Copy this
whole folder and drag the GeoTIFFs into QGIS. Everything is EPSG:4326 at 30 m on
one grid, so the layers overlay exactly with no reprojection.

## Two versions of every date, and which to use

Each date has a `_near` and a `_gapfilled` version, and the difference matters.

`_near` composites Landsat within 24 days of the photograph. It is the honest
record of that fortnight. In February 2012 only four Landsat 7 scenes fall in
that window, so its scan-line wedges are plainly visible as diagonal striping.

`_gapfilled` composites within 75 days. The wedges fill, because they move
between passes, at the cost of the image no longer belonging to a single date.

**For recognising a shoreline, a channel or a reed edge from a photograph, use
`_gapfilled`.** You are matching landforms, and a complete image does that
better. **For reading water extent on the day, use `_near`**, and read its gaps
as missing rather than as land.

## Layers

| File | What it is | Display |
|:---|:---|:---|
| `DATE_rgb_near.tif` | True colour, 24-day composite | Already byte-scaled, load as RGB |
| `DATE_rgb_gapfilled.tif` | True colour, 75-day composite | Already byte-scaled, load as RGB |
| `DATE_mndwi_near.tif` | Water index, 24-day | Int16, MNDWI x 10000. Water is above 0 |
| `DATE_mndwi_gapfilled.tif` | Water index, 75-day | Int16, MNDWI x 10000. Water is above 0 |
| `water_occurrence_1984_2021.tif` | Share of the record each pixel held water | Byte 0 to 100, use a single-band pseudocolour ramp |

A sixth layer sits outside this folder and is worth adding:
`../harmonic_2012/mndwi_harmonic_2012-02-22.tif` is the modelled water surface
for 22 February 2012, gap free by construction and validated against Water
Observations from Space at 98.6 to 99.6 per cent agreement. It is the best
available picture of the lake on the day of the 43-frame flooded-vegetation
session.

## Dates and what the frames show

| Date | Frames | Content |
|:---|---:|:---|
| 2011-12-21 | 2 | dry vegetation and exposed lakebed |
| 2012-02-20 | 10 | open water at Phaloni and the Sombani inflow |
| 2012-02-22 | 43 | flooded vegetation and open water |
| 2012-03-29 | 3 | flooded vegetation |
| 2012-04-19 | 5 | open water |
| 2012-05-20 | 4 | flooded vegetation at a landing site |
| 2012-09-26 | 4 | open water |

## The classification scheme

These are the classes the mapping tool offers, and until now they existed only
inside that tool's code. They are written out as
`../../CSV/ground_truth_class_scheme.csv`, with a QGIS categorised style beside
this file as `ground_truth_points.qml`. Load your points, open Properties,
Symbology, Load Style, and pick that file; it keys on an integer field named
`code`.

| Code | Class | Status | Definition | Frames with this evidence |
|---:|:---|:---|:---|---:|
| 1 | open water | in the study | water open to the sky, no emergent stems breaking the surface | 37 |
| 2 | flooded vegetation | in the study | standing water beneath or between live emergent macrophytes | 57 |
| 3 | dry vegetation | in the study | rooted vegetation with no standing water at the surface | 4 |
| 4 | bare soil | in the study | exposed lakebed or unvegetated ground | 2 |
| 8 | rice paddy | excluded | flooded cropland, sparse cover, spectrally ambiguous with bare soil | 1 |
| 9 | settlement | excluded | built ground, market, landing site or processing area | 72 |
| 0 | unassigned | working | placed but not yet classified | 0 |

Codes 1 to 4 are the study's four classes, closed on 5 July 2026 and unchanged.
Codes 8 and 9 exist so that a frame which is plainly none of the four can be
recorded as itself and excluded, rather than forced into a class it does not
belong to. The rice paddy session was flagged at labelling for exactly that
reason.

### What the frame counts say, and it is not comfortable

Of 243 photographed frames, 173 carry land cover evidence and 70 do not, being
archival documents, interiors, gear details and specimens. That is expected.

The problem is the split within the four study classes. Open water and flooded
vegetation are well evidenced at 37 and 57 frames. **Dry vegetation has 4 frames
and bare soil has 2**, and the two bare soil frames come from a single session
on 21 December 2011. A four-class scheme cannot be trained or validated on two
frames from one session, whatever their positions turn out to be. Either those
two classes are evidenced from somewhere else, or the reference set supports two
classes and not four, and the manuscript should say which.

## Two cautions

The water occurrence layer is a locating aid only. It does not enter the
manuscript, per the closed decision that public global products stay out of its
text. Its use here is that it shows the permanent channel network and the lake's
usual margin, which is what a photograph is recognised against.

The Phaloni reference at 35.6984, -15.5721 came from a place name read out of a
filename and geocoded. It carries about 5 km of uncertainty. Treat it as a hint
about where to look, never as a control point.
