#!/usr/bin/env python3
"""Derive candidate coordinates for the field photographs.

The photographs carry no position. This was verified exhaustively rather than
assumed: all 245 frames were read for EXIF, XMP and IPTC, and none carries a
latitude or longitude. Four Nikon frames hold a GPS IFD, but it contains only
GPSVersionID, so a receiver was attached and never wrote a fix.

Position must therefore be reconstructed. Four independent routes exist, in
descending order of reliability:

  A. Differential GPS survey records held by the author. Authoritative; this
     script leaves those rows for them and never invents a coordinate.
  B. Toponyms in the filenames the author assigned. "Sombani River at Phaloni"
     places Phaloni on the Sombani, whose course is mapped, so the river mouth
     bounds it. Reliable to a few kilometres, which is inside a MODIS cell but
     not inside a Landsat pixel.
  C. Horizon matching. Several frames show the hill profile west of the lake
     against open water. A viewshed over the 30 m terrain constrains the camera
     to the locus reproducing that skyline. Not implemented here.
  D. Participatory re-elicitation: show the frames to the communities and have
     them place them. Native to this study's own method.

Confidence is recorded per row and never silently upgraded. Rows without a
defensible coordinate are emitted empty rather than filled with a guess.

Writes 03.outputs/CSV/photo_coordinate_candidates.csv
"""

import csv
import os
import re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CSVD = os.path.join(ROOT, "03.outputs", "CSV")

# Toponyms recoverable from the filenames the author assigned.
TOPONYMS = {
    "phaloni": "Phaloni",
    "kachulu": "Kachulu",
    "sombani": "Sombani River",
    "zambia": "Zambia fishing camp",
}

# Gazetteer. Sources are named per entry so nothing is anonymous.
#   Sombani mouth : OSM way "Sombani", lakeward end of the mapped course
#   Ntila, Nchisi : GeoNames MW gazetteer
# Kachulu is absent from GeoNames and from OSM place nodes in the basin, so it
# carries no coordinate here despite being the lake's main landing site.
GAZETTEER = {
    "Sombani River": (-15.5721, 35.6984, "OSM way Sombani, lakeward end", "3 km"),
    "Phaloni": (-15.5721, 35.6984, "on the Sombani; river mouth as proxy", "5 km"),
    "Ntila": (-14.9833, 35.6333, "GeoNames MW", "2 km"),
    "Nchisi Island": (-15.3333, 35.6333, "GeoNames MW", "2 km"),
    "Kachulu": (None, None, "absent from GeoNames and OSM; needs dGPS or elicitation", ""),
    "Zambia fishing camp": (None, None, "local camp name, not in any gazetteer", ""),
}


def toponyms_in(name):
    low = name.lower()
    return [v for k, v in TOPONYMS.items() if k in low]


def main():
    meta = list(csv.DictReader(open(os.path.join(CSVD, "field_photo_metadata.csv"))))
    labelled = {}
    path = os.path.join(CSVD, "photo_sessions_labelled.csv")
    if os.path.exists(path):
        for r in csv.DictReader(open(path)):
            labelled[r["session_id"]] = r

    out = []
    for m in meta:
        tops = toponyms_in(m["photo"])
        sid = m.get("session_id", "")
        sess = labelled.get(sid, {})
        lon = lat = ""
        method = "unresolved"
        conf = ""
        src = ""
        if tops:
            # Prefer the most specific toponym that carries a coordinate.
            for t in tops:
                g = GAZETTEER.get(t)
                if g and g[0] is not None:
                    lat, lon = g[0], g[1]
                    src, conf = g[2], g[3]
                    method = "toponym in filename"
                    break
            else:
                method = "toponym found, no gazetteer entry"
                src = "; ".join(tops)
        out.append(
            {
                "photo": m["photo"],
                "session_id": sid,
                "date": m["date"],
                "model": m["model"],
                "cover_class_reviewed": sess.get("cover_class_reviewed", ""),
                "toponym": "; ".join(tops),
                "lon": lon,
                "lat": lat,
                "method": method,
                "positional_uncertainty": conf,
                "source": src,
                "dgps_lon": "",
                "dgps_lat": "",
            }
        )

    dest = os.path.join(CSVD, "photo_coordinate_candidates.csv")
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)

    res = sum(1 for r in out if r["lat"] != "")
    print(f"{len(out)} frames written to {dest}")
    print(f"  {res} carry a derived coordinate")
    print(f"  {sum(1 for r in out if r['toponym'])} carry a toponym")
    print(f"  {len(out) - res} await dGPS, horizon matching, or elicitation")


if __name__ == "__main__":
    main()
