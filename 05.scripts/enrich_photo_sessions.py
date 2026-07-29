#!/usr/bin/env python3
"""
Narrow the field photograph sessions and mine the filenames for placement hints.

`build_photo_sessions.py` clusters the 244 photographs at a six-hour gap, which
is too loose to stand for one place: the largest cluster holds 152 frames across
an eight-hour day containing a 5.8-hour break, so a single coordinate would
misplace most of it. This pass re-clusters at a tighter gap and reports the
internal structure so the author can see which clusters are genuinely one stop.

It also reads the filenames. Sixty-four of the photographs were renamed by hand
during fieldwork, and those names carry the two things the EXIF lacks: where the
frame was taken and what is in it. Phaloni, the Sombani River, and the zimbowera
platforms all appear by name. The script matches those against a gazetteer of
sites named in the project record and against cover-class keywords, then writes
one row per session with the hints attached and the placement columns left
blank for the author.

Nothing here invents a coordinate. The hints narrow the search; the author still
places each session by hand.

Outputs
  03.outputs/photo_sessions_enriched.csv   one row per session, hints attached
  03.outputs/photo_placement_hints.csv     one row per hinted photograph

Run: python3 05.scripts/enrich_photo_sessions.py
"""
import csv
import datetime as dt
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
OUT = os.path.join(ROOT, "03.outputs")

# A stop, not a day. Ninety minutes holds a single visit together while breaking
# the morning boat trip apart from the afternoon landing.
GAP = dt.timedelta(minutes=90)

# Sites named in the ethnographic record (Methods 2.1, Section 1.3) and in the
# photograph filenames. Matched case-insensitively on word boundaries.
GAZETTEER = [
    "Kachulu", "Phaloni", "Sombani", "Napali", "Andere", "Chambwalu",
    "Lingoni", "Manda Manjeza", "Chisi", "Mapila", "Zomba", "Phalombe",
    "Machinga", "Mpoto", "Zambia",
]

# Filename keywords that point at one of the four classification classes.
# These are hints for the author to confirm against the frame, not labels.
COVER_KEYWORDS = {
    "emergent/flooded vegetation": [
        "reed", "reeds", "marsh", "typha", "swamp", "zimbo", "zimbowera",
        "zimbowere", "fence", "fences",
    ],
    "open water": [
        "boat", "boatman", "canoe", "river", "lake", "water", "net",
        "trap", "traps", "fishing", "fish", "dock", "paddl",
    ],
    "dry vegetation": ["sugarcane", "garden", "maize", "field", "crop"],
    "bare soil": ["beach", "shore", "sand", "lakebed", "mud"],
}

# Camera-default naming. Anything else was renamed by hand and may carry meaning.
DEFAULT_NAME = re.compile(
    r"^(CIMG\d+|DSC_?\d+|IMG_?\d+|P\d{7}|SAM_\d+|PICT\d+|\d{8}_\d{6})",
    re.IGNORECASE,
)

# The fieldwork window given in Methods 2.1. Frames outside it are not
# necessarily wrong, but the camera clock is the first thing to suspect.
FIELD_START = dt.date(2012, 9, 1)
FIELD_END = dt.date(2014, 3, 31)


def site_hits(name):
    hits = []
    for site in GAZETTEER:
        if re.search(r"\b" + re.escape(site) + r"\b", name, re.IGNORECASE):
            hits.append(site)
    return hits


def cover_hits(name):
    low = name.lower()
    hits = []
    for cls, words in COVER_KEYWORDS.items():
        if any(re.search(r"\b" + w, low) for w in words):
            hits.append(cls)
    return hits


def main():
    path = os.path.join(OUT, "field_photo_metadata.csv")
    rows = []
    with open(path) as fh:
        for r in csv.DictReader(fh):
            if not r["datetime"]:
                continue
            r["_t"] = dt.datetime.fromisoformat(r["datetime"])
            rows.append(r)
    rows.sort(key=lambda r: r["_t"])

    sessions, current = [], []
    for r in rows:
        if current and r["_t"] - current[-1]["_t"] > GAP:
            sessions.append(current)
            current = []
        current.append(r)
    if current:
        sessions.append(current)

    hint_rows = []
    with open(os.path.join(OUT, "photo_sessions_enriched.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([
            "session_id", "start", "end", "duration_min", "n_photos", "devices",
            "named_frames", "site_hint", "cover_hint", "date_flag",
            "first_photo", "last_photo",
            "site_name", "lon", "lat", "cover_class", "confidence", "notes",
        ])
        for i, s in enumerate(sessions, 1):
            sites, covers, named = [], [], 0
            for r in s:
                name = os.path.splitext(r["photo"])[0]
                if DEFAULT_NAME.match(r["photo"]):
                    continue
                named += 1
                for h in site_hits(name):
                    if h not in sites:
                        sites.append(h)
                for c in cover_hits(name):
                    if c not in covers:
                        covers.append(c)
                if site_hits(name) or cover_hits(name):
                    hint_rows.append([
                        r["photo"], r["datetime"], i,
                        "; ".join(site_hits(name)), "; ".join(cover_hits(name)),
                    ])

            d = s[0]["_t"].date()
            if d.year == 2011 and d.month == 1:
                flag = "date unusable: Casio epoch default after a clock reset"
            elif d < FIELD_START:
                flag = "predates the Sept 2012 fieldwork start given in Methods 2.1"
            elif d > FIELD_END:
                flag = "postdates the Mar 2014 fieldwork end given in Methods 2.1"
            else:
                flag = ""

            devices = "; ".join(sorted({r["model"] for r in s if r["model"]}))
            dur = round((s[-1]["_t"] - s[0]["_t"]).total_seconds() / 60)
            w.writerow([
                i, s[0]["datetime"], s[-1]["datetime"], dur, len(s), devices,
                named, "; ".join(sites), "; ".join(covers), flag,
                s[0]["photo"], s[-1]["photo"], "", "", "", "", "", "",
            ])

    with open(os.path.join(OUT, "photo_placement_hints.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["photo", "datetime", "session_id", "site_hint", "cover_hint"])
        w.writerows(hint_rows)

    sizes = [len(s) for s in sessions]
    placed = sum(1 for s in sessions
                 if any(site_hits(os.path.splitext(r["photo"])[0]) for r in s))
    covered = sum(len(s) for s in sessions
                  if any(site_hits(os.path.splitext(r["photo"])[0]) for r in s))
    print(f"photographs clustered: {len(rows)}")
    print(f"sessions at a {GAP} gap: {len(sessions)} "
          f"(was 18 at six hours; largest now {max(sizes)} frames, "
          f"{sum(1 for n in sizes if n == 1)} singletons)")
    print(f"sessions carrying a site name: {placed} ({covered} photographs)")
    print(f"photographs with any hint: {len(hint_rows)}")
    flagged = sum(1 for s in sessions if s[0]["_t"].date() == dt.date(2011, 1, 1)
                  or not (FIELD_START <= s[0]["_t"].date() <= FIELD_END))
    print(f"sessions with a suspect date: {flagged}")
    print("\nwrote 03.outputs/photo_sessions_enriched.csv "
          "and 03.outputs/photo_placement_hints.csv")


if __name__ == "__main__":
    main()
