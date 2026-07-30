#!/opt/local/bin/python
"""Rebuild the field-photograph metadata on EXIF capture time.

The earlier extraction read EXIF tag 0x0132 (DateTime), which records the last
file modification, not the exposure. For 143 of 245 files that value is a
Photoshop resave or a device transfer, so the session structure derived from it
grouped 140 frames into a 13-minute "session" on the evening of 19 June 2012
that never happened in the field.

This script reads tag 0x9003 (DateTimeOriginal), falls back to 0x9004
(DateTimeDigitized), and only then to DateTime, recording which field supplied
each timestamp. It collapses derivative copies to distinct capture instants,
rebuilds sessions from exposure times, and flags the camera-clock reset block.

Outputs, all to 03.outputs/CSV:
  field_photo_metadata.csv   one row per file, corrected timestamps
  photo_captures.csv         one row per distinct capture instant
  photo_sessions.csv         one row per rebuilt session
"""

import csv
import os
import re
from collections import Counter, defaultdict
from datetime import datetime, timedelta

from PIL import Image
from PIL.ExifTags import TAGS

REPO = "/Volumes/PortableSSD/Github/mapping-wetland-inundation-lake-chilwa"
PNG = os.path.join(REPO, "02.inputs", "PNG")
OUT = os.path.join(REPO, "03.outputs", "CSV")

# A new session starts after this long a break in exposures, or on a date change.
SESSION_GAP = timedelta(minutes=90)

# The Casio EX-Z1000 clock was reset to a 2011-01-01 00:00 epoch and then ran
# freely for some days before being corrected. Frames stamped 2011-01-01 are the
# obvious cases, but CIMG1362 to CIMG1368, stamped 10 and 11 January 2011,
# continue the same filename sequence and precede CIMG1804 on 22 February 2012
# with its clock correct. The whole January 2011 block from that camera is
# therefore reset-era and its dates are unrecoverable.
CLOCK_RESET_MONTH = "2011-01"
CLOCK_RESET_MODEL = "EX-Z1000"

# Flatbed scanners digitise prints. The EXIF date is the scan, not the exposure,
# so these frames carry no field date at all.
SCANNER_MODELS = ("CanoScan",)

# Some descriptive filenames carry a date the EXIF has lost.
FILENAME_DATE = re.compile(r"(\d{1,2})-(\d{1,2})-(\d{2,4})")

# Toponyms that recur in the descriptive filenames.
TOPONYMS = [
    "Phaloni", "Sombani", "Kachulu", "Zomba", "Machinga", "Phalombe",
    "Mtengo", "Ntila", "Andere", "Chambwalu", "Lingoni", "Zambia",
    "Likangala", "Namasika", "Mposa", "Swang'oma", "Chisi",
]

COVER_TERMS = [
    (r"zimbo|zimbowere|zimbowera|reed|marsh|typha", "emergent/flooded vegetation"),
    (r"fish\s*fence|fish\s*trap|fishtrap|boat|canoe|water|river|lake", "open water"),
    (r"garden|maize|field|crop|tree|bank|grass", "dry vegetation"),
    (r"lakebed|mud|sand|dry\s*bed|salt", "bare soil"),
]

EXIF_FIELDS = ("DateTimeOriginal", "DateTimeDigitized", "DateTime")


def read_exif(path):
    """Return a dict of the EXIF tags this script needs, empty on failure."""
    try:
        with Image.open(path) as im:
            size = im.size
            raw = im._getexif() or {}
    except Exception:
        return {"width": 0, "height": 0}
    tags = {TAGS.get(k, k): v for k, v in raw.items()}
    out = {"width": size[0], "height": size[1]}
    for f in EXIF_FIELDS + ("Make", "Model", "Software"):
        v = tags.get(f)
        out[f] = str(v).strip() if v is not None else ""
    # A GPSInfo block is not a fix. Four Nikon frames carry a GPSVersionID stub
    # with no coordinates, because the body writes the block whether or not a
    # receiver is attached. Require an actual latitude and longitude.
    gps = tags.get("GPSInfo")
    lat = lon = None
    if isinstance(gps, dict):
        lat, lon = gps.get(2), gps.get(4)
    out["has_gps_fix"] = bool(lat and lon)
    out["gps_block_only"] = bool(gps) and not (lat and lon)
    return out


def parse_exif_dt(s):
    """EXIF datetimes are 'YYYY:MM:DD HH:MM:SS'. Null stamps use zeros."""
    if not s or s.startswith("0000"):
        return None
    try:
        return datetime.strptime(s[:19], "%Y:%m:%d %H:%M:%S")
    except ValueError:
        return None


def pick_timestamp(ex):
    """Prefer exposure time; report which field was used so the fallback is visible."""
    for field in EXIF_FIELDS:
        dt = parse_exif_dt(ex.get(field, ""))
        if dt:
            return dt, field
    return None, ""


def filename_date(name):
    """Recover a date the EXIF has lost, e.g. 'Phaloni 5-3-13.JPG'. Ambiguous in
    day-month order, so it is reported as a hint for the author to resolve."""
    m = FILENAME_DATE.search(name)
    if not m:
        return ""
    a, b, y = m.groups()
    y = int(y) + 2000 if len(y) == 2 else int(y)
    return f"{y}-{int(b):02d}-{int(a):02d} or {y}-{int(a):02d}-{int(b):02d}"


def site_hint(name):
    found = [t for t in TOPONYMS if re.search(re.escape(t), name, re.I)]
    return "; ".join(dict.fromkeys(found))


def cover_hint(name):
    found = [label for pat, label in COVER_TERMS if re.search(pat, name, re.I)]
    return "; ".join(dict.fromkeys(found))


def main():
    files = sorted(
        f for f in os.listdir(PNG)
        if f.lower().endswith((".jpg", ".jpeg", ".png")) and not f.startswith("._")
    )

    rows = []
    for f in files:
        ex = read_exif(os.path.join(PNG, f))
        dt, src = pick_timestamp(ex)
        orig = parse_exif_dt(ex.get("DateTimeOriginal", ""))
        modt = parse_exif_dt(ex.get("DateTime", ""))
        rows.append({
            "photo": f,
            "datetime": dt.isoformat(sep=" ") if dt else "",
            "date": dt.date().isoformat() if dt else "",
            "timestamp_source": src,
            "modify_time": modt.isoformat(sep=" ") if modt else "",
            "modify_differs": bool(orig and modt and orig != modt),
            "make": ex.get("Make", ""),
            "model": ex.get("Model", ""),
            "software": ex.get("Software", ""),
            "width": ex.get("width", 0),
            "height": ex.get("height", 0),
            "has_gps": ex.get("has_gps_fix", False),
            "gps_block_only": ex.get("gps_block_only", False),
            "clock_reset": bool(dt)
                and dt.strftime("%Y-%m") == CLOCK_RESET_MONTH
                and CLOCK_RESET_MODEL in ex.get("Model", ""),
            "scanned": any(s in ex.get("Model", "") for s in SCANNER_MODELS),
            "filename_date_hint": filename_date(f),
            "site_hint": site_hint(f),
            "cover_hint": cover_hint(f),
        })

    # Group into distinct capture instants. Two files sharing an exposure second
    # and a camera model are the same frame, one a derivative of the other.
    groups = defaultdict(list)
    undated = []
    for r in rows:
        if r["datetime"]:
            groups[(r["datetime"], r["model"])].append(r)
        else:
            undated.append(r)

    # Canonical frame per instant: the largest by pixel area, i.e. the least
    # cropped or downsampled. Ties break on the shortest filename.
    captures = []
    for (dtstr, model), members in sorted(groups.items()):
        members.sort(key=lambda r: (-(r["width"] * r["height"]), len(r["photo"])))
        canon = members[0]
        aliases = [m["photo"] for m in members[1:]]
        sites = [m["site_hint"] for m in members if m["site_hint"]]
        covers = [m["cover_hint"] for m in members if m["cover_hint"]]
        captures.append({
            "capture_id": 0,
            "datetime": dtstr,
            "date": canon["date"],
            "model": model,
            "canonical_photo": canon["photo"],
            "n_files": len(members),
            "derivative_files": "; ".join(aliases),
            "width": canon["width"],
            "height": canon["height"],
            "timestamp_source": canon["timestamp_source"],
            "clock_reset": canon["clock_reset"],
            "scanned": canon["scanned"],
            "date_reliable": bool(canon["date"]) and not canon["clock_reset"]
                and not canon["scanned"]
                and canon["timestamp_source"] == "DateTimeOriginal",
            "site_hint": "; ".join(dict.fromkeys("; ".join(sites).split("; "))) if sites else "",
            "cover_hint": "; ".join(dict.fromkeys("; ".join(covers).split("; "))) if covers else "",
        })
    for i, c in enumerate(captures, start=1):
        c["capture_id"] = i

    # Map every file back to its capture, so the two tables join.
    file_to_capture = {}
    for c in captures:
        file_to_capture[c["canonical_photo"]] = c["capture_id"]
        for a in [x for x in c["derivative_files"].split("; ") if x]:
            file_to_capture[a] = c["capture_id"]
    for r in rows:
        r["capture_id"] = file_to_capture.get(r["photo"], "")
        r["is_canonical"] = r["photo"] in {c["canonical_photo"] for c in captures}

    # Rebuild sessions from exposure times. The clock-reset block is held apart:
    # its internal ordering is real but its date is not, so it cannot be placed
    # in the chronology alongside the dated sessions.
    dated = [c for c in captures if not c["clock_reset"]]
    reset = [c for c in captures if c["clock_reset"]]

    sessions = []

    def close(members, label):
        if not members:
            return
        ts = [datetime.fromisoformat(m["datetime"]) for m in members]
        sites = [m["site_hint"] for m in members if m["site_hint"]]
        covers = []
        for m in members:
            covers += [x for x in m["cover_hint"].split("; ") if x]
        sessions.append({
            "session_id": len(sessions) + 1,
            "kind": label,
            "start": min(ts).isoformat(sep=" "),
            "end": max(ts).isoformat(sep=" "),
            "duration_min": round((max(ts) - min(ts)).total_seconds() / 60, 1),
            "n_captures": len(members),
            "n_files": sum(m["n_files"] for m in members),
            "devices": "; ".join(sorted({m["model"] for m in members if m["model"]})),
            "first_capture": members[0]["canonical_photo"],
            "last_capture": members[-1]["canonical_photo"],
            "site_hint": "; ".join(dict.fromkeys("; ".join(sites).split("; "))) if sites else "",
            "cover_hint": "; ".join(dict.fromkeys(covers)),
            "site_name": "", "lon": "", "lat": "",
            "cover_class": "", "position_uncertainty_m": "", "confidence": "", "notes": "",
        })

    run = []
    prev = None
    for c in sorted(dated, key=lambda x: x["datetime"]):
        t = datetime.fromisoformat(c["datetime"])
        if prev is not None and (t - prev > SESSION_GAP or t.date() != prev.date()):
            close(run, "field")
            run = []
        run.append(c)
        prev = t
    close(run, "field")
    close(sorted(reset, key=lambda x: x["datetime"]), "clock-reset (date unrecoverable)")

    for c in captures:
        c["session_id"] = ""
    by_dt = {c["datetime"]: c for c in captures}
    for s in sessions:
        for c in captures:
            if s["start"] <= c["datetime"] <= s["end"] and (
                (c["clock_reset"] and s["kind"].startswith("clock"))
                or (not c["clock_reset"] and s["kind"] == "field")
            ):
                c["session_id"] = s["session_id"]

    os.makedirs(OUT, exist_ok=True)

    def write(name, records, fields):
        with open(os.path.join(OUT, name), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields, quoting=csv.QUOTE_MINIMAL)
            w.writeheader()
            for r in records:
                w.writerow({k: r.get(k, "") for k in fields})

    write("field_photo_metadata.csv", rows, [
        "photo", "capture_id", "session_id", "is_canonical", "datetime", "date",
        "timestamp_source", "modify_time", "modify_differs", "clock_reset",
        "scanned", "filename_date_hint", "make", "model", "software",
        "width", "height", "has_gps", "gps_block_only",
        "site_hint", "cover_hint"])

    write("photo_captures.csv", captures, [
        "capture_id", "session_id", "datetime", "date", "model", "canonical_photo",
        "n_files", "derivative_files", "width", "height", "timestamp_source",
        "clock_reset", "scanned", "date_reliable", "site_hint", "cover_hint"])

    write("photo_sessions.csv", sessions, [
        "session_id", "kind", "start", "end", "duration_min", "n_captures", "n_files",
        "devices", "first_capture", "last_capture", "site_hint", "cover_hint",
        "site_name", "lon", "lat", "cover_class", "position_uncertainty_m",
        "confidence", "notes"])

    # Report
    print(f"files scanned                 {len(rows)}")
    print(f"  timestamped                 {sum(1 for r in rows if r['datetime'])}")
    print(f"  no usable timestamp         {len(undated)}")
    print(f"  modify time differs         {sum(1 for r in rows if r['modify_differs'])}")
    print(f"  fell back off DateTimeOriginal "
          f"{sum(1 for r in rows if r['timestamp_source'] not in ('DateTimeOriginal', ''))}")
    print(f"  with a GPS fix              {sum(1 for r in rows if r['has_gps'])}")
    print(f"  empty GPS block, no fix     {sum(1 for r in rows if r['gps_block_only'])}")
    print(f"  scanner-digitised prints    {sum(1 for r in rows if r['scanned'])}")
    print()
    print(f"distinct capture instants     {len(captures)}")
    print(f"  with derivative copies      {sum(1 for c in captures if c['n_files'] > 1)}")
    print(f"  redundant files collapsed   {sum(c['n_files'] - 1 for c in captures)}")
    print(f"  clock-reset captures        {sum(1 for c in captures if c['clock_reset'])}")
    print(f"  scanned prints              {sum(1 for c in captures if c['scanned'])}")
    print(f"  RELIABLE field dates        {sum(1 for c in captures if c['date_reliable'])}")
    print()
    print(f"sessions rebuilt              {len(sessions)}")
    print()
    print("captures by month (exposure time)")
    for k, v in sorted(Counter(c["date"][:7] for c in captures).items()):
        print(f"  {k}  {v:4d}")
    print()
    wet = sum(1 for c in captures if c["date"][:7] in
              ("2012-02", "2012-03", "2012-04", "2012-05"))
    dry = sum(1 for c in captures if c["date"][:7] in
              ("2012-06", "2012-07", "2012-08", "2012-09", "2012-10"))
    print(f"2012 wet season Feb-May       {wet}")
    print(f"2012 dry season Jun-Oct       {dry}")
    print()
    fv = [c for c in captures if "emergent" in c["cover_hint"]]
    print(f"captures hinting flooded veg  {len(fv)}")
    print(f"  in dated sessions           {sum(1 for c in fv if not c['clock_reset'])}")
    print(f"  in the clock-reset block    {sum(1 for c in fv if c['clock_reset'])}")


if __name__ == "__main__":
    main()
