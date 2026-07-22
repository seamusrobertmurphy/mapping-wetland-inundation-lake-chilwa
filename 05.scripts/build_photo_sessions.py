#!/usr/bin/env python3
"""
Group the field photographs into locatable sessions.

None of the 244 photographs carries a GPS fix. Four Nikon frames hold a GPS
IFD, but it is empty, which is what the D7000 writes when no GPS unit is
attached. Every photograph does carry a camera timestamp, so location can be
recovered a session at a time rather than a frame at a time: consecutive
photographs separated by less than a set gap were almost certainly taken at one
place, so one coordinate placed by hand fixes the whole group.

This script reads the EXIF timestamps, clusters them into sessions, and writes
a worksheet with one row per session for the author to fill in. It also writes
the per-photograph metadata that `03.outputs/field_photo_metadata.csv` should
have carried; that file reports every photograph as having no GPS IFD, which is
true of the coordinates but hides the device split that decides what is
recoverable.

Outputs
  03.outputs/photo_sessions_worksheet.csv   one row per session, to be completed
  03.outputs/field_photo_metadata.csv       per-photograph EXIF, rewritten

Run: /opt/local/bin/python 05.scripts/build_photo_sessions.py
"""
import csv
import datetime as dt
import glob
import os

from PIL import Image, ExifTags

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
PHOTOS = os.path.join(ROOT, "02.inputs", "PNG")
OUT = os.path.join(ROOT, "03.outputs")

# A session breaks when the gap between consecutive frames exceeds this. Six
# hours keeps a morning and an afternoon at different sites apart while holding
# a single visit together.
GAP = dt.timedelta(hours=6)

# Cameras with no GPS receiver at all: no amount of reprocessing recovers a fix
# from these, so they can only be placed by session.
NO_GPS_HARDWARE = {"EX-Z1000", "DSC-S5000", "CanoScan LiDE 110",
                   "Canon DIGITAL IXUS 430", "C2Z,D520Z,C220Z"}


def exif(path):
    """Timestamp, device, and whether a usable GPS fix is present."""
    try:
        ex = Image.open(path).getexif()
    except Exception:
        return None
    tags = {ExifTags.TAGS.get(k, k): v for k, v in ex.items()}
    gi = ex.get_ifd(0x8825)
    gps = {ExifTags.GPSTAGS.get(k, k): v for k, v in gi.items()} if gi else {}
    lat = lon = None
    if "GPSLatitude" in gps and "GPSLongitude" in gps:
        def dms(v):
            d, m, s = (float(x) for x in v)
            return d + m / 60 + s / 3600
        lat = dms(gps["GPSLatitude"]) * (-1 if gps.get("GPSLatitudeRef") == "S" else 1)
        lon = dms(gps["GPSLongitude"]) * (-1 if gps.get("GPSLongitudeRef") == "W" else 1)
    stamp = None
    for key in ("DateTimeOriginal", "DateTime"):
        if tags.get(key):
            try:
                stamp = dt.datetime.strptime(str(tags[key]), "%Y:%m:%d %H:%M:%S")
                break
            except ValueError:
                pass
    model = str(tags.get("Model", "") or "").strip()
    return {"photo": os.path.basename(path), "datetime": stamp,
            "make": str(tags.get("Make", "") or "").strip(), "model": model,
            "gps_lat": lat, "gps_lon": lon,
            "has_gps": int(lat is not None),
            "gps_capable": int(model not in NO_GPS_HARDWARE and model != "")}


def main():
    files = sorted(glob.glob(os.path.join(PHOTOS, "*.jpg")) +
                   glob.glob(os.path.join(PHOTOS, "*.JPG")))
    rows = [r for r in (exif(f) for f in files) if r]
    dated = sorted((r for r in rows if r["datetime"]), key=lambda r: r["datetime"])
    undated = [r for r in rows if not r["datetime"]]

    # cluster on the time gap
    sessions, current = [], []
    for r in dated:
        if current and r["datetime"] - current[-1]["datetime"] > GAP:
            sessions.append(current)
            current = []
        current.append(r)
    if current:
        sessions.append(current)

    with open(os.path.join(OUT, "field_photo_metadata.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["photo", "datetime", "date", "make", "model",
                    "gps_lat", "gps_lon", "has_gps", "gps_capable", "session_id"])
        sid = {r["photo"]: i + 1 for i, s in enumerate(sessions) for r in s}
        for r in dated + undated:
            d = r["datetime"]
            w.writerow([r["photo"], d.strftime("%Y-%m-%d %H:%M:%S") if d else "",
                        d.date() if d else "", r["make"], r["model"],
                        r["gps_lat"] or "", r["gps_lon"] or "",
                        r["has_gps"], r["gps_capable"], sid.get(r["photo"], "")])

    with open(os.path.join(OUT, "photo_sessions_worksheet.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["session_id", "start", "end", "n_photos", "devices",
                    "first_photo", "last_photo",
                    "site_name", "lon", "lat", "cover_class", "confidence", "notes"])
        for i, s in enumerate(sessions, 1):
            devices = "; ".join(sorted({r["model"] for r in s if r["model"]}))
            w.writerow([i, s[0]["datetime"], s[-1]["datetime"], len(s), devices,
                        s[0]["photo"], s[-1]["photo"], "", "", "", "", "", ""])

    span = [len(s) for s in sessions]
    print(f"photographs: {len(rows)}  dated: {len(dated)}  undated: {len(undated)}")
    print(f"sessions at a {GAP} gap: {len(sessions)}  "
          f"(largest {max(span)} frames, {sum(1 for n in span if n == 1)} singletons)")
    print(f"with a usable GPS fix: {sum(r['has_gps'] for r in rows)}")
    print(f"taken on GPS-capable hardware: {sum(r['gps_capable'] for r in rows)}")
    print("\nsessions by year:")
    years = {}
    for i, s in enumerate(sessions, 1):
        years.setdefault(s[0]["datetime"].year, []).append(len(s))
    for y in sorted(years):
        print(f"  {y}: {len(years[y]):3d} sessions, {sum(years[y]):3d} photographs")
    print("\nwrote 03.outputs/photo_sessions_worksheet.csv "
          "and rewrote 03.outputs/field_photo_metadata.csv")


if __name__ == "__main__":
    main()
