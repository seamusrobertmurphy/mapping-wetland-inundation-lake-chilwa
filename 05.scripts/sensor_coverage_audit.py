import ee, json

ee.Initialize(project="murphys-deforisk")

aoi = ee.FeatureCollection(
    json.load(open("/Volumes/PortableSSD/Github/mapping-wetland-inundation-lake-chilwa/03.outputs/JSON/chilwa_basin.geojson"))
).geometry()

for cid, label in [("LANDSAT/LT05/C02/T1_L2", "Landsat 5 TM"),
                   ("LANDSAT/LE07/C02/T1_L2", "Landsat 7 ETM+")]:
    a = ee.ImageCollection(cid).filterBounds(aoi)
    print(f"\n=== {label} over Chilwa basin ===", flush=True)
    print("  total scenes (all dates, no cloud filter):", a.size().getInfo(), flush=True)
    last = a.aggregate_max("system:time_start").getInfo()
    print("  LAST acquisition over Chilwa:", ee.Date(last).format("YYYY-MM-dd").getInfo(), flush=True)

    win = a.filterDate("2011-06-01", "2013-06-30").sort("system:time_start")
    n = win.size().getInfo()
    print(f"  scenes 2011-06-01 to 2013-06-30 (NO cloud filter): {n}", flush=True)
    if n:
        rows = win.reduceColumns(
            ee.Reducer.toList(3),
            ["system:time_start", "CLOUD_COVER", "system:index"]
        ).getInfo()["list"]
        for t, cc, idx in rows:
            import datetime
            d = datetime.datetime.utcfromtimestamp(t / 1000).strftime("%Y-%m-%d")
            print(f"     {d}  cloud={cc}  {idx}", flush=True)
