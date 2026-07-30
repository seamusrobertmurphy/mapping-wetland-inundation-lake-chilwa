import ee, json, urllib.request, os

ee.Initialize(project="murphys-deforisk")

OUT = "/Volumes/PortableSSD/Github/mapping-wetland-inundation-lake-chilwa/02.inputs/L7_2012_RGB"
os.makedirs(OUT, exist_ok=True)

basin = ee.FeatureCollection(
    json.load(open("/Volumes/PortableSSD/Github/mapping-wetland-inundation-lake-chilwa/03.outputs/JSON/chilwa_basin.geojson"))
).geometry()

SCENES = [
    "LE07_166071_20120327",
    "LE07_167070_20120419",
    "LE07_166071_20120428",
    "LE07_167071_20120505",
    "LE07_167070_20120521",
    "LE07_167071_20120521",
    "LE07_166071_20120514",
    "LE07_167071_20120606",
]

log = []
for sid in SCENES:
    img = ee.Image(f"LANDSAT/LE07/C02/T1_L2/{sid}")
    cc = img.get("CLOUD_COVER").getInfo()
    region = basin.intersection(img.geometry(), 1)

    # fraction of pixels carrying data inside the basin/scene overlap
    valid = img.select("SR_B4").mask().rename("v")
    f = valid.reduceRegion(ee.Reducer.mean(), region, 30, maxPixels=1e10).getInfo()["v"]
    lost = 100 - f * 100

    rgb = img.select(["SR_B3", "SR_B2", "SR_B1"]).multiply(0.0000275).add(-0.2)
    url = rgb.getThumbURL({"min": 0.0, "max": 0.22, "region": region, "dimensions": 1600})
    name = f"{sid}_cloud{int(cc)}pc_dataloss{lost:.0f}pc.png"
    urllib.request.urlretrieve(url, os.path.join(OUT, name))

    print(f"{sid}  cloud={cc}%  data lost to SLC-off gaps = {lost:.1f}%  -> {name}", flush=True)
    log.append((sid, cc, lost, name))

with open(os.path.join(OUT, "README.md"), "w") as fh:
    fh.write("# Landsat 7 ETM+ 2012 true-colour renders, Lake Chilwa basin\n\n")
    fh.write("Bands SR_B3, SR_B2, SR_B1 as red, green, blue, surface reflectance scaled by 0.0000275 with an offset of -0.2, ")
    fh.write("stretched 0 to 0.22, clipped to the intersection of each scene footprint with `03.outputs/JSON/chilwa_basin.geojson`, ")
    fh.write("rendered at 1600 px on the long edge.\n\n")
    fh.write("`data lost` is the percentage of pixels inside that overlap carrying no data, measured as the mean of the ")
    fh.write("SR_B4 mask at 30 m. It counts the Scan Line Corrector gaps, which the `CLOUD_COVER` field does not record. ")
    fh.write("The two columns are independent: a scene can be almost cloud-free and still lose a fifth of its pixels.\n\n")
    fh.write("| scene | cloud cover | data lost to gaps | file |\n|---|---|---|---|\n")
    for sid, cc, lost, name in log:
        fh.write(f"| {sid} | {cc} % | {lost:.1f} % | `{name}` |\n")
print("\nwrote", os.path.join(OUT, "README.md"))
