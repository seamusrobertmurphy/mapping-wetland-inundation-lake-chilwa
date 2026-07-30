# Landsat 7 ETM+ 2012 true-colour renders, Lake Chilwa basin

Bands SR_B3, SR_B2, SR_B1 as red, green, blue, surface reflectance scaled by 0.0000275 with an offset of -0.2, stretched 0 to 0.22, clipped to the intersection of each scene footprint with `03.outputs/JSON/chilwa_basin.geojson`, rendered at 1600 px on the long edge.

`data lost` is the percentage of pixels inside that overlap carrying no data, measured as the mean of the SR_B4 mask at 30 m. It counts the Scan Line Corrector gaps, which the `CLOUD_COVER` field does not record. The two columns are independent: a scene can be almost cloud-free and still lose a fifth of its pixels.

| scene | cloud cover | data lost to gaps | file |
|---|---|---|---|
| LE07_166071_20120327 | 9 % | 48.3 % | `LE07_166071_20120327_cloud9pc_dataloss48pc.png` |
| LE07_167070_20120419 | 7 % | 10.1 % | `LE07_167070_20120419_cloud7pc_dataloss10pc.png` |
| LE07_166071_20120428 | 12 % | 45.8 % | `LE07_166071_20120428_cloud12pc_dataloss46pc.png` |
| LE07_167071_20120505 | 3 % | 15.5 % | `LE07_167071_20120505_cloud3pc_dataloss16pc.png` |
| LE07_167070_20120521 | 0 % | 13.8 % | `LE07_167070_20120521_cloud0pc_dataloss14pc.png` |
| LE07_167071_20120521 | 1 % | 16.2 % | `LE07_167071_20120521_cloud1pc_dataloss16pc.png` |
| LE07_166071_20120514 | 8 % | 43.1 % | `LE07_166071_20120514_cloud8pc_dataloss43pc.png` |
| LE07_167071_20120606 | 14 % | 13.4 % | `LE07_167071_20120606_cloud14pc_dataloss13pc.png` |
