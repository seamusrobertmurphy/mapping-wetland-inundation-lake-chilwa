# Literature Review: Spectral Water Indices for Inland Wetland Mapping

## Context

This review evaluates the five multispectral water-extraction indices under consideration for mapping inundation dynamics in the Lake Chilwa Basin: NDWI, MNDWI, AWEIsh, WRI, and NDPI. For each index the review covers the original formulation, demonstrated strengths, and known limitations, with particular attention to the conditions that define Lake Chilwa: shallow, turbid, vegetated, endorheic, and subject to dramatic recession-refilling cycles. The review builds the methodological case for why no single optical index suffices in this environment and why SAR integration is necessary.

------------------------------------------------------------------------

## 1. Normalized Difference Water Index (NDWI)

**Formulation:** (Green - NIR) / (Green + NIR) **Reference:** McFeeters (1996), *International Journal of Remote Sensing*, 17(7), 1425-1432.

McFeeters designed NDWI to delineate open water features by exploiting the strong absorption of near-infrared radiation by water. The index suppresses vegetation and soil noise while enhancing water signals and remains widely used for mapping lakes, rivers, and reservoirs where water is deep and relatively clear. Ji, Zhang, and Wylie (2009) confirmed its reliability for surface water extent detection in clear-water lakes across the US Great Plains.

Three limitations bear directly on inland wetland applications. First, NDWI struggles with shallow, turbid water because suspended sediment increases NIR reflectance, compressing index values toward zero and producing false negatives (Xu, 2006; Li et al., 2013). Second, in vegetated wetlands the green-NIR ratio cannot separate water beneath emergent vegetation from surrounding marshland, since canopy reflectance dominates both bands (Ozesmi and Bauer, 2002). Third, endorheic lakes with fluctuating shorelines and extensive littoral zones present mixed pixels where NDWI thresholds misclassify wet soil and shallow margins as dry land (Deus and Gloaguen, 2013).

In African contexts, Malahlela (2016) applied NDWI to South African wetlands and found it underperformed in seasonal floodplain systems compared to SWIR-based indices. Ogilvie et al. (2015) mapped West African floodplain inundation and reported that MNDWI improved discrimination of turbid and shallow water by 15-20% in overall accuracy over NDWI.

**Assessment for Lake Chilwa:** Useful as a baseline measure for permanent open water. Will systematically underestimate inundation extent in the shallow, turbid, vegetated littoral zones and Typha marshes that characterise the lake's margins.

------------------------------------------------------------------------

## 2. Modified Normalized Difference Water Index (MNDWI)

**Formulation:** (Green - SWIR1) / (Green + SWIR1) **Reference:** Xu (2006), *International Journal of Remote Sensing*, 27(14), 3025-3033.

Xu developed MNDWI to address a specific failure of NDWI: confusion between open water and built-up surfaces or bare soil, both of which reflect similarly in the NIR. SWIR reflectance from water is near zero while built-up areas and bare soil reflect strongly, so the substitution suppresses that confusion and yields higher positive values for water bodies. MNDWI excels at separating water from impervious surfaces, performs better than NDWI in turbid or sediment-laden water, and maintains sensitivity across a range of water depths (Li et al., 2013; Acharya et al., 2018).

The index falters in shallow vegetated wetlands. Emergent macrophytes dominate the spectral signal in mixed pixels, driving MNDWI negative even where standing water persists beneath the canopy (Amani et al., 2020). In endorheic and saline systems, dissolved salts and algal blooms alter water's spectral properties, shifting MNDWI thresholds unpredictably (Pekel et al., 2016). The standard zero threshold rarely holds in these contexts; adaptive or Otsu-based thresholding is typically required.

Heimhuber et al. (2018) demonstrated in the Murray-Darling system that MNDWI outperformed NDWI in turbid floodwaters but required SAR fusion where vegetation obscured the water surface.

**Assessment for Lake Chilwa:** The strongest single optical index for this study's conditions. Will capture open water and moderately turbid zones better than NDWI. Will still fail in the dense Typha marshes and at the fluctuating land-water boundary during recession periods.

------------------------------------------------------------------------

## 3. Automated Water Extraction Index (AWEIsh)

**Formulation:** Blue + 2.5 x Green - 1.5 x (NIR + SWIR1) - 0.25 x SWIR2 **Reference:** Feyisa et al. (2014), *Remote Sensing of Environment*, 140, 23-35.

Feyisa designed AWEIsh as a five-band linear combination with coefficients optimised to suppress commission errors from shadow, dark surfaces, and built-up areas that confound simpler indices. The multi-band formulation captures spectral contrasts that two-band normalised indices cannot resolve, particularly the distinction between shadow and open water in the SWIR region. Feyisa et al. (2014) reported consistently higher overall accuracy than NDWI and MNDWI across diverse test sites.

The index performs less reliably in shallow, vegetated, or turbid water. Ji et al. (2015) found that subpixel mixing of emergent vegetation and water reduces AWEIsh values below detection thresholds, a problem shared with NDWI and MNDWI but not resolved by the additional bands. Fisher et al. (2016) compared multiple water indices globally and noted that no single index consistently mapped shallow or seasonal water bodies. In endorheic systems where salinity raises turbidity and alters spectral signatures, AWEIsh can misclassify exposed lakebed or brine as water (Pekel et al., 2016).

Oliphant et al. (2019) mapped surface water across sub-Saharan Africa and found AWEIsh effective for permanent open water but less so for seasonal or marshy areas.

**Assessment for Lake Chilwa:** Strongest performance where shadow or topographic interference affects classification (basin margins, elevated terrain). In the shallow, vegetated lake centre and littoral zones, AWEIsh offers no clear advantage over MNDWI and may introduce noise from mixed vegetation pixels through its additional coefficients.

------------------------------------------------------------------------

## 4. Water Ratio Index (WRI)

**Formulation:** (Green + Red) / (NIR + SWIR1) **Reference:** Shen and Li (2010), *Proceedings of the 18th International Conference on Geoinformatics*, Beijing.

WRI is a ratio index rather than a normalised difference, yielding values greater than 1.0 for water pixels and below 1.0 for most land covers. It suppresses cloud and shadow noise more effectively than NDWI or AWEIsh because clouds and shadows produce low, undifferentiated values in the ratio image (Li et al., 2023). Demelash et al. (2025) reported producer accuracy of 97.98% and kappa of 0.88 for surface water extraction from Landsat 8 in Ethiopia, outperforming both MNDWI and NDWI in landscapes dominated by bare soil with sparse vegetation.

The principal limitation is mixed pixels: where water, vegetation, and sediment share a pixel, the visible-band numerator is inflated by soil and vegetation reflectance, reducing contrast. In shallow turbid water the red band reflectance rises sharply from suspended sediment, compressing the ratio toward ambiguity. WRI also lacks a standard threshold scale, forcing site-specific calibration. It has not been widely applied to endorheic or fluctuating lake systems, and no published application to Lake Chilwa or comparable African endorheic wetlands was found.

**Assessment for Lake Chilwa:** Potentially useful for the bare-soil and sparse-vegetation zones of the exposed lakebed during recession periods. Untested in this specific environment and will require site-specific threshold calibration. The red-band sensitivity to suspended sediment is a liability in Lake Chilwa's turbid waters.

------------------------------------------------------------------------

## 5. Normalized Difference Pond Index (NDPI)

**Formulation:** (SWIR1 - Green) / (SWIR1 + Green) **Reference:** Lacaux et al. (2007), *Remote Sensing of Environment*, 106, 66-80.

Lacaux designed NDPI using SPOT-5 imagery (10 m) to classify temporary ponds in the Ferlo region of Senegal as part of a Rift Valley Fever early warning system. The index exploits the strong SWIR absorption of water against the moderate green reflectance of vegetation, making it sensitive to the water-vegetation boundary. It was designed to separate small water bodies (surface area greater than 100 sq m) from surrounding vegetation, a task where NDWI and MNDWI routinely fail because green vegetation and shallow water produce similar NIR responses.

NDPI has been applied in sub-Saharan epidemiological studies tracking mosquito breeding habitat and in the Okavango endorheic basin for land use change analysis. Campos, Sillero, and Brito (2012) confirmed that SWIR-based indices outperform NIR-based ones for ephemeral ponds in the Sahara-Sahel transition zone. Its limitations mirror its design: it performs best on small, shallow, vegetated ponds and loses sensitivity in deep or turbid water where SWIR absorption saturates. Dense emergent vegetation can mask the water signal entirely.

Note that NDPI is the algebraic inverse of MNDWI (the bands are swapped in the numerator). The two indices carry identical information; NDPI simply reverses the sign so that water yields negative values and vegetation yields positive ones.

**Assessment for Lake Chilwa:** The most directly relevant index for the seasonal, vegetated margins of the lake and surrounding ponds. Its African provenance (Senegal, Okavango) and design for temporary water in semi-arid savanna make it ecologically appropriate. However, it will lose sensitivity in the deeper open water and during high-turbidity conditions that follow recession events.

------------------------------------------------------------------------

## Synthesis: The Case for Multi-Index and SAR Integration

The literature reveals a consistent pattern across all five indices. Each performs well within a defined niche but fails at the boundaries that define Lake Chilwa's hydrology:

| Condition                  | NDWI | MNDWI    | AWEIsh   | WRI      | NDPI     |
|----------------------------|------|----------|----------|----------|----------|
| Deep open water            | Good | Good     | Good     | Good     | Poor     |
| Shallow turbid water       | Poor | Moderate | Moderate | Poor     | Moderate |
| Vegetated wetland          | Poor | Poor     | Poor     | Poor     | Moderate |
| Exposed lakebed            | Poor | Moderate | Moderate | Moderate | Moderate |
| Mixed pixels (littoral)    | Poor | Poor     | Poor     | Poor     | Moderate |
| Shadow/topographic effects | Poor | Moderate | Good     | Moderate | Poor     |

No single index captures the full range of conditions that Lake Chilwa presents across its recession-refilling cycle. The shallow, turbid, and vegetated margins that constitute the most ecologically and socially significant zones of the basin are precisely where all optical indices underperform. This is the fundamental justification for SAR integration: C-band backscatter detects surface water beneath vegetation canopies and is insensitive to the atmospheric and spectral interference that degrades optical indices in this environment.

The multi-index approach adopted in this study (testing all five indices against SAR-derived water maps) provides both a comparative evaluation of optical methods and a quantitative measure of what optical sensors miss in endorheic wetland systems.

------------------------------------------------------------------------

## Key References

-   Acharya, T.D., Subedi, A. and Lee, D.H. (2018) Evaluation of water indices for surface water extraction in a Landsat 8 scene of Nepal. *Sensors*, 18(8), 2580.
-   Amani, M. et al. (2020) Canadian wetland inventory using Google Earth Engine. *Remote Sensing*, 12, 1190.
-   Campos, J.C., Sillero, N. and Brito, J.C. (2012) Normalized difference water indexes have dissimilar performances in detecting seasonal and permanent water in the Sahara-Sahel transition zone. *Journal of Hydrology*, 464-465, 438-446.
-   Demelash, T. et al. (2025) Evaluation of water extraction indices for spatial mapping of surface water bodies using Sentinel-2: the case of Ethiopia. *H2Open Journal*, 8(5), 402 ff.
-   Deus, D. and Gloaguen, R. (2013) Remote sensing analysis of lake dynamics in semi-arid regions: implication for water resource management. *Remote Sensing*, 5, 5765-5781.
-   Feyisa, G.L. et al. (2014) Automated water extraction index: a new technique for surface water mapping using Landsat imagery. *Remote Sensing of Environment*, 140, 23-35.
-   Fisher, A., Flood, N. and Danaher, T. (2016) Comparing Landsat water index methods for automated water classification in eastern Australia. *Remote Sensing of Environment*, 175, 167-177.
-   Heimhuber, V. et al. (2018) Modelling 25 years of spatio-temporal surface water and inundation dynamics on large river basin scale using time series of Earth observation data. *Remote Sensing of Environment*, 209, 27-42.
-   Ji, L., Zhang, L. and Wylie, B. (2009) Analysis of dynamic thresholds for the normalized difference water index. *Photogrammetric Engineering and Remote Sensing*, 75(11), 1307-1317.
-   Ji, L. et al. (2015) Estimating aquatic vegetation fraction in a desert lake. *International Journal of Applied Earth Observation and Geoinformation*, 41, 109-117.
-   Lacaux, J.P. et al. (2007) Classification of ponds from high-spatial resolution remote sensing: application to Rift Valley Fever epidemics in Senegal. *Remote Sensing of Environment*, 106, 66-80.
-   Li, J. et al. (2023) Comparing water indices for Landsat data for automated surface water body extraction under complex ground background. *Remote Sensing*, 15(6), 1678.
-   Li, W. et al. (2013) A comparison of land surface water mapping using the normalized difference water index from TM, ETM+ and ALI. *Remote Sensing*, 5(11), 5530-5549.
-   Mahdavi, S. et al. (2018) Remote sensing for wetland classification: a comprehensive review. *Remote Sensing of Environment*, 206, 1-21.
-   Malahlela, O. (2016) Inland waterbody mapping: towards improving discrimination and extraction of inland surface water features. *Applied Geography*, 73, 77-88.
-   McFeeters, S.K. (1996) The use of the normalized difference water index (NDWI) in the delineation of open water features. *International Journal of Remote Sensing*, 17(7), 1425-1432.
-   Ogilvie, A. et al. (2015) Surface water monitoring in small water bodies: potential and limits of multi-sensor Landsat time series. *Remote Sensing*, 7, 16955-16982.
-   Oliphant, A.J. et al. (2019) Mapping cropland extent of Southeast and Northeast Asia using multi-year time-series Landsat 30-m data. *Remote Sensing of Environment*, 228, 1-13.
-   Ozesmi, S.L. and Bauer, M.E. (2002) Satellite remote sensing of wetlands. *Wetlands Ecology and Management*, 10, 381-402.
-   Pekel, J.-F. et al. (2016) High-resolution mapping of global surface water and its long-term changes. *Nature*, 540, 418-422.
-   Shen, L. and Li, C. (2010) Water body extraction from Landsat ETM+ imagery using AdaBoost algorithm. *Proceedings of the 18th International Conference on Geoinformatics*, Beijing, 1-4.
-   Xu, H. (2006) Modification of normalised difference water index (NDWI) to enhance open water features in remotely sensed imagery. *International Journal of Remote Sensing*, 27(14), 3025-3033.
