# Wetland Inundation Mapping through Remote Sensing and Ethnographic Validation: A Socio-Ecological Study of Lake Chilwa Basin

Seamus Murphy ^a\*^, John Wilson ^b^

^a\*^ Corresponding author: [seamusrobertmurphy\@gmail.com](mailto:seamusrobertmurphy@gmail.com){.email} ^b^ Community Natural Resource Management Consultant

**Abstract:**

Remote sensing of wetland landscapes in African conservation areas typically proceeds without reference to the people who inhabit them. This study integrates multi-sensor remote sensing with ethnographic fieldwork to map inundation dynamics in the endorheic Lake Chilwa Basin, Malawi. We evaluate five water-extraction indices (NDWI, MNDWI, AWEIsh, WRI, NDPI) derived from a multi-decadal Landsat time series and combine these with Sentinel-1 C-band SAR backscatter analysis to characterise recession-refilling cycles that occur at roughly 20-year intervals. Spectral mixture analysis enables sub-pixel estimation of water fraction across the shallow, turbid, and vegetated littoral zones where standard optical thresholds fail. Ethnographic data collected between 2012 and 2014 across lakeshore villages in Zomba, Phalombe, and Machinga districts provide ground validation and reveal socio-ecological patterns invisible to satellite analysis: the spatial structure of fishing regulations, seasonal migration timing, and enforcement conflicts across distinct territorial jurisdictions. Results show significant spatiotemporal variation in water extent, with historical recession events documented in 1879, 1900, 1914-15, 1922, 1931-32, 1934, 1954, 1960-61, 1967, 1973, 1995, and 2012. SAR integration addresses the specific failures of optical indices in this environment, particularly cloud-cover gaps during the wet season, spectral confusion in turbid shallow water, and inability to detect sub-canopy inundation beneath dense Typha marshes. The socio-ecological systems framework developed here demonstrates that effective conservation mapping requires both technical precision and local knowledge, and offers a transferable methodology for monitoring dynamic wetland ecosystems.

**Keywords:** wetland inundation mapping, spectral mixture analysis, SAR backscatter, participatory mapping, socio-ecological systems, endorheic lakes, fishery co-management, Lake Chilwa, Malawi

------------------------------------------------------------------------

## 1. Introduction

Endorheic watersheds in Africa present distinct challenges for remote sensing. Their shallow, turbid waters confound optical classification. Their margins shift seasonally and inter-annually across kilometres of flat terrain. Their ecological and economic significance is mediated by social systems that satellites cannot observe. Remote sensing provides spatial data essential to wetland monitoring (Eva and Lambin, 2000; Ozesmi and Bauer, 2002; Philippe and Karume, 2019), yet it captures neither the adaptive strategies of local resource users nor the political ecology that governs access to fluctuating resources (Yiran et al., 2012; Sulieman and Ahmed, 2013; Demichelis et al., 2020).

Technical limitations compound this problem. Inland waters differ fundamentally from the terrestrial surfaces for which most satellite sensors were designed. Coloured dissolved organic matter (CDOM), suspended sediment, and phytoplankton blooms alter the spectral properties of water in ways that atmospheric correction algorithms handle poorly (Matthews, 2011; Ogashawara et al., 2017). Higher-resolution sensors introduce their own constraints in spectral coverage, radiometric sensitivity, and temporal frequency (Palmer et al., 2015; Kutser, 2012). The challenge is not merely technical. It is that remote sensing alone produces maps that are technically defensible but ecologically and socially incomplete.

This study addresses both limitations simultaneously. We combine multi-sensor remote sensing with sustained ethnographic fieldwork to map inundation dynamics in the Lake Chilwa Basin, one of Africa's most productive and volatile endorheic systems. The remote sensing component evaluates optical water-extraction indices against SAR-derived water maps to quantify what each sensor family detects and what it misses. The ethnographic component, conducted over 18 months across three districts, validates the satellite-derived classifications and reveals the socio-ecological structures that determine how the lake's resources are governed, contested, and used.

### 1.1 The Problem with Optical Indices in Inland Wetlands

Water-extraction indices derived from multispectral imagery are the standard tools for mapping surface water extent. Each exploits a different combination of spectral bands to distinguish water from land, but each fails under specific conditions that Lake Chilwa routinely presents.

The Normalized Difference Water Index (NDWI; McFeeters, 1996) uses the ratio of green to near-infrared reflectance to delineate open water. It performs well where water is deep and clear (Ji et al., 2009) but struggles in shallow, turbid conditions because suspended sediment raises NIR reflectance, compressing index values toward zero (Xu, 2006; Li et al., 2013). In vegetated wetlands, canopy reflectance dominates both bands and the index cannot distinguish water beneath emergent vegetation from surrounding marshland (Ozesmi and Bauer, 2002). Endorheic lakes with fluctuating shorelines present mixed pixels at the land-water boundary where NDWI thresholds systematically misclassify wet soil as dry land (Deus and Gloaguen, 2013).

The Modified NDWI (MNDWI; Xu, 2006) substitutes shortwave infrared for NIR, suppressing the confusion between water and built-up or bare-soil surfaces that afflicts NDWI. It performs better in turbid and sediment-laden water (Li et al., 2013; Acharya et al., 2018) and is the strongest single optical index for the conditions this study addresses. It still fails beneath emergent vegetation, where macrophyte reflectance drives values negative even where standing water persists (Amani et al., 2020). In saline and endorheic systems, dissolved salts and algal blooms shift MNDWI thresholds unpredictably, requiring adaptive or Otsu-based methods (Pekel et al., 2016). Heimhuber et al. (2018) found that MNDWI required SAR fusion to map water extent beneath vegetation in the Murray-Darling system.

The Automated Water Extraction Index (AWEIsh; Feyisa et al., 2014) combines five bands to suppress shadow and dark-surface commission errors. Its multi-band formulation captures contrasts that two-band indices cannot resolve, but subpixel mixing of emergent vegetation and water still reduces values below detection thresholds (Ji et al., 2015). Fisher et al. (2016) found that no single index, AWEIsh included, consistently mapped shallow or seasonal water bodies in a global comparison. In endorheic systems, AWEIsh can misclassify exposed lakebed and brine as water (Pekel et al., 2016).

The Water Ratio Index (WRI; Shen and Li, 2010) uses the ratio of visible to infrared bands, yielding values above 1.0 for water. It suppresses cloud and shadow noise effectively (Li et al., 2023) and has shown strong performance in bare-soil landscapes in Ethiopia (Demelash et al., 2025). But where water, vegetation, and sediment share a pixel, the visible-band numerator inflates from non-water reflectance. The red band is particularly vulnerable to suspended sediment, which raises its reflectance sharply in turbid conditions. WRI has not been applied to endorheic or fluctuating lake systems.

The Normalized Difference Pond Index (NDPI; Lacaux et al., 2007) was designed in and for semi-arid Africa, using SPOT-5 imagery to classify temporary ponds in Senegal's Ferlo region. It exploits the SWIR absorption of water against green vegetation reflectance, making it sensitive to the water-vegetation boundary where other indices fail. Campos et al. (2012) confirmed that SWIR-based indices outperform NIR-based ones for ephemeral ponds in the Sahel. NDPI is the algebraic inverse of MNDWI, carrying identical spectral information with reversed sign convention. It performs best on small, shallow, vegetated water bodies and loses sensitivity in deep or highly turbid water.

The pattern across all five indices is consistent: each occupies a spectral niche but fails at the boundaries that define Lake Chilwa's hydrology. Shallow turbid water, vegetated wetland margins, fluctuating shorelines, and mixed pixels at the land-water transition are precisely the conditions where optical classification underperforms. These are also the zones of greatest ecological and socio-economic significance.

### 1.2 SAR as a Complement to Optical Methods

Synthetic aperture radar addresses the specific failures of optical indices in inland wetlands. SAR operates at microwave frequencies, independent of cloud cover and solar illumination, two properties that make it indispensable for monitoring tropical wetlands where persistent cloud obscures optical sensors for months during the critical wet season (Mahdavi et al., 2018).

Smooth open water produces near-specular reflection, returning very little backscatter to the sensor. C-band SAR typically records open water at -20 to -30 dB, well below surrounding land surfaces, enabling detection through simple thresholding (Martinis et al., 2015). Where vegetation stands in water, the signal follows a double-bounce pathway, reflecting from the water surface and then from vertical plant structures. This produces backscatter paradoxically higher than from the same vegetation when dry, because the smooth water surface enhances the coherent specular return (Hess et al., 1995; Tsyganskaya et al., 2018). The contrast between double-bounce returns from flooded vegetation and volume scattering from dry canopies is the physical basis for detecting sub-canopy inundation.

Sentinel-1, operational since 2014, provides C-band imagery at 10 m resolution with 6 to 12 day revisit. Its dual-polarisation mode (VV/VH) has been applied to wetland inundation mapping across diverse environments: the St. Lucia wetlands in South Africa (Clement et al., 2018), the Amazon lowlands (Hardy et al., 2019), the Okavango Delta (Luca et al., 2025), and seasonal lake cycles on the Tibetan Plateau (Zhang et al., 2019). Multi-temporal approaches exploit wet-dry season contrast to characterise hydroperiod without ground data. The Sentinel-1 Global Flood Monitoring service now processes all incoming acquisitions for near-real-time detection (Roth et al., 2025).

C-band has limitations. Speckle noise degrades classification accuracy. Wind roughens water surfaces, raising backscatter above detection thresholds and causing false negatives (Martinis et al., 2015). Most critically, C-band's 5.6 cm wavelength cannot penetrate dense emergent vegetation such as Typha and Phragmites. Clement et al. (2018) reported poor classification in heavily vegetated zones. L-band SAR (approximately 23 cm wavelength), carried by ALOS PALSAR and the forthcoming NISAR mission, penetrates canopies far more effectively. Hess et al. (2003) showed that L-band detected sub-canopy flooding at 50 cm depth in Amazonian floodplains, versus 80 cm for C-band.

SAR interferometry offers additional capability. InSAR measures phase differences between acquisitions to detect centimetre-scale vertical displacement. In wetlands, the double-bounce mechanism preserves coherence, enabling water level change detection at resolutions exceeding gauge networks (Wdowinski et al., 2008; Hong and Wdowinski, 2017). Kim et al. (2021) mapped water level gradients across Tonle Sap at sub-monthly intervals. Temporal decorrelation and atmospheric phase delays limit the technique, particularly at C-band over intervals exceeding two weeks (Zebker and Villasenor, 1992). InSAR water level retrieval has not been attempted for African endorheic systems.

Fusing SAR and optical data exploits their complementarity. Pixel-level stacking, decision-level fusion, and machine learning classifiers trained on combined feature spaces consistently outperform single-source models, with reported accuracies of 85 to 95 percent for wetland classes (Amani et al., 2019; Whyte et al., 2018). Xu et al. (2025) mapped surface water dynamics across East Africa at 10 m resolution by integrating Sentinel-1 and Sentinel-2 time series. Lubala et al. (2023) fused Sentinel-1, Sentinel-2, and ALOS PALSAR to map small inland wetlands in the Democratic Republic of Congo. These studies confirm that multi-sensor integration captures dynamics that either sensor family misses alone.

### 1.3 Study Area

Lake Chilwa Basin is one of Africa's most productive endorheic ecosystems. Located in southern Malawi, this shallow terminal basin supports one of the continent's most densely populated regions (500 km^-1^) and exhibits remarkable ecological productivity during wet periods (Njaya et al., 2011).

The basin's hydrology is driven by unimodal rainfall (November to April) and sporadic chiperone rains (May to August). Annual fluctuations follow precipitation patterns closely (Ngongondo et al., 2011; Nicholson et al., 2014). Longer-term cycles of approximately 15 to 20 years produce dramatic lake recessions and varying degrees of complete desiccation (Kambombe et al., 2021). Historical records document major recessions in 1879, 1900, 1914-15, 1922, 1931-32, 1934, 1954, 1960-61, 1967, 1973, 1995, and 2012.

During recession, aquatic species take refuge in residual swamps dominated by salt-hardy Typha domingensis Pers. (Howard-Williams and Walker, 1974; Howard-Williams and Gaudet, 1985). Refilling initiates complex succession dynamics, with emergent food webs driven by detritus and bacterial processes in alkaline, nutrient-rich sediments (Kalk and Schulten-Senden, 1977; Furse et al., 1979). Peak productivity reaches 159 kg ha^-1^ (1979) and 113 kg ha^-1^ (1990), surpassing Lake Malawi (40 kg ha^-1^), Lake Tanganyika (90 kg ha^-1^), and Lake Victoria (116 kg ha^-1^) (Njaya et al., 2011; Vanden Bossche and Bernacsek, 1990). This boom-and-bust pattern sustains permanent lakeshore communities and seasonal migrant fishing populations.

### 1.4 Research Objectives

This study presents an integrated remote sensing and participatory mapping approach to characterise lacustrine transgression-regression dynamics and associated socio-ecological patterns in the Lake Chilwa Basin. The objectives are:

1.  Evaluate five optical water-extraction indices (NDWI, MNDWI, AWEIsh, WRI, NDPI) and Sentinel-1 C-band SAR for mapping inundation extent across the basin's full range of conditions, from deep open water through shallow turbid zones to vegetated marshland.

2.  Integrate SAR and optical time series using spectral mixture analysis to quantify hydroperiod fluctuations and littoral zone dynamics across multiple recession-refilling cycles.

3.  Validate satellite-derived classifications against ethnographic data collected through key informant interviews, focus group discussions, and participatory mapping with migrant fishing communities.

4.  Document the spatiotemporal structure of fishing regulations, population movements, and resource use patterns that formal monitoring frameworks have not captured.

------------------------------------------------------------------------

## 2. Methods

### 2.1 Socio-Ecological Systems Framework

The study employs a two-stage socio-ecological systems (SES) framework integrating biophysical remote sensing with social research methods. The framework targets the perspectives of marginalised groups, particularly migrant fishers, who form a principal economic segment of the Lake Chilwa system yet are overlooked in conventional management frameworks.

Data collection occurred between September 2012 and March 2014 across lakeshore villages in Zomba, Phalombe, and Machinga districts within the Lake Chilwa Ramsar zones. Study site boundaries were defined through participatory mapping workshops with multi-stakeholder groups and Department of Fisheries officers.

**Qualitative data collection** comprised four methods. Key informant interviews (n=45) were semi-structured conversations with village leaders, fishing camp chairmen, Department of Fisheries officers, and long-term residents, focused on historical lake dynamics, fishing regulations, and seasonal migration. Focus group discussions (n=18) convened separate sessions with migrant fishers, women fish processors, boat owners, and service providers to capture competing perspectives on resource access, livelihood strategies, and enforcement disputes. Participatory rural appraisals employed seasonal calendars, resource mapping, and historical timelines to document collective knowledge of lake dynamics and fishery management. Extended participatory observation in fishing camps documented daily practices, operational logistics, social networks, and adaptive strategies during different hydrological phases.

**Geographic data collection** used differential GPS receivers to record landscape structure and lakeshore dynamics during field visits. Participatory workshops enabled community identification of fishing infrastructure (permanent and seasonal camps, landing sites, processing areas), ecological zones (wetland boundaries, vegetation transitions, spawning areas), cultural landscapes (sacred sites, traditional fishing territories, conflict zones), and seasonal patterns (water level indicators, migration routes, market locations). Community knowledge was integrated with remote sensing through iterative validation workshops where preliminary satellite-derived maps were ground-truthed against local observations. This process revealed discrepancies between technical classifications and actual resource use patterns, leading to refined mapping approaches.

### 2.2 Remote Sensing Framework

The remote sensing workflow integrates multi-temporal SAR backscatter analysis with optical spectral indices to map surface water extent variability. The temporal window extends as far back as Landsat data quality permits, driven by the need to capture multiple recession events that recur at approximately 20-year intervals.

#### SAR Processing

SAR processing exploits the sensitivity of C-band radar to backscatter differences between smooth water surfaces and rough terrestrial features. Calm water yields low backscatter (-20 to -30 dB) while vegetated areas exhibit higher returns from volume scattering and surface roughness interactions. Wet soils produce higher backscatter than dry soils due to their increased dielectric constant. VV polarisation provides greatest sensitivity to soil moisture; cross-polarisation (VH) differentiates woody from herbaceous vegetation (Tsyganskaya et al., 2018).

Sentinel-1 data were processed using a standardised workflow:

1.  Radiometric calibration: conversion of digital numbers to sigma-0 backscatter coefficients.
2.  Speckle filtering: Lee Sigma filter (7 x 7 window) to reduce multiplicative noise.
3.  Geometric correction: Range-Doppler terrain correction using SRTM 30 m DEM.
4.  Multi-temporal coregistration: sub-pixel alignment of wet/dry season image pairs.
5.  Conversion to dB scale for threshold-based water detection.

Multi-temporal gradient analysis, adapted from sea ice monitoring methodologies, enhanced detection of dynamic water boundaries through comparative analysis of seasonal backscatter patterns. This approach proved effective for identifying transitions between open water, flooded vegetation, and terrestrial surfaces.

#### Landsat Processing

Optical analysis used Analysis Ready Data products from Landsat Collection 2, specifically Level-2 surface reflectance from the Thematic Mapper (L5-TM), Enhanced Thematic Mapper Plus (L7-ETM+), and Operational Land Imager (L8-OLI). Earlier sensors (Landsat 3 and 4) were evaluated but present gaps, cloud interference, sensor degradation, and archival quality issues that reduce usable coverage, particularly before 1984. The temporal window is constrained by these data quality limitations rather than by methodological choice.

To address atmospheric overcorrection and geometric artefacts in aquatic pixels, raster pixels were scanned for open water voids using an inverted elastic-net, treated with dark spectrum fitting, then sorted into spectral clusters and assessed for outliers using a time-weighted hierarchical clustering function.

#### Spectral Index Variables

Five water-extraction indices were evaluated, selected to span the range of spectral approaches available for inland water mapping and to test their relative performance under the specific conditions Lake Chilwa presents:

**NDWI** (McFeeters, 1996): (Green - NIR) / (Green + NIR). The original water index, sensitive to deep open water but prone to false negatives in turbid, shallow, or vegetated conditions due to elevated NIR reflectance from suspended sediment and canopy effects. Serves as the baseline against which other indices are compared.

**MNDWI** (Xu, 2006): (Green - SWIR1) / (Green + SWIR1). Substitutes SWIR for NIR, improving discrimination of water from built-up and bare-soil surfaces. The strongest single optical index for turbid water detection, though it fails beneath emergent vegetation and requires adaptive thresholding in saline systems.

**AWEIsh** (Feyisa et al., 2014): Blue + 2.5 x Green - 1.5 x (NIR + SWIR1) - 0.25 x SWIR2. A five-band combination optimised for shadow suppression and complex landscape discrimination. Outperforms simpler indices where topographic or shadow effects confound classification but offers no clear advantage in shallow vegetated wetlands.

**WRI** (Shen and Li, 2010): (Green + Red) / (NIR + SWIR1). A ratio index that suppresses cloud and shadow noise effectively. Competitive in bare-soil landscapes but vulnerable to inflation from suspended sediment in the red band and untested in endorheic systems.

**NDPI** (Lacaux et al., 2007): (SWIR1 - Green) / (SWIR1 + Green). Designed for temporary pond detection in semi-arid Africa using SPOT-5 data. The algebraic inverse of MNDWI, sensitive to the water-vegetation boundary where other indices fail. Most ecologically appropriate for Lake Chilwa's seasonal vegetated margins but loses sensitivity in deep or turbid water.

The multi-index approach tests each against SAR-derived water maps to quantify what optical sensors detect and what they miss across the basin's full range of conditions.

#### Spectral Mixture Analysis

Spectral mixture analysis enabled sub-pixel water fraction estimation, critical for monitoring gradual transitions between terrestrial and aquatic habitats. The approach was selected over object-based methods based on demonstrated superior performance in delineating turbid waters, shallow wetlands, and mixed vegetation-water pixels in lakeshore environments (Halabisky et al., 2016; Huang et al., 2014; Shanmugam et al., 2006). The spectral characteristics of Lake Chilwa, with its dense marshlands, shallow waters, extensive detritus, phytoplankton blooms, and shoreline shadowing, make sub-pixel estimation essential.

Endmember selection followed standard protocols, with training samples drawn from spectrally pure pixels identified through iterative refinement and community validation:

-   Open water: 150 samples, purity threshold \>95%
-   Flooded vegetation: 200 samples, purity threshold \>85%
-   Dry vegetation: 180 samples, purity threshold \>90%
-   Bare soil: 120 samples, purity threshold \>85%
-   Urban/built: 80 samples, purity threshold \>90%

#### Training Samples and Classification

Training data collection integrated remote sensing requirements with community knowledge validation. Participatory workshops enabled local experts to identify spectrally similar but functionally different landscape units, such as seasonal versus permanent wetlands and distinct fishing zones, that satellite imagery alone could not distinguish. Random Forest classification was implemented with 500 trees and temporal weighting to account for seasonal variations in water extent and vegetation phenology.

#### Accuracy Assessment

Classification accuracy was assessed through both conventional metrics and community validation. Overall accuracy reached 81% with a kappa coefficient of 0.77. Producer's accuracy ranged from 71% (bare soil) to 89% (open water); user's accuracy from 68% (bare soil) to 92% (open water). Community validation showed 74 to 92% agreement on boundary delineation and seasonal timing, with disagreements concentrated in mixed-pixel zones and areas of rapid temporal change.

------------------------------------------------------------------------

## 3. Results

[Results section to be expanded with analysis outputs.]

The analysis revealed significant spatiotemporal variation in water extent across the study period. Maximum wet-season extent ranged from 1,500 to 2,200 km^2^ (January to March); minimum dry-season extent from 400 to 800 km^2^ (August to October). The inter-annual coefficient of variation was 0.42. Both the 1995 and 2012 recession events were detected within the observation period.

The integrated approach revealed socio-ecological patterns invisible to technical analysis alone. Migration timing correlated strongly with satellite-detected water level changes (r=0.73, p\<0.001) with a 3-month anticipatory lag. Fishing regulation compliance reached 65% during active enforcement, with selective compliance based on traditional territorial boundaries. Elder mediation resolved 89% of resource disputes where traditional boundaries received formal recognition. Market participation showed distance decay (beta=-0.45), with road quality as the critical determinant.

Community mapping identified specific locations where formal fisheries regulations conflicted with traditional practices. Recognition of traditional territorial boundaries proved essential for improving compliance during spawning periods, challenging the assumption that formal and traditional management systems are inherently conflicting.

------------------------------------------------------------------------

## 4. Discussion

[Discussion section retained from previous draft, to be revised with deeper literature grounding in subsequent iterations.]

The integration of remote sensing and participatory methods revealed critical challenges and demonstrated solutions for wetland conservation research. The primary methodological challenge involved reconciling temporal mismatches between 16-day Landsat revisit cycles and daily community observations. This required seasonal aggregation methods that preserved both satellite data precision and local knowledge temporality. Thirty-metre Landsat pixels proved inadequate for capturing fishers' spatial knowledge of specific fishing grounds, necessitating sub-pixel analysis validated through community mapping.

The integrated approach uncovered socio-ecological patterns invisible to technical analysis alone. Community mapping identified locations where formal fisheries regulations conflicted with traditional practices, enabling policy adjustments that improved conservation outcomes and compliance. Local observers consistently identified environmental shifts days or weeks before satellite detection, demonstrating the potential for collaborative monitoring networks.

This study demonstrates that purely technical remote sensing approaches miss the social dimensions that determine conservation success. The socio-ecological systems framework shows transferability through standardised remote sensing protocols, replicable participatory methods, and documented integration workflows. Future research should extend temporal coverage using earlier Landsat sensors where data quality permits, integrate L-band SAR for improved vegetation penetration, develop coupled hydrological-social models for scenario planning, and test transferability across other African endorheic systems.

------------------------------------------------------------------------

## 5. Conclusions

[To be developed.]

------------------------------------------------------------------------

## 6. References

[References to be compiled from both the existing references.bib and the new citations introduced in this draft. New citations requiring addition to the bibliography are listed below.]

### New citations introduced in this draft (not yet confirmed in references.bib):

-   Acharya, T.D., Subedi, A. and Lee, D.H. (2018) Sensors, 18(8), 2580.
-   Amani, M. et al. (2019) GIScience and Remote Sensing, 56(7), 1023-1045.
-   Amani, M. et al. (2020) Remote Sensing, 12, 1190.
-   Campos, J.C., Sillero, N. and Brito, J.C. (2012) Journal of Hydrology, 464-465, 438-446.
-   Clement, M.A. et al. (2018) Journal of Flood Risk Management, 11(2), 152-168.
-   Demelash, T. et al. (2025) H2Open Journal, 8(5), 402 ff.
-   Deus, D. and Gloaguen, R. (2013) Remote Sensing, 5, 5765-5781.
-   Feyisa, G.L. et al. (2014) Remote Sensing of Environment, 140, 23-35.
-   Fisher, A., Flood, N. and Danaher, T. (2016) Remote Sensing of Environment, 175, 167-177.
-   Halabisky, M. et al. (2016) Remote Sensing of Environment, 177, 171-183.
-   Hardy, A. et al. (2019) Remote Sensing, 11(6), 720.
-   Heimhuber, V. et al. (2018) Remote Sensing of Environment, 209, 27-42.
-   Hess, L.L. et al. (1995) IEEE Trans. Geoscience and Remote Sensing, 33(4), 896-904.
-   Hess, L.L. et al. (2003) Remote Sensing of Environment, 87, 404-428.
-   Hong, S.-H. and Wdowinski, S. (2017) IEEE Geoscience and Remote Sensing Letters, 14(8), 1370-1374.
-   Huang, C. et al. (2014) Remote Sensing of Environment, 141, 231-242.
-   Ji, L., Zhang, L. and Wylie, B. (2009) Photogrammetric Engineering and Remote Sensing, 75(11), 1307-1317.
-   Ji, L. et al. (2015) Int. J. Applied Earth Observation and Geoinformation, 41, 109-117.
-   Kim, S. et al. (2021) Remote Sensing of Environment, 252, 112150.
-   Lacaux, J.P. et al. (2007) Remote Sensing of Environment, 106, 66-80.
-   Li, J. et al. (2023) Remote Sensing, 15(6), 1678.
-   Li, W. et al. (2013) Remote Sensing, 5(11), 5530-5549.
-   Lubala, N. et al. (2023) Ecological Indicators, 155, 110962.
-   Luca, G. et al. (2025) Remote Sensing of Environment, 310, 114474.
-   Mahdavi, S. et al. (2018) Remote Sensing of Environment, 206, 1-21.
-   Martinis, S. et al. (2015) Natural Hazards and Earth System Sciences, 15(5), 1069-1085.
-   McFeeters, S.K. (1996) Int. J. Remote Sensing, 17(7), 1425-1432.
-   Pekel, J.-F. et al. (2016) Nature, 540, 418-422.
-   Roth, F. et al. (2025) Remote Sensing of Environment, 310, 114493.
-   Shanmugam, P., Ahn, Y.-H. and Sanjeevi, S. (2006) Ecological Modelling, 194(4), 379-394.
-   Shen, L. and Li, C. (2010) Proceedings 18th Int. Conf. Geoinformatics, Beijing, 1-4.
-   Tsyganskaya, V. et al. (2018) Int. J. Applied Earth Observation and Geoinformation, 73, 205-218.
-   Wdowinski, S. et al. (2008) Remote Sensing of Environment, 112, 681-696.
-   Whyte, A. et al. (2018) Environmental Modelling and Software, 104, 40-54.
-   Xu, H. (2006) Int. J. Remote Sensing, 27(14), 3025-3033.
-   Xu, Y. et al. (2025) ISPRS J. Photogrammetry and Remote Sensing, 215, 101-118.
-   Zebker, H.A. and Villasenor, J. (1992) IEEE Trans. Geoscience and Remote Sensing, 30(5), 950-959.
-   Zhang, G. et al. (2019) Science of the Total Environment, 655, 1390-1401.
