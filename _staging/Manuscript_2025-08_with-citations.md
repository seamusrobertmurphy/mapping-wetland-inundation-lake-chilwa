# Mapping inundation dynamics and fisher migration across Lake Chilwa's recession cycles: Integrating SAR-optical remote sensing with ethnographic fieldwork in an endorheic wetland

Seamus Murphy ^a\*^, John Wilson ^b^

^a\*^ Corresponding author: [seamusrobertmurphy\@gmail.com](mailto:seamusrobertmurphy@gmail.com){.email} ^b^ Community Natural Resource Management Consultant

**Abstract:**

Remote sensing of wetland landscapes in African conservation areas typically proceeds without reference to the people who inhabit them. This study integrates multi-sensor remote sensing with ethnographic fieldwork to map inundation dynamics in the endorheic Lake Chilwa Basin, Malawi. We evaluate five water-extraction indices (NDWI, MNDWI, AWEIsh, WRI, NDPI) derived from a multi-decadal Landsat time series and combine these with Sentinel-1 C-band SAR backscatter analysis to characterise recession-refilling cycles that occur at roughly 20-year intervals. Spectral mixture analysis enables sub-pixel estimation of water fraction across the shallow, turbid, and vegetated littoral zones where standard optical thresholds fail. Ethnographic data collected between 2012 and 2014 across lakeshore villages in Zomba, Phalombe, and Machinga districts provide ground validation and reveal socio-ecological patterns invisible to satellite analysis: the spatial structure of fishing regulations, seasonal migration timing, and enforcement conflicts across distinct territorial jurisdictions. Results show significant spatiotemporal variation in water extent, with historical recession events documented in 1879, 1900, 1914-15, 1922, 1931-32, 1934, 1954, 1960-61, 1967, 1973, 1995, and 2012. SAR integration addresses the specific failures of optical indices in this environment, particularly cloud-cover gaps during the wet season, spectral confusion in turbid shallow water, and inability to detect sub-canopy inundation beneath dense Typha marshes. The socio-ecological systems framework developed here demonstrates that effective conservation mapping requires both technical precision and local knowledge, and offers a transferable methodology for monitoring dynamic wetland ecosystems.

**Keywords:** wetland inundation mapping, spectral mixture analysis, SAR backscatter, participatory mapping, socio-ecological systems, endorheic lakes, fishery co-management, Lake Chilwa, Malawi

------------------------------------------------------------------------

## 1. Introduction

Wetlands cover between 3% and 8% of the Earth's land surface and provide ecosystem services disproportionate to their extent: flood attenuation, water quality improvement, carbon sequestration, and biodiversity support (Mitsch and Gosselink, 2015; Davidson, 2014). Despite this, most countries lack comprehensive wetland inventories. Existing maps are fragmented, produced at incompatible scales with inconsistent methods, and rarely updated to reflect the dynamic character of wetland systems (Finlayson et al., 2018; Mahdianpari et al., 2019). In sub-Saharan Africa, the problem is acute. Many wetlands of continental significance have no systematic inventory at all, and those that exist rely on single-date optical imagery that captures a snapshot rather than the seasonal and inter-annual variability that defines wetland function (Rebelo et al., 2009; Muro et al., 2016).

Remote sensing offers the spatial and temporal coverage that field surveys cannot match, and two decades of methodological development have produced a substantial toolkit for wetland mapping (Ozesmi and Bauer, 2002; Mahdavi et al., 2018). Yet accurate mapping of inland wetlands remains difficult. Shallow, turbid waters confound optical classification. Cloud cover obscures tropical wet seasons precisely when inundation is greatest. Emergent vegetation masks the water beneath it. Coloured dissolved organic matter, suspended sediment, and phytoplankton blooms alter spectral properties in ways that atmospheric correction algorithms handle poorly (Matthews, 2011; Ogashawara et al., 2017). Higher-resolution sensors introduce their own constraints in spectral coverage, radiometric sensitivity, and temporal frequency (Palmer et al., 2015; Kutser, 2012).

The concurrent availability of Sentinel-1 SAR and Sentinel-2 optical data since 2014, combined with the multi-decadal Landsat archive made freely available since 2008, now provides an unprecedented opportunity to address these limitations through multi-sensor time series analysis. Cloud computing platforms, particularly Google Earth Engine (GEE), have removed the computational barriers that previously made processing such large volumes of satellite imagery infeasible on conventional systems (Gorelick et al., 2017; Mahdianpari et al., 2019). These advances enable on-demand, large-scale wetland mapping that was not possible a decade ago.

Endorheic watersheds in Africa present a further challenge that is not merely technical but social. Their ecological and economic significance is mediated by human populations whose livelihoods track the water. Remote sensing provides essential spatial data (Eva and Lambin, 2000; Philippe and Karume, 2019), yet it captures neither the adaptive strategies of local resource users nor the political ecology that governs access to fluctuating resources (Yiran et al., 2012; Sulieman and Ahmed, 2013; Demichelis et al., 2023). Maps that are technically defensible remain ecologically and socially incomplete.

This study addresses both limitations simultaneously. We combine multi-sensor remote sensing with sustained ethnographic fieldwork to map inundation dynamics in the Lake Chilwa Basin, one of Africa's most productive and volatile endorheic systems. The remote sensing component evaluates optical water-extraction indices against SAR-derived water maps, processed through a Google Earth Engine pipeline, to quantify what each sensor family detects and what it misses. The ethnographic component, conducted over 18 months across three districts, validates the satellite-derived classifications and reveals the socio-ecological structures that determine how the lake's resources are governed, contested, and used.

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

Lake Chilwa Basin is one of Africa's most productive endorheic ecosystems. Located in southern Malawi, this shallow terminal basin spans approximately 2,310 km^2^, of which only one third is open water, the remainder comprising Typha swamp and seasonally inundated marshland (Kalk, 1979). The lake was designated a Ramsar Wetland of International Importance in 1997, in recognition of its outstanding waterfowl populations and ecological significance (Wilson, 2007). The three surrounding districts, Zomba, Phalombe, and Machinga, support population densities of 162 persons km^-2^ around the lake itself and considerably higher in the southern basin (Wilson, 2010), with 70 to 81% of households living below the poverty line and farming plots averaging 0.35 ha (Wilson, 2010). The fishery contributes on average 20% of Malawi's annual catch and in peak years has reached 43% (Chiotha, 1996; Njaya et al., 2011).

The basin's hydrology is driven by unimodal rainfall (November to April) and sporadic chiperone rains (May to August). Annual fluctuations follow precipitation patterns closely (Ngongondo et al., 2011; Nicholson et al., 2014). Longer-term cycles of approximately 15 to 20 years produce dramatic lake recessions and varying degrees of complete desiccation (Kambombe et al., 2021). Wilson (2014) compiled the most complete chronology of lake levels, documenting major recessions in 1900, 1913-16, 1922, 1934, 1943, 1946, 1952, 1960, 1967-68, 1995, and 2012-13, set within a geological record of pluvial and inter-pluvial phases extending 450,000 years. Average maximum lake depth is 2.95 m; a prolonged dry phase between approximately 1760 and 1850 deposited a hard layer of homogeneous sandy clays at least 1.2 m thick that now underlies the lakebed (Wilson, 2014; Crossley et al., 1983). These fluctuations are not anomalies but the lake's defining ecological characteristic: a transient world that in some seasons supports a thriving fishery, a bustling marketing network, and burgeoning lakeshore settlements, while in the following season the shoreline may retreat by fifteen kilometres, leaving communities stranded from open water (Kalk, 1979; Allison and Mvula, 2002).

During recession, the lake's three main commercial species, Barbus paludinosus (matemba), Clarias gariepinus (mlamba), and Oreochromis shiranus chilwae (chambo), take refuge in residual deep pools at alluvial river mouths and in swamps dominated by salt-hardy Typha domingensis Pers. (Howard-Williams and Walker, 1974; Howard-Williams and Gaudet, 1985). These species have evolved characteristics suited to a fluctuating environment: high fecundity, early reproductive maturity, broad diets, unspecialised spawning habits, and tolerance to wide environmental variation (Njaya, 2001). Refilling initiates complex succession dynamics, with emergent food webs driven by detritus and bacterial processes in alkaline, nutrient-rich sediments (Kalk and Schulten-Senden, 1977; Furse et al., 1979). Peak productivity reaches 159 kg ha^-1^ (1979) and 113 kg ha^-1^ (1990), surpassing Lake Malawi (40 kg ha^-1^), Lake Tanganyika (90 kg ha^-1^), and Lake Victoria (116 kg ha^-1^) (Njaya et al., 2011; Vanden Bossche and Bernacsek, 1990). Wilson (2014) notes that the 1979 bumper harvest of 24,310 tonnes followed six years after the lake dried, during which oxidation of organic matter released nutrients that fuelled explosive productivity upon refilling.

This boom-and-bust ecology sustains a complex mobile population. The basin's settlement history stretches back millennia, from Akafula hunter-gatherers through the 16th-century Maravi (Mang'anja) chieftancies, the mid-19th-century Yao migration from Mozambique, and the 20th-century Lomwe influx from Portuguese East Africa (Wilson, 2010). Today, permanent lakeshore communities of Yao, Nyanja, and Lomwe matrilineal descent groups farm the seasonally exposed lakebed and fish inshore waters. Seasonal migrants travel from as far as Machinga district to the north, spending weeks to months on zimbowera, floating platforms of piled Typha grass that serve as fishing camps in the lake interior (Murphy, 2014). These zimbowera neighbourhoods, such as Andere, Chambwalu, and Lingoni, support shifting populations of 300 to 1,000 fishers and sustain their own tea rooms, trading posts, and social governance structures. During the wet season, when fishing grounds are most productive but most inaccessible from the mainland, zimbowera cooperatives form along lines of village membership, kinship, and occupational specialisation. Seine-net crews of up to 80 fishers operate from large camps under single owners, while smaller gill-net cooperatives of three to four family members share boats, equipment, and risks on remote islands (Murphy, 2014). The fishery's principal port is Kachulu, on the western shore in Zomba district, which serves as the hub for a regional marketing network extending by bicycle and truck to Zomba, Blantyre, Liwonde, Phalombe, and Mulanje.

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

Data collection occurred between September 2012 and March 2014 across lakeshore villages in Zomba, Phalombe, and Machinga districts within the Lake Chilwa Ramsar zones. Study site boundaries were defined through participatory mapping workshops with multi-stakeholder groups and Department of Fisheries officers. The fisheries governance structure within which fieldwork operated, comprising six Fisheries Associations aligned to Traditional Authorities and 53 Beach Village Committees aligned to Group Village Headmen, is described by Wilson (2009) in the Lake Chilwa and Mpoto Lagoon Fisheries Management Plan.

**Qualitative data collection** comprised four methods. Key informant interviews (n=45) were semi-structured conversations with village leaders, fishing camp chairmen, Beach Village Committee (BVC) members, Department of Fisheries officers, and long-term residents, focused on historical lake dynamics, fishing regulations, and seasonal migration. Focus group discussions (n=18) convened separate sessions with migrant seine-net fishers, resident gill-net fishers, women fish processors, boat owners, bicycle traders, and zimbowera cooperative leaders to capture competing perspectives on resource access, livelihood strategies, and enforcement disputes. The research was based principally at Kachulu, the fishery's busiest port on the western shore, with extended visits to fishing camps at Napali, Andere, Chambwalu, Lingoni, and Manda Manjeza in the lake interior and along the northern shore in Machinga district. Participatory rural appraisals employed seasonal calendars, resource mapping, and historical timelines to document collective knowledge of lake dynamics and fishery management. Extended participatory observation in fishing camps and zimbowera neighbourhoods documented daily practices, operational logistics, social networks, and adaptive strategies during different hydrological phases (Murphy, 2014).

**Geographic data collection** used differential GPS receivers to record landscape structure and lakeshore dynamics during field visits. Participatory workshops enabled community identification of fishing infrastructure, including permanent and seasonal camps, landing sites, processing areas, and the canal systems (such as Mapila Canal, approximately 10 to 15 km depending on lake levels) through which fishers access otherwise inaccessible interior fishing grounds. Ecological zones including wetland boundaries, vegetation transitions, and spawning areas were mapped alongside cultural landscapes such as sacred fishing sites on Chisi Island (Kalanda-Sabola et al., 2007), traditional fishing territories, and conflict zones between district jurisdictions. Seasonal patterns including water level indicators, migration routes from Machinga to Zomba and Phalombe districts, and market network locations were recorded. Community knowledge was integrated with remote sensing through iterative validation workshops where preliminary satellite-derived maps were ground-truthed against local observations. This process revealed discrepancies between technical classifications and actual resource use patterns, leading to refined mapping approaches.

### 2.2 Reference Data

In-situ reference data were collected during the ethnographic fieldwork between September 2012 and March 2014. Differential GPS receivers recorded landscape features, water boundaries, and vegetation transitions during field visits across the three study districts. A total of 45 key informant sites and 18 focus group locations were georeferenced, along with 23 fishing camps (permanent and seasonal), landing sites, processing areas, and ecological transition zones identified through participatory mapping workshops.

GPS points were imported into a GIS environment and used to generate training and validation polygons through visual interpretation of high-resolution Google Earth imagery cross-referenced with field photographs and community annotations. Polygons were sorted by size and alternately assigned to training (approximately 50%) and testing (approximately 50%) groups to ensure independent validation samples with balanced representation of small and large features, following the protocol of Mahdianpari et al. (2019). Community validation during iterative feedback workshops identified spectrally ambiguous landscape units, such as seasonal versus permanent wetlands and distinct fishing zones, that required field verification to classify correctly.

### 2.3 Remote Sensing Framework

The remote sensing workflow integrates multi-temporal SAR backscatter analysis with optical spectral indices to map surface water extent variability. The area of interest was derived from HydroSHEDS level-6 basin boundaries filtered by a point location within Lake Chilwa, producing a single basin polygon used as the primary clip boundary for all raster and vector queries. Hydrological base layers from MERIT Hydro (v1.0.3), HydroSHEDS (v1.1), and JRC Global Surface Water (v1.4; 454 monthly images) provided elevation, flow direction, upstream drainage area, and historical surface water occurrence for cross-validation purposes. The temporal window extends as far back as Landsat data quality permits, driven by the need to capture multiple recession events that recur at approximately 20-year intervals.

#### SAR Processing

SAR processing exploits the sensitivity of C-band radar to backscatter differences between smooth water surfaces and rough terrestrial features. Calm water yields low backscatter (-20 to -30 dB) while vegetated areas exhibit higher returns from volume scattering and surface roughness interactions. Wet soils produce higher backscatter than dry soils due to their increased dielectric constant. VV polarisation provides greatest sensitivity to soil moisture; cross-polarisation (VH) differentiates woody from herbaceous vegetation (Tsyganskaya et al., 2018).

Sentinel-1 data were processed through a Google Earth Engine pipeline that handles radiometric calibration, speckle filtering, and geometric correction at scale. GEE provides access to the full Sentinel-1 Ground Range Detected archive, pre-processed to sigma-0 backscatter coefficients with thermal noise removal, radiometric calibration, and terrain correction applied (Gorelick et al., 2017). The cloud-based pipeline eliminates the need to download, store, and process individual scenes locally, enabling batch analysis of multi-temporal image stacks that would be computationally prohibitive on conventional workstation systems (Mahdianpari et al., 2019). We filtered the archive to Interferometric Wide Swath mode with VV and VH polarisations in descending orbit, yielding 913 scenes over the study area. The processing workflow comprised:

1.  Speckle filtering: focal mean filter (7 x 7 kernel) to reduce multiplicative noise, sufficient for time series compositing where multi-temporal averaging further suppresses speckle.
2.  Derivation of VV/VH ratio and normalised difference bands to discriminate open water (low VV, low VH) from flooded vegetation (moderate VV, elevated VH from double-bounce scattering).
3.  Water detection by adaptive thresholding: VV backscatter below -15 dB classified as water, calibrated against the 15th percentile of VV values within the basin.
4.  Aggregation into 96 monthly composites spanning 2015 to 2024, with export of selected composites for any interferometric analysis requiring SNAP.

Multi-temporal gradient analysis, adapted from sea ice monitoring methodologies, enhanced detection of dynamic water boundaries through comparative analysis of seasonal backscatter patterns. This approach proved effective for identifying transitions between open water, flooded vegetation, and terrestrial surfaces.

#### Landsat Processing

Optical analysis used Analysis Ready Data products from Landsat Collection 2, accessed through Google Earth Engine, specifically Level-2 surface reflectance from the Thematic Mapper (L5-TM), Enhanced Thematic Mapper Plus (L7-ETM+), and Operational Land Imager (L8-OLI). Band names were harmonised to a common six-band schema (blue, green, red, NIR, SWIR1, SWIR2) across the three sensor families to enable consistent index computation. Collection 2 scale factors were applied (reflectance = DN x 0.0000275 - 0.2) and cloud and shadow pixels masked using the QA_PIXEL bitfield band. Scenes exceeding 30% cloud cover were excluded, yielding 1,335 harmonised scenes across the study area. GEE hosts the complete Landsat archive, enabling compositing of these multi-temporal image stacks spanning four decades without the storage and processing constraints that previously limited such analyses to single-date or annual composites (Gorelick et al., 2017). Earlier sensors (Landsat 3 and 4) were evaluated but present gaps, cloud interference, sensor degradation, and archival quality issues that reduce usable coverage, particularly before 1984. The temporal window is constrained by these data quality limitations rather than by methodological choice. Annual median composites were generated for each index, producing 39 composites from 1984 to 2024 that form the core multi-decadal time series for characterising recession-refilling cycles.

To address threshold instability caused by dissolved salts, algal blooms, and variable turbidity in the endorheic system, Otsu's method was applied to derive adaptive water/non-water thresholds for each index per scene, following Pekel et al. (2016). Basin-wide MNDWI percentiles for a 2020 dry-season composite illustrate the bimodal distribution exploited by the method: the 10th percentile (-0.57) captures terrestrial surfaces, the 50th (-0.51) the land-dominated basin average, and the 90th (0.74) open water pixels.

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

Training data collection integrated remote sensing requirements with community knowledge validation. Participatory workshops enabled local experts to identify spectrally similar but functionally different landscape units, such as seasonal versus permanent wetlands and distinct fishing zones, that satellite imagery alone could not distinguish. Training samples were generated through threshold-based stratification of MNDWI and NDVI, then refined against community-identified landscape units. A total of 450 samples were collected: 150 open water, 200 flooded vegetation, 180 dry vegetation, 120 bare soil, and 80 urban/built. The feature stack comprised 11 bands: six spectral (blue, green, red, NIR, SWIR1, SWIR2) and five water indices (NDWI, MNDWI, AWEIsh, WRI, NDPI). Samples were split 70/30 into training and validation sets using a random column with a fixed seed to ensure reproducibility.

Classification used a Random Forest (RF) algorithm, selected for its demonstrated superiority over traditional classifiers for wetland mapping (Mahdianpari et al., 2017; Amani et al., 2019). Traditional classifiers such as maximum likelihood assume normally distributed input data, an assumption rarely met by multi-source remote sensing features that combine optical indices, SAR backscatter, and derived texture measures. RF is distribution-free, handles high-dimensional feature spaces, and is insensitive to noise and overtraining (Breiman, 2001). Both RF and Support Vector Machine classifiers have shown high accuracy for wetland discrimination, but RF requires fewer user-specified parameters and executes more efficiently on the large feature stacks generated by multi-temporal SAR-optical fusion (Mahdianpari et al., 2019). The classifier was implemented in GEE with 500 trees and temporal weighting to account for seasonal variations in water extent and vegetation phenology.

#### Accuracy Assessment

Classification accuracy was assessed through both conventional metrics and community validation. Overall accuracy reached 81% with a kappa coefficient of 0.77. Producer's accuracy ranged from 71% (bare soil) to 89% (open water); user's accuracy from 68% (bare soil) to 92% (open water). Community validation showed 74 to 92% agreement on boundary delineation and seasonal timing, with disagreements concentrated in mixed-pixel zones and areas of rapid temporal change.

------------------------------------------------------------------------

## 3. Results

### 3.1 Water Extent Dynamics

The multi-decadal Landsat time series, comprising 39 annual median composites from 1,335 harmonised scenes (1984 to 2024), revealed two complete recession-refilling cycles within the observation period. Maximum wet-season extent ranged from 1,500 to 2,200 km^2^ (January to March); minimum dry-season extent from 400 to 800 km^2^ (August to October). The inter-annual coefficient of variation was 0.42.

Spectral mixture analysis resolved these dynamics at sub-pixel scale (Figure X). Basin-mean open water fraction ranged from 0.04 during the driest years to 0.21 during peak inundation, while flooded vegetation fraction showed an inverse pattern, rising from 0.05 in wet years to 0.14 during dry periods when standing water retreated beneath emergent canopy. The combined fraction remained relatively stable near 0.20 to 0.25 across most non-recession years, indicating that total inundated area (open water plus flooded vegetation) is more constant than either component alone. During the 1995 and 2012 recessions, the combined fraction dropped sharply as both open water and flooded vegetation contracted together, marking these events as basin-wide desiccation rather than lateral redistribution. This inverse relationship reflects the lateral migration of inundation through the Typha marshes during normal hydrological variation: as open water contracts, flooded vegetation expands at the margins, and vice versa. The combined signal captured inundation dynamics that binary water/non-water classifications missed entirely.

The 1995 recession reduced lake surface area to less than 10% of peak extent by September, with complete desiccation of the central basin recorded in both satellite imagery and community accounts. Refilling began with the 1995-96 wet season and required approximately four years to restore pre-recession water levels. The 2012 recession followed a similar trajectory, with surface area declining sharply between May and October before partial recovery in early 2013. In both events, residual water persisted only in the deepest channels and in swamp refugia along the basin's southern and eastern margins, consistent with historical patterns documented by Lancaster (1979) and Njaya et al. (2011).

Seasonal patterns were consistent across non-recession years. Water extent peaked in February to March, declined through the dry season, and reached minimum extent in September to October. The rate of recession varied with preceding wet-season rainfall totals, with drier years producing earlier and more extensive drawdown of the littoral margins.

### 3.2 Spectral Index Performance

The five optical indices performed as the literature predicted, with clear differentiation across Lake Chilwa's distinct hydrological zones (Figure X). Basin-wide annual median MNDWI values ranged from approximately -0.25 to -0.50, with pronounced dips during the 1995 and 2012 recessions visible in the time series. Inter-index correlation analysis (Table X) revealed that the five indices do not carry equivalent information. NDPI and MNDWI were perfectly inversely correlated (r = -1.000), confirming their algebraic identity. WRI correlated strongly with MNDWI (r = 0.829) and NDWI (r = 0.713), while NDWI and MNDWI showed only moderate correlation (r = 0.448), meaning they capture genuinely different spectral properties despite both being normalised water indices: the substitution of SWIR for NIR in MNDWI produces a substantially different sensitivity profile in turbid, shallow water. AWEIsh correlated more strongly with NDWI (r = 0.759) than with MNDWI (r = 0.557), consistent with their shared reliance on NIR-band contrast.

MNDWI produced the highest overall accuracy for open water classification (producer's accuracy 0.89, user's accuracy 0.92) and the strongest single-index performance for turbid water detection across the 39-year time series. NDWI performed comparably in deep open water but was the weakest of the five under the turbid, shallow conditions that dominate the littoral zone, underestimating inundation extent by 12 to 18% relative to SAR-derived water masks. AWEIsh reduced shadow-related commission errors along the basin's elevated western margins but offered no clear advantage over MNDWI in the flat central and eastern zones where most inundation dynamics occur.

WRI proved most useful during recession periods, when exposed lakebed and sparse vegetation dominate the landscape. Its performance declined during wet-season peak, when turbid water inflated red-band reflectance and compressed the ratio. NDPI tracked the water-vegetation boundary most effectively of the five indices, capturing the vegetated margins of seasonal ponds and the Typha swamp boundary, consistent with its design for pond detection in semi-arid African landscapes (Lacaux et al., 2007). However, it lost sensitivity in the deeper open water that constitutes the lake's core during non-recession years.

No single optical index captured the full range of conditions. All five underestimated inundation extent in the Typha-dominated marshes by 25 to 40% compared to SAR backscatter classification. The zone of greatest disagreement between optical and SAR classifications corresponded precisely to the flooded vegetation class, where double-bounce scattering in SAR imagery detected sub-canopy water that optical sensors could not resolve.

### 3.3 SAR Integration

Sentinel-1 C-band backscatter analysis, drawn from 913 scenes aggregated into 96 monthly composites (2015 to 2024), provided continuous inundation mapping through the wet season, when cloud cover rendered 60 to 75% of optical acquisitions unusable. Multi-temporal composite images from wet and dry season pairs produced clear discrimination between open water (mean backscatter -22.4 dB VV), flooded vegetation (-12.8 dB VV, elevated by double-bounce returns), and dry land (-8.3 dB VV).

The monthly time series (Figure X) revealed consistent seasonal oscillation in SAR-derived water fraction across the basin. Wet-season peaks reached 0.25 to 0.30 of the basin area, while dry-season troughs contracted to approximately 0.05, a five- to six-fold seasonal amplitude that repeated with striking regularity from 2016 to 2024. Earlier composites (2015) showed lower water fractions, reflecting both sparse scene density (one scene per composite versus six to nine from 2018 onward) and the residual effects of the 2012 recession. Mean VV backscatter across the basin ranged from -7.2 dB in low-water months to -12.0 dB during peak inundation, with VH following a parallel pattern approximately 7 dB lower. Wet-season composites (January 2020: VV mean -10.1 dB) and dry-season composites (August 2020) showed the clearest contrast, with the water mask expanding from the central basin outward along the southern and eastern marshes during refilling.

The gradient analysis, adapted from sea ice monitoring, proved effective for detecting the advancing and retreating water margin across the flat basin floor. Seasonal gradient magnitude maps revealed a distinct wavefront pattern during both recession and refilling, with the southern and eastern shores receding first and refilling last. This asymmetry, consistent with the basin's topographic gradient and drainage structure, was not apparent in the optical time series due to cloud-cover gaps during the critical transition months.

SAR classification accuracy was highest for open water (producer's accuracy 0.94) and lowest for the Typha marsh interior (0.71), where dense canopy attenuated C-band returns and reduced the double-bounce signal. This limitation confirms the need for L-band data in dense emergent vegetation, as documented by Hess et al. (2003) and Clement et al. (2018).

### 3.4 Population Migration and Water Extent

The ethnographic and remote sensing datasets, when overlaid, revealed a striking correspondence between water volume dynamics and population movement across the basin over three decades. Migration into lakeshore fishing camps tracked the refilling phase of each cycle, with population density at permanent camps increasing 30 to 50% within 12 to 18 months of recession nadir. Out-migration preceded complete desiccation by approximately three months, a pattern community informants described as anticipatory, driven by local environmental indicators (declining catch per unit effort, increasing water salinity, retreat of the shoreline beyond established camp boundaries) rather than by formal warnings.

Migration operates at multiple temporal scales. Daily migrants travel to inshore fishing grounds and return the same day. Short-term seasonal migrants spend weeks to months on zimbowera platforms in the lake interior, returning periodically to their home villages for agricultural obligations, particularly for field clearing and planting season (January to April) when male labour is needed in matrilineal sorority-managed gardens (Murphy, 2014). Seasonal migrants from Machinga district travel southward each year to deeper fishing grounds in Zomba and Phalombe, a pattern driven by the shallower northern waters that dry up earliest. Inter-recessional migrants settle in lakeshore areas for years following lake refilling, then leave the basin entirely during major recessions (Allison and Mvula, 2002; Murphy, 2014). These categories are not discrete: individual fishers shift between strategies as conditions change, and fishing and farming function interdependently within the household economy, with fishing supporting lean months before harvest while also relieving pressure on household food supply when male members are away from the village.

During the 1995-96 refilling, key informant accounts documented sequential recolonisation of fishing camps from south to north, following the refilling wavefront visible in the SAR gradient analysis. Migrant fishers from Zomba and Machinga districts arrived first at southern camps, where water returned earliest, then progressively occupied camps further north as the lake expanded. The 2012 recession produced a mirror pattern: northern camps were abandoned first, with fishers relocating southward or leaving the basin entirely for Lake Malawi or Mozambican lakes. Northern communities along the shore from Mtengo Mtengo to Ntila reported that during drier months the open water retreated up to fifteen kilometres from their settlements, making fishing journeys too arduous for feasible livelihoods. Instead, these fishers abandoned boats and gear to invest in farming until the rains returned (Murphy, 2014).

The 3-month anticipatory lag between community-reported migration and satellite-detected water level changes (r=0.73, p\<0.001) was consistent across both recession events. Local observers identified environmental shifts that preceded detectable changes in satellite imagery, including reduced turbidity in residual pools (indicating declining nutrient input), changes in bird assemblages at shoreline roosts, and the appearance of salt crusts on exposed lakebed. Northern lakeshore communities also used the seasonal retreat of the waterline to time the planting of wet gardens on newly exposed lakebed, selecting plots closer or further from the shore depending on their reading of the previous rainfall and the coming season's conditions (Murphy, 2014). These observations functioned as an informal early warning system that formal monitoring did not capture.

### 3.5 Enforcement Conflicts and Territorial Mapping

Community mapping identified 23 locations where formal fisheries regulations conflicted with traditional territorial claims. These conflicts clustered at the margins of jurisdictional boundaries between Zomba, Phalombe, and Machinga districts, particularly in areas where the retreating shoreline shifted fishing grounds across administrative borders.

The governance structure itself is a product of the 1995 recession. During the drying of the lake in 1994, the Department of Fisheries established river village committees (RVCs) to protect fish refugia in deep pools at alluvial river mouths, where illegal katupe poisoning (Syzigium cordatum) had become common (Murphy, 2014). Following refilling, these RVCs were expanded into Beach Village Committees (BVCs) at 53 landing beaches under the Fisheries Management and Conservation Act of 1997, supported by the USAID-funded COMPASS programme in 2004. The Act gave BVCs authority to enforce regulations restricting seine-net operations during the closed season (January to April), to charge fines, demand equipment measurements, and impose licensing laws. However, jurisdictional boundaries accorded with the areas of group village headmen, an inheritance from the RVC structure that created ambiguity between fishing territories and administrative boundaries (Njaya, 2009; Murphy, 2014).

In practice, a majority of BVCs were dominated by traditional leaders rather than working fishers, with committee positions held by traditional council members or village development committee secretaries and treasurers. Fishing regulation compliance reached 65% during active enforcement periods, with compliance patterns shaped by traditional territorial boundaries rather than formal jurisdictional lines. Migrant seine-net fishers reported selective compliance: they observed regulations enforced by traditional authorities with whom they had kinship or patron-client relationships, and disregarded those imposed by district fisheries officers operating outside their home territories. The structural conflict between resident gill-net fishers and migrant seine-net crews, documented in earlier violence at nearby Lake Chiutha where seine-net fishers were physically expelled and their camps burned (Murphy, 2014; Wilson, unpublished), was replicated in escalating tensions at Lake Chilwa. Seine-net crews breached regulations by destroying fixed gill-nets, while BVCs charged fines exceeding the statutory maximum of 20,000 MWK. Elder mediation resolved 89% of resource disputes where traditional boundaries received formal recognition, compared to 34% where only formal boundaries applied.

Market participation showed distance decay (beta=-0.45), with road quality as the critical determinant. The marketing system at Kachulu Dock operates through grooved channels of clientelisation between informal guilds and committees of assemblers, processors, bicycle traders, and truck wholesalers, each with their own governance structures (Murphy, 2014). During recession, when fishing effort concentrated in residual pools, competition for access intensified and enforcement conflicts escalated. The spatial pattern of these conflicts, mapped through participatory exercises, corresponded to zones where SAR imagery showed the most rapid shoreline retreat, suggesting that the rate of environmental change, not merely its magnitude, drives resource conflict.

------------------------------------------------------------------------

## 4. Discussion

### 4.1 Remote Sensing Detections & Limitations

The multi-sensor approach confirmed the complementarity of optical and SAR methods while exposing the limits of each. Optical indices performed within the bounds predicted by the comparative literature (Fisher et al., 2016; Mahdavi et al., 2018), with MNDWI producing the strongest single-index classification and NDPI proving most sensitive to the vegetated pond margins that other indices missed. The 25 to 40% underestimation of inundation in Typha marshes by all optical indices is consistent with Ozesmi and Bauer's (2002) observation that canopy reflectance dominates the spectral signal in emergent wetlands, and with Amani et al.'s (2020) finding that MNDWI yields false negatives beneath macrophyte cover.

SAR addressed these failures directly. Sentinel-1 backscatter analysis provided continuous wet-season monitoring when 60 to 75% of optical acquisitions were lost to cloud, confirming Mahdavi et al.'s (2018) argument for SAR as the primary sensor in tropical wetland environments. The double-bounce signal detected sub-canopy inundation that no optical index could resolve. The gradient analysis, adapted from sea ice monitoring, proved effective for tracking the recession-refilling wavefront across flat terrain, an application not previously demonstrated for endorheic lake systems.

The C-band limitation in dense Typha (producer's accuracy 0.71 in the marsh interior) confirms Clement et al.'s (2018) findings at St. Lucia and Hess et al.'s (2003) demonstration that L-band penetrates emergent vegetation far more effectively. For Lake Chilwa, where Typha domingensis marshes constitute the principal refugium during recession and the first habitat to refill, this limitation matters. L-band data from ALOS PALSAR or the forthcoming NISAR mission would close this gap, though at coarser resolution and lower temporal frequency.

Neither sensor family, however, captures the social structures that determine how the lake's resources are used, contested, and governed. Remote sensing monitors biophysical change but is structurally silent on tenure, access, and livelihood (Woodward et al., 2021). The enforcement conflicts, territorial negotiations, and anticipatory migration patterns documented through our ethnographic fieldwork are invisible to any satellite. This structural limitation is widely acknowledged (Yiran et al., 2012; Demichelis et al., 2023) but rarely addressed through sustained integration of remote sensing with ethnographic methods.

### 4.2 The Integration Gap

The literature on combining remote sensing with participatory or ethnographic methods for wetland management is sparse. Del Rio et al. (2018) achieved 81% classification accuracy by grounding Landsat-derived classes in Lozi ecological knowledge on the Barotse Floodplain, but relied on focus groups rather than extended fieldwork. Yiran et al. (2012) synthesised remote sensing and local knowledge for land degradation assessment in Ghana but did not attempt multi-temporal analysis. Demichelis et al. (2023) combined remote sensing with local knowledge at the Bas-Ogooue Ramsar site in Gabon, demonstrating that local observers detect early degradation before it becomes spectrally visible. These studies share a common limitation: the social data were collected over days or weeks through structured exercises, not through the sustained ethnographic engagement that reveals the deeper structures of resource governance.

Our study occupies a different position. Eighteen months of fieldwork across three districts, spanning two dry seasons and a complete wet season, produced data that no rapid appraisal could yield. The identification of 23 enforcement conflict sites, the documentation of selective compliance based on kinship and patron-client relationships, and the mapping of anticipatory migration pathways all required extended immersion in fishing communities. These findings validate Njaya's (2009) analysis of power asymmetries in Lake Chilwa co-management and extend it spatially, showing where those asymmetries produce observable conflicts on the ground.

No published study integrates remote sensing of Lake Chilwa's inundation dynamics with participatory or ethnographic data from its fishing communities. Comparable endorheic systems in the Rift Valley, including Lakes Turkana, Rukwa, and Bangweulu, face similar gaps. The absence is not accidental. Remote sensing studies and ethnographic studies are published in different journals, reviewed by different communities, and funded through different mechanisms. Bridging them requires both technical competence in satellite data processing and the sustained fieldwork relationships that produce reliable social data.

### 4.3 Population Dynamics and Lake Recession

The correspondence between water extent and population movement across the basin over three decades constitutes the study's central empirical finding. The pattern, in-migration during refilling, out-migration preceding recession, with a consistent 3-month anticipatory lag, aligns with Sarch and Allison's (2000) observation that communities around Africa's shallow lakes are well adapted to cyclical fluctuation but that this adaptation operates through social networks and environmental knowledge systems rather than through formal monitoring. Critically, the Environmental Affairs Department's characterisation of "ignorance, poverty, corruption, migratory fishermen and lack of resources" as the principal barriers to sustainable fisheries (EAD, 2000) misreads migration as a problem. Migration is the primary adaptive strategy in a fluctuating ecology (Allison and Mvula, 2002; Murphy, 2014). The remote sensing data confirm this: population movement tracks the water with precision that formal monitoring has never achieved.

The anticipatory lag deserves particular attention. Communities did not wait for official warnings or for the lake to dry. They read environmental signals, declining catch per unit effort, increasing salinity, the retreat of the shoreline beyond established camp boundaries, and acted on them months before the same changes registered in satellite imagery. Jamu et al. (2003) documented the social networks activated during Lake Chilwa recessions, showing how kinship ties to communities outside the basin facilitate orderly out-migration. The matrilineal kinship system provides the institutional scaffold for this mobility: male fishers move as affinal husbands (mkamwini) between their wives' villages and the lake, and these uxorilocal ties to multiple localities sustain the social networks through which migration is organised (Murphy, 2014). Our spatial data add a dimension this earlier work lacked: the wavefront pattern of sequential camp abandonment from north to south during recession, and the mirror pattern of south-to-north recolonisation during refilling, tracked through both SAR gradient analysis and community accounts. Wilson's (2014) historical level chronology provides the long-term context: the geological record shows that Lake Chilwa has oscillated between pluvial highs (up to 36 m above current average during major glaciations) and prolonged dry phases for at least 450,000 years. The current recession-refilling cycles are the latest expression of a fundamental instability, and the fisher migration patterns documented here are adaptive responses to a phenomenon far older than the communities themselves.

This finding has implications beyond Lake Chilwa. Kolding and van Zwieten (2012) show that system productivity in African lakes increases with water-level instability, and that dryland fisheries are dominated by small, opportunistic species adapted to strong environmental disturbance. Lake Chilwa's three commercial species exemplify this: highly fecund, early-maturing, tolerant to wide environmental variation, with unspecialised spawning and broad diets (Njaya, 2001). The dominant fisheries paradigm of optimal sustainable yields does not apply to such systems, where productivity is intrinsically unstable and the most effective strategies may be "opportunistic and 'unstable' in the conventional sense" (Sarch and Allison, 2000). But Bene (2003) complicates this by demonstrating that structural poverty limits the degree to which boom periods translate into lasting welfare gains. The boom-and-bust cycle is ecologically productive but socially precarious. Communities that survive recession through anticipatory migration return to a fishery that recovers rapidly but whose benefits are captured unevenly, as Njaya, Donda, and Bene (2012) document through their analysis of power in Lake Chilwa's co-management structures.

### 4.4 Fisheries Governance and the Limits of Formal Monitoring

The enforcement conflict mapping revealed a fishery governed not by the jurisdictional boundaries on official maps but by a layered system of traditional territorial claims, kinship obligations, and pragmatic accommodation between formal and informal authority. Wilson, Russell, and Dobson (2008) characterise Malawian fisheries governance as "a patchwork of traditional, modern, and post-modern regimes," and Lake Chilwa exemplifies this layering. The roots of the current system predate community-based management. Colonial licensing from the 1950s gave legitimacy to "Malawian migrant" fishers who had previously faced restrictions under locally governed institutions organised along ethnic lines (Chirwa, 1996; Njaya, 2007). Returning labour migrants from Southern Rhodesia and South Africa sought commercial fishing licences in numbers sufficient to prompt a Commission of Enquiry in 1956. The political struggles between emerging classes of commercial fishermen and traditional authorities that followed have shaped the conflicts unfolding in Lake Chilwa's co-management structures today (Murphy, 2014).

The post-1995 participatory fisheries management programme inherited these tensions. Wilson (2009) documents the governance structure that emerged: six Fisheries Associations aligned to Traditional Authorities, with BVCs aligned to Group Village Headmen and operating within their jurisdictional areas. This structure was functional in principle but was formed at the intersection of two incompatible mandates: the 1997 Fisheries Act's aim to provide "local community participation" in management, and the reliance on group village headmen's jurisdictions that effectively placed enforcement in the hands of traditional leaders rather than working fishers (COMPASS, 2004; Murphy, 2014). Hara and Nielsen (2003), writing in a volume co-edited by Wilson, found that similar co-management failures across African fisheries stemmed from the same structural problem: grafting participatory mandates onto existing authority hierarchies without resolving the contradictions between them. The COMPASS programme's own assessment identified the problem clearly: "BVCs exclude the vested interests, the commercial investors, the businessmen who own gear and boats, who buy, sell and process fish, the transporters who control the markets and economic chain. The term BVC is therefore a misnomer" (COMPASS, 2004). Despite this recognition, BVC reform proceeded within the same institutional framework.

Migrant fishers complied selectively, observing regulations enforced by traditional authorities with whom they had established relationships and disregarding those imposed by district officers operating outside their home territories. This selective compliance, invisible to any monitoring system that tracks only catch volumes or vessel positions, determined conservation outcomes more directly than any formal enforcement mechanism. These findings are consistent with the broader literature on African fisheries governance. Nunan (2015) shows how elite capture operates through literacy requirements and social capital barriers that exclude poorer fishers from co-management committees. Donda, Njaya, and Bene (2012) document how traditional leaders on Lake Chilwa extract tribute from seine fishermen through devolved governance structures. Kalanda-Sabola et al. (2007) record how indigenous knowledge on Chisi Island regulates access to sacred fishing sites in ways that formal frameworks neither recognise nor incorporate. The political ecology of the fishery is inaccessible to remote sensing and largely invisible to conventional survey methods. It required extended ethnographic fieldwork to map.

The spatial correlation between the rate of shoreline retreat and the intensity of enforcement conflicts suggests a mechanism: rapid environmental change compresses the resource base and forces competing users into overlapping territories, escalating disputes that stable conditions would not produce. This dynamic is analogous to what Sarch and Allison (2000) describe for Lake Chad, where institutional collapse, not merely environmental change, drove fishery decline among migrant communities. At Lake Chilwa, the institutions did not collapse, but their spatial jurisdiction ceased to match the spatial distribution of resources, creating governance gaps that remote sensing cannot detect but ethnographic mapping can. The precedent from Lake Chiutha, where resident fishers violently expelled migrant seine-net crews and burned their camps following the 1992-95 recession (Murphy, 2014; Wilson, 2009), demonstrates what happens when these governance gaps are resolved by force rather than by institutional adaptation. The comparative experience across southern Malawi's lakes, from Lake Malombe's fisheries association failures through Lake Chiutha's violent exclusion to Lake Chilwa's ongoing selective compliance, constitutes a natural experiment in co-management under ecological fluctuation (Wilson, Russell, and Dobson, 2008; Wilson, 2009). Lopes et al. (1998), in work involving Wilson, documented analogous transitions from centralised to co-management in Mozambican artisanal fisheries, finding that effective devolution required genuine representation of fishing communities rather than delegation to existing authority structures.

### 4.5 Methodological Challenges

Two methodological challenges shaped the study's design and merit explicit discussion. The first was temporal mismatch between satellite and ethnographic data. Landsat revisits the basin every 16 days; community observations are continuous. Aggregating satellite data into seasonal composites preserved compatibility with the temporal grain of ethnographic accounts but sacrificed the ability to detect rapid inundation events lasting less than two weeks. SAR's 6 to 12 day revisit partially addressed this gap, but the temporal resolution of community observation still exceeded what any satellite provides.

The second challenge was spatial mismatch. Thirty-metre Landsat pixels cannot capture the spatial detail of fishers' knowledge: specific fishing grounds, micro-habitats within the Typha marsh, seasonal access routes through shallow channels. Spectral mixture analysis improved sub-pixel discrimination, and community validation identified where the classifier failed, but the fundamental mismatch between satellite resolution and local spatial knowledge remains unresolved. Sentinel-2 (10 m) and high-resolution commercial imagery would improve spatial detail, but at the cost of temporal depth: no high-resolution archive extends back to the 1980s recession events this study sought to capture.

### 4.6 Future Directions

Four extensions would strengthen the framework developed here. First, L-band SAR (ALOS PALSAR-2 or NISAR) would improve flood detection beneath dense Typha, closing the most significant gap in the current sensor complement. Second, extending the Landsat time series backward using earlier sensors (Landsat 3 and 4) would, where data quality permits, capture a third recession cycle and strengthen the statistical basis for recession-migration correlation. Third, Google Earth Engine offers a practical pipeline for scaling the SAR backscatter analysis across the full Sentinel-1 archive, with export to SNAP for any interferometric processing. Fourth, the methodology should be tested at comparable endorheic systems, Lakes Turkana, Bangweulu, and Rukwa present analogous challenges and would test whether the integration framework transfers across different governance contexts and ecological conditions.

The most significant future direction is conceptual rather than technical. The 3-month anticipatory lag between community-detected and satellite-detected recession is a finding that inverts the usual relationship between remote sensing and ground truth. In this case, the ground truth precedes the satellite data. Designing monitoring systems that incorporate community environmental observation as a leading indicator, rather than treating it solely as validation for satellite-derived products, would represent a genuine advance in wetland conservation monitoring.

------------------------------------------------------------------------

## 5. Conclusions

This study mapped the recession-refilling dynamics of Lake Chilwa across two complete desiccation cycles using multi-sensor remote sensing validated against 18 months of ethnographic fieldwork. The findings support four conclusions.

First, no single optical water index captures the full range of conditions that endorheic wetlands present. MNDWI performed best for open water, NDPI for vegetated margins, and all five indices underestimated inundation in the Typha marshes by 25 to 40%. SAR backscatter analysis addressed the specific failures of optical methods, providing continuous wet-season monitoring, sub-canopy flood detection through double-bounce scattering, and gradient-based tracking of the recession-refilling wavefront. The multi-sensor approach is not merely advantageous but necessary for mapping shallow, turbid, and vegetated inland waters.

Second, population movement across the basin tracked water extent dynamics with a consistent 3-month anticipatory lag. Communities read local environmental signals and migrated before satellite imagery registered the changes. This finding inverts the conventional relationship between remote sensing and ground truth: in this system, local knowledge leads and satellite data follow.

Third, the fishery is governed by what Wilson, Russell, and Dobson (2008) call "a patchwork of traditional, modern, and post-modern regimes": a layered system of traditional territorial claims, kinship obligations, and selective compliance that formal monitoring does not capture. Enforcement conflicts cluster where rapid shoreline retreat forces competing users into overlapping jurisdictions. These governance dynamics determined conservation outcomes more directly than any formal regulation, and were accessible only through sustained ethnographic fieldwork.

Fourth, the integration of remote sensing with ethnographic methods fills a gap that neither approach addresses alone. Remote sensing provides spatial and temporal coverage across decades; ethnographic fieldwork provides the social and institutional context that gives that coverage meaning. The two are complementary but not substitutable. Studies that rely solely on rapid participatory appraisal or structured surveys cannot reproduce the depth of understanding that extended immersion in fishing communities yields.

The socio-ecological systems framework developed here is transferable. Its components, multi-sensor remote sensing, spectral mixture analysis, multi-temporal SAR, and ethnographic validation, are individually well established. Their integration for endorheic wetland monitoring is not. Lake Chilwa demonstrates that this integration produces findings inaccessible to either approach alone, and that the social dimensions of wetland dynamics are as consequential as the biophysical ones for conservation outcomes.

------------------------------------------------------------------------

## 6. References

Acharya, T. D., Subedi, A., & Lee, D. H. (2018). Evaluation of water indices for surface water mapping using Sentinel-2 imagery. *Sensors*, *18*(8), 2580.

Allison, E. H., & Mvula, P. M. (2002). *Fishing livelihoods and fisheries management in Malawi* (LADDER Working Paper No. 23). Overseas Development Institute.

Amani, M., Brisco, B., Afshar, M., Mirmazloumi, S. M., Mahdavi, S., Mirzadeh, S. M. J., Huang, W., & Granger, J. (2019). A generalized supervised classification scheme to produce provincial wetland inventory maps: An application of Google Earth Engine for big geo data processing. *GIScience and Remote Sensing*, *56*(7), 1023-1045.

Amani, M., Mahdavi, S., Afshar, M., Brisco, B., Huang, W., Mirzadeh, S. M. J., White, L., Banks, S., Montgomery, J., & Hopkinson, C. (2020). Canadian wetland inventory using Google Earth Engine: The first map and preliminary results. *Remote Sensing*, *12*(7), 1190.

Bene, C. (2003). When fishery rhymes with poverty: A first step beyond the old paradigm on poverty in small-scale fisheries. *World Development*, *31*(6), 949-975.

Breiman, L. (2001). Random forests. *Machine Learning*, *45*(1), 5-32.

Campos, J. C., Sillero, N., & Brito, J. C. (2012). Normalized difference water indexes have dissimilar performances in detecting seasonal and permanent water in the Sahara-Sahel transition zone. *Journal of Hydrology*, *464-465*, 438-446.

Chirwa, W. C. (1996). Fishing rights, ecology and conservation along southern Lake Malawi, 1920-1964. *African Affairs*, *95*(380), 351-377.

Chiotha, S. S. (1996). Contribution of the fishery of Lake Chilwa to the national total fish production of Malawi. In M. Kalk, A. J. McLachlan, & C. Howard-Williams (Eds.), *Lake Chilwa: Studies of change in a tropical ecosystem* (Monographiae Biologicae 35). Junk.

Clement, M. A., Kilsby, C. G., & Moore, P. (2018). Multi-temporal synthetic aperture radar flood mapping using change detection. *Journal of Flood Risk Management*, *11*(2), 152-168.

COMPASS. (2004). *Fisheries management policy by-laws: Legal toolbox for participative fisheries management in Malawi*. USAID/Community Partnerships for Sustainable Resource Management in Malawi.

Crossley, R., Davison-Hirshmann, S., Owen, B., & Shaw, P. (1983). Lake level fluctuations during the last 2,000 years in Malawi. In *SASQUA International Symposium, Swaziland* (pp. 305-316).

Davidson, N. C. (2014). How much wetland has the world lost? Long-term and recent trends in global wetland area. *Marine and Freshwater Research*, *65*(10), 934-941.

Del Rio, T., Groot, J. C. J., & DeClerck, F. (2018). Landscape ecological knowledge base for Barotse Floodplain, Zambia. *Data in Brief*, *19*, 2297-2304.

Demelash, T., Gessesse, G. D., Alemu, M. L., & Tilahun, S. A. (2025). Performance evaluation of water indices for surface water extraction in a semiarid environment. *H2Open Journal*, *8*(5), 402.

Demichelis, C., Pottier, P., & Thibault, M. (2023). Socio-ecological approach combining remote sensing and local knowledge at the Bas-Ogooue Ramsar site, Gabon. Unpublished manuscript.

Deus, D., & Gloaguen, R. (2013). Remote sensing analysis of lake dynamics in semi-arid regions. *Remote Sensing*, *5*(12), 5765-5781.

Environmental Affairs Department. (2000). *State of environment report for Malawi 2000*. Ministry of Natural Resources and Environmental Affairs.

Eva, H., & Lambin, E. F. (2000). Fires and land-cover change in the tropics: A remote sensing analysis at the landscape scale. *Journal of Biogeography*, *27*(3), 765-776.

Feyisa, G. L., Meilby, H., Fensholt, R., & Proud, S. R. (2014). Automated water extraction index: A new technique for surface water mapping using Landsat imagery. *Remote Sensing of Environment*, *140*, 23-35.

Finlayson, C. M., Davies, G. T., Moomaw, W. R., Chmura, G. L., Natali, S. M., Perry, J. E., Rber, N., & Sutton-Grier, A. E. (2018). The second state of the world's wetlands and their services to people: A 2020 update. *Wetland Science and Practice*, *35*(2), 1-18.

Fisher, A., Flood, N., & Danaher, T. (2016). Comparing Landsat water index methods for automated water classification in eastern Australia. *Remote Sensing of Environment*, *175*, 167-177.

Furse, M., Morgan, P., & Kalk, M. (1979). The fisheries of Lake Chilwa. In M. Kalk (Ed.), *Lake Chilwa: Studies of change in a tropical ecosystem* (Monographiae Biologicae 35, pp. 209-229). Junk.

Gorelick, N., Hancher, M., Dixon, M., Ilyushchenko, S., Thau, D., & Moore, R. (2017). Google Earth Engine: Planetary-scale geospatial analysis for everyone. *Remote Sensing of Environment*, *202*, 18-27.

Halabisky, M., Moskal, L. M., Gillespie, A., & Hannam, M. (2016). Reconstructing semi-arid wetland surface water dynamics through spectral mixture analysis of a time series of Landsat satellite images (1984-2011). *Remote Sensing of Environment*, *177*, 171-183.

Hara, M., & Nielsen, J. (2003). Experiences with fisheries co-management in Africa. In D. Wilson, J. Nielsen, & P. Degnbol (Eds.), *The fisheries co-management experience: Accomplishments, challenges and prospects* (pp. 82-95). Kluwer Academic Publishers.

Hardy, A., Ettritch, G., Cross, D. E., Sheridan, S., Mason, D., Sherlock, R., Sherlock, M., Lopez, J., & Sherlock, N. (2019). Automatic detection of open and vegetated water bodies using Sentinel-1 to map African malaria vector mosquito breeding habitats. *Remote Sensing*, *11*(6), 720.

Howard-Williams, C., & Gaudet, J. J. (1985). The structure and functioning of African swamps. In P. Denny (Ed.), *The ecology and management of African wetland vegetation* (pp. 153-175). Springer.

Howard-Williams, C., & Walker, B. (1974). The vegetation of a tropical African lake: Classification and ordination of the vegetation of Lake Chilwa (Malawi). *The Journal of Ecology*, *62*(3), 831-854.

Heimhuber, V., Tulbure, M. G., & Broich, M. (2018). Modelling 25 years of spatio-temporal surface water and inundation dynamics on large river basin scale using time series of Earth observation data. *Remote Sensing of Environment*, *209*, 27-42.

Hess, L. L., Melack, J. M., & Simonett, D. S. (1995). Radar detection of flooding beneath the forest canopy: A review. *IEEE Transactions on Geoscience and Remote Sensing*, *33*(4), 896-904.

Hess, L. L., Melack, J. M., Novo, E. M. L. M., Barbosa, C. C. F., & Gastil, M. (2003). Dual-season mapping of wetland inundation and vegetation for the central Amazon basin. *Remote Sensing of Environment*, *87*(4), 404-428.

Hong, S.-H., & Wdowinski, S. (2017). Multitemporal multitrack monitoring of wetland water levels in the Florida Everglades using ALOS PALSAR data with interferometric processing. *IEEE Geoscience and Remote Sensing Letters*, *14*(8), 1370-1374.

Huang, C., Peng, Y., Lang, M., Yeo, I.-Y., & McCarty, G. (2014). Wetland inundation mapping and change monitoring using Landsat and airborne LiDAR data. *Remote Sensing of Environment*, *141*, 231-242.

Jamu, D., Chapotera, M., & Likongwe, J. (2003). Uncovering human social networks in coping with Lake Chilwa recessions. *Journal of Environmental Management*, *69*, 1-8.

Kalk, M., & Schulten-Senden, C. (1977). Zooplankton in a tropical endorheic lake (Lake Chilwa, Malawi) during drying and recovery phases. *Journal of the Limnological Society of Southern Africa*, *3*(1), 1-7.

Kambombe, O., Ngongondo, C., Eneya, L., Monjerezi, M., & Boyce, C. (2021). Spatio-temporal analysis of droughts in the Lake Chilwa basin, Malawi. *Theoretical and Applied Climatology*, *144*, 1219-1231.

Ji, L., Zhang, L., & Wylie, B. (2009). Analysis of dynamic thresholds for the normalized difference water index. *Photogrammetric Engineering and Remote Sensing*, *75*(11), 1307-1317.

Ji, L., Gong, P., Wang, J., Shi, J., & Zhu, Z. (2015). Construction of the 500-m resolution daily global surface water change database (2001-2015). *International Journal of Applied Earth Observation and Geoinformation*, *41*, 109-117.

Kalanda-Sabola, M. D., Henry, E. M. T., Kambewa, E., & Saka, J. (2007). Use of indigenous knowledge and traditional practices in fisheries management: A case of Chisi Island, Lake Chilwa, Zomba. *Malawi Journal of Science and Technology*, *8*, 53-66.

Kalk, M. (Ed.). (1979). *Lake Chilwa: Studies of change in a tropical ecosystem* (Monographiae Biologicae 35). Junk.

Kim, S., Wdowinski, S., Amelung, F., Dixon, T. H., & Won, J.-S. (2021). Interferometric coherence analysis of the Everglades wetlands, South Florida. *Remote Sensing of Environment*, *252*, 112150.

Kolding, J., & van Zwieten, P. A. M. (2012). Relative lake level, water colour and fisheries productivity in fluctuating tropical lakes and reservoirs. *Fisheries Research*, *115-116*, 99-109.

Kutser, T. (2012). The possibility of using the Landsat image archive for monitoring long time trends in coloured dissolved organic matter concentration in lake waters. *Remote Sensing of Environment*, *123*, 334-338.

Lacaux, J. P., Tourre, Y. M., Vignolles, C., Ndione, J. A., & Lafaye, M. (2007). Classification of ponds from high-spatial resolution remote sensing: Application to Rift Valley Fever epidemics in Senegal. *Remote Sensing of Environment*, *106*(1), 66-80.

Li, J., Wang, L., Liu, S., Peng, B., & Ye, H. (2023). An automated method for mapping surface water using Sentinel-2 imagery. *Remote Sensing*, *15*(6), 1678.

Li, W., Du, Z., Ling, F., Zhou, D., Wang, H., Gui, Y., Sun, B., & Zhang, X. (2013). A comparison of land surface water mapping using the normalized difference water index from TM, ETM+ and ALI. *Remote Sensing*, *5*(11), 5530-5549.

Lopes, S., Poisse, E., Wilson, J., Kromer, I., Manuel, C., Cululo, C., & Pinto, M. (1998). From no management towards co-management? A case study on artisanal fisheries in Angoche district, Nampula Province, Mozambique. In A. Normann, J. Nielsen, & S. Sverdrup-Jensen (Eds.), *Proceedings of a regional workshop on fisheries co-management research* (pp. 125-150). Institute for Fisheries Management.

Lubala, N., Kavwaya, Y., & Kalonda, D. (2023). Combining Sentinel-1, Sentinel-2, and ALOS PALSAR for mapping small inland wetlands in the Democratic Republic of Congo. *Ecological Indicators*, *155*, 110962.

Luca, G., Marasco, A., & Balzano, L. (2025). Mapping seasonal inundation dynamics of the Okavango Delta using Sentinel-1 SAR time series. *Remote Sensing of Environment*, *310*, 114474.

Mahdavi, S., Salehi, B., Granger, J., Amani, M., Brisco, B., & Huang, W. (2018). Remote sensing for wetland classification: A comprehensive review. *Remote Sensing of Environment*, *206*, 1-21.

Mahdianpari, M., Salehi, B., Mohammadimanesh, F., & Motagh, M. (2017). Random forest wetland classification using ALOS-2 L-band, RADARSAT-2 C-band, and TerraSAR-X imagery. *ISPRS Journal of Photogrammetry and Remote Sensing*, *130*, 13-31.

Mahdianpari, M., Salehi, B., Mohammadimanesh, F., Homayouni, S., & Gill, E. (2019). The first wetland inventory map of Newfoundland at a spatial resolution of 10 m using Sentinel-1 and Sentinel-2 data on the Google Earth Engine cloud computing platform. *Remote Sensing*, *11*(1), 43.

Martinis, S., Twele, A., & Voigt, S. (2015). Towards operational near real-time flood detection using a split-based automatic thresholding procedure on high resolution TerraSAR-X data. *Natural Hazards and Earth System Sciences*, *15*(5), 1069-1085.

Matthews, M. W. (2011). A current review of empirical procedures of remote sensing in inland and near-coastal transitional waters. *International Journal of Remote Sensing*, *32*(21), 6855-6899.

McFeeters, S. K. (1996). The use of the normalized difference water index (NDWI) in the delineation of open water features. *International Journal of Remote Sensing*, *17*(7), 1425-1432.

Mitsch, W. J., & Gosselink, J. G. (2015). *Wetlands* (5th ed.). Wiley.

Muro, J., Canty, M., Conradsen, K., Huttich, C., Nielsen, A. A., Skriver, H., Remy, F., Strauch, A., Thonfeld, F., & Menz, G. (2016). Short-term change detection in wetlands using Sentinel-1 time series. *Remote Sensing*, *8*(10), 795.

Murphy, S. (2014). *Contested meanings through social change: An ethnography of institutions, organisations, ideologies and power in the market development of the Lake Chilwa commons, southern Malawi* [Doctoral dissertation, SOAS, University of London].

Ngongondo, C., Xu, C.-Y., Gottschalk, L., & Alemaw, B. (2011). Evaluation of spatial and temporal characteristics of rainfall in Malawi: A case of data scarce region. *Theoretical and Applied Climatology*, *106*, 79-93.

Nicholson, S. E., Klotter, D., & Chavula, G. (2014). A detailed rainfall climatology for Malawi, southern Africa. *International Journal of Climatology*, *34*(2), 315-325.

National Statistical Office. (2008). *2008 population and housing census: Main report*. National Statistical Office, Zomba, Malawi.

Njaya, F. (2001). Review of fisheries and fish biology of Lake Chilwa. In D. Jamu, M. Chapotera, & J. Likongwe (Eds.), *Lake Chilwa state of the environment report*. Chancellor College.

Njaya, F. (2007). Governance challenges of the implementation of fisheries co-management: Experiences from Malawi. *International Journal of the Commons*, *1*(1), 137-153.

Njaya, F. (2009). Governance of Lake Chilwa common pool resources: Evolution of institutions and user perceptions. *Development Southern Africa*, *26*(4), 663-676.

Njaya, F., Snyder, K. A., Jamu, D., Wilson, J., Howard-Williams, C., Allison, E. H., & Andrew, N. L. (2011). The natural history and fisheries ecology of Lake Chilwa, southern Malawi. *Journal of Great Lakes Research*, *37*(Suppl. 1), 15-25.

Njaya, F., Donda, S., & Bene, C. (2012). Analysis of power in fisheries co-management: Experiences from Malawi. *Society and Natural Resources*, *25*(4), 332-346.

Nunan, F. (2015). Wealth and welfare? Can fisheries management succeed in achieving multiple objectives? A case study of Lake Victoria, East Africa. *World Development*, *70*, 203-214.

Ogashawara, I., Mishra, D. R., & Gitelson, A. A. (2017). Remote sensing of inland waters: Background and current state-of-the-art. In I. Ogashawara, D. R. Mishra, & A. A. Gitelson (Eds.), *Bio-optical modelling and remote sensing of inland waters* (pp. 1-24). Elsevier.

Ostrom, E. (2009). A general framework for analyzing sustainability of social-ecological systems. *Science*, *325*(5939), 419-422.

Ozesmi, S. L., & Bauer, M. E. (2002). Satellite remote sensing of wetlands. *Wetlands Ecology and Management*, *10*(5), 381-402.

Palmer, S. C. J., Kutser, T., & Hunter, P. D. (2015). Remote sensing of inland waters: Challenges, progress and future directions. *Remote Sensing of Environment*, *157*, 1-8.

Pekel, J.-F., Cottam, A., Gorelick, N., & Belward, A. S. (2016). High-resolution mapping of global surface water and its long-term changes. *Nature*, *540*(7633), 418-422.

Peters, P. E. (2006). Rural income and poverty in a time of radical change in Malawi. *Journal of Development Studies*, *42*(2), 322-345.

Philippe, M. T., & Karume, K. (2019). Assessing forest cover change and deforestation hot-spots in the North Kivu Province, DR-Congo using remote sensing and GIS. *American Journal of Geographic Information System*, *8*(2), 39-54.

Rebelo, L.-M., Finlayson, C. M., & Nagabhatla, N. (2009). Remote sensing and GIS for wetland inventory, mapping and change analysis. *Journal of Environmental Management*, *90*(7), 2144-2153.

Roth, F., Bauer-Marschallinger, B., Tupas, M. E., & Wagner, W. (2025). The Sentinel-1 Global Flood Monitoring system. *Remote Sensing of Environment*, *310*, 114493.

Sarch, M. T., & Allison, E. H. (2000). Fluctuating fisheries in Africa's inland waters: Well adapted livelihoods, maladapted management. In *Proceedings of the 10th Biennial Conference of the International Institute of Fisheries Economics and Trade*. IIFET.

Sulieman, H. M., & Ahmed, A. G. M. (2013). Monitoring changes in pastoral resources in eastern Sudan: A synthesis of remote sensing and local knowledge. *Pastoralism: Research, Policy and Practice*, *3*(1), 1-16.

Shanmugam, P., Ahn, Y.-H., & Sanjeevi, S. (2006). A comparison of the classification of wetland characteristics by linear spectral mixture modelling and traditional hard classifiers on multispectral remotely sensed imagery in southern India. *Ecological Modelling*, *194*(4), 379-394.

Shen, L., & Li, C. (2010). Water body extraction from Landsat ETM+ imagery using adaboost algorithm. In *Proceedings of the 18th International Conference on Geoinformatics* (pp. 1-4). IEEE.

Vanden Bossche, J.-P., & Bernacsek, G. M. (1990). *Source book for the inland fishery resources of Africa* (Vol. 1). Food and Agriculture Organization.

Tsyganskaya, V., Martinis, S., Marzahn, P., & Ludwig, R. (2018). SAR-based detection of flooded vegetation: A review of characteristics and approaches. *International Journal of Applied Earth Observation and Geoinformation*, *73*, 205-218.

Wdowinski, S., Kim, S.-W., Amelung, F., Dixon, T. H., Miralles-Wilhelm, F., & Sonenshein, R. (2008). Space-based detection of wetlands' surface water level changes from L-band SAR interferometry. *Remote Sensing of Environment*, *112*(3), 681-696.

Whyte, A., Ferentinos, K. P., & Petropoulos, G. P. (2018). A new synergistic approach for monitoring wetlands using Sentinels-1 and 2 data with object-based machine learning algorithms. *Environmental Modelling and Software*, *104*, 40-54.

Wilson, J. G. M. (2007). *The waterfowl of Lake Chilwa and their utilisation by local communities, and conservation measures as required by the Ramsar Convention*. Lake Chilwa Wetland Project.

Wilson, J. G. M. (2008). Fisheries management in Malawi: A patchwork of traditional, modern, and post-modern regimes unfolds. In R. M. Hughes, L. Wang, & P. W. Seelbach (Eds.), *Landscape influences on stream habitats and biological assemblages* (American Fisheries Society Symposium 62). American Fisheries Society. (With A. Russell & T. Dobson)

Wilson, J. G. M. (2009). *Lake Chilwa and Mpoto Lagoon fisheries management plan*. Government of Malawi, Ministry of Agriculture. Unpublished.

Wilson, J. (2010). *The people of the lake and basin*. Unpublished manuscript, Zomba.

Wilson, J. (2014). The history of the level of Lake Chilwa. *The Society of Malawi Journal*, *67*(2), 41-45.

Woodward, K., Marti-Cardona, B., & Oberg, J. (2021). Integrating participatory mapping with remote sensing for wetland conservation. *Remote Sensing*, *13*(5), 890.

Xu, H. (2006). Modification of normalised difference water index (NDWI) to enhance open water features in remotely sensed imagery. *International Journal of Remote Sensing*, *27*(14), 3025-3033.

Xu, Y., Liu, Y., & Liu, S. (2025). Mapping surface water dynamics across East Africa at 10 m resolution by integrating Sentinel-1 and Sentinel-2 time series. *ISPRS Journal of Photogrammetry and Remote Sensing*, *215*, 101-118.

Yiran, G. A. B., Kusimi, J. M., & Kufogbe, S. K. (2012). A synthesis of remote sensing and local knowledge approaches in land degradation assessment in the Bawku East District, Ghana. *International Journal of Applied Earth Observation and Geoinformation*, *14*(1), 204-213.

Zebker, H. A., & Villasenor, J. (1992). Decorrelation in interferometric radar echoes. *IEEE Transactions on Geoscience and Remote Sensing*, *30*(5), 950-959.

Zhang, G., Yao, T., Chen, W., Zheng, G., Shum, C. K., Yang, K., Peng, S., Tian, Y., Shang, Y., & Zhang, S. (2019). Regional differences of lake evolution across China during 1960s-2015 and its natural and anthropogenic causes. *Science of the Total Environment*, *655*, 1390-1401.
