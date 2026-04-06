# Manuscript Update Changelog

Working file: `_staging/Manuscript_2025-08_wilson-citations.md`
Source (your working copy): `_staging/Manuscript_2025-08_integrated.md`

All future edits by Claude will be made ONLY in the `_wilson-citations.md` copy.
Line numbers reference the wilson-citations copy.

---

## Session 1 updates (already applied to integrated.md)

These changes were made directly to `Manuscript_2025-08_integrated.md` earlier in this session. You may already be integrating these into your own working copy.

### 1. Title updated (line 1)
- **Old:** "Wetland Inundation Mapping through Remote Sensing and Ethnographic Validation: A Socio-Ecological Study of Lake Chilwa Basin"
- **New:** "Mapping inundation dynamics and fisher migration across Lake Chilwa's recession cycles: Integrating SAR-optical remote sensing with ethnographic fieldwork in an endorheic wetland"

### 2. Introduction opening restructured (Section 1, paragraphs 1-4)
- **What changed:** Replaced the original 3-paragraph opening with 5 paragraphs following Mahdianpari et al. (2019) structure.
- **Added:** Global wetland significance and inventory gaps (para 1), RS methodological challenges (para 2), Sentinel/Landsat/GEE enabling context (para 3), socio-ecological framing narrowing to Lake Chilwa (para 4), study scope paragraph updated to mention GEE pipeline (para 5).
- **New citations:** Mitsch and Gosselink (2015), Davidson (2014), Finlayson et al. (2018), Rebelo et al. (2009), Muro et al. (2016), Gorelick et al. (2017), Mahdianpari et al. (2019).

### 3. SAR Processing rewritten around GEE pipeline (Section 2.3, SAR Processing)
- **What changed:** Replaced generic SNAP-style workflow with GEE pipeline description. Sentinel-1 GRD archive accessed via GEE with pre-applied calibration and terrain correction. Processing steps reduced from 5 to 4 (radiometric calibration now handled by GEE).
- **New citations:** Gorelick et al. (2017), Mahdianpari et al. (2019).

### 4. Landsat Processing updated with GEE context (Section 2.3, Landsat Processing)
- **What changed:** Added "accessed through Google Earth Engine" and sentence on GEE hosting the complete archive.
- **New citation:** Gorelick et al. (2017).

### 5. RF classification context expanded (Section 2.3, Training Samples and Classification)
- **What changed:** Added paragraph justifying Random Forest over traditional classifiers (maximum likelihood) and SVM, citing distribution-free properties, noise insensitivity, and efficiency on multi-temporal feature stacks. Classifier now noted as "implemented in GEE."
- **New citations:** Breiman (2001), Mahdianpari et al. (2017), Mahdianpari et al. (2019).

### 6. New Section 2.2 Reference Data added
- **What changed:** Inserted new subsection between SES Framework (2.1) and Remote Sensing Framework (now 2.3). Describes GPS recording, training/testing polygon generation, 50/50 split protocol following Mahdianpari et al. (2019), and community validation workflow.
- **New citation:** Mahdianpari et al. (2019).

### 7. Study Area (1.3) expanded from 2 to 5 paragraphs
- **What changed:** Added basin surface area (2,310 km2), population density (321/km2), fishery share of national catch (20-43%), three commercial species with evolutionary adaptations, matrilineal Yao-Nyanja-Lomwe context, zimbowera floating cooperatives (Andere, Chambwalu, Lingoni), seine-net camp structure, gill-net family cooperatives, Kachulu port hub, bicycle trader network.
- **New citations:** Kalk (1979), NSO (2008), Njaya et al. (2014), Peters (2006), Chiotha (1996), Njaya (2001), Allison and Mvula (2002), Murphy (2014).

### 8. Ethnographic Methods (2.1) enriched with thesis field sites
- **What changed:** Named BVC members and occupational groups in focus group descriptions. Added specific field sites: Kachulu, Napali, Andere, Chambwalu, Lingoni, Manda Manjeza. Added Mapila Canal (10-15 km), Chisi Island sacred sites, district-level migration routes.
- **New citations:** Murphy (2014), Kalanda-Sabola et al. (2007).

### 9. Results 3.4 (Population Migration) enriched with thesis migration typology
- **What changed:** Added four-category migration typology (daily, short-term seasonal, seasonal, inter-recessional). Added matrilineal agricultural calendar linkage. Added northern shore retreat distances (15 km from Mtengo Mtengo to Ntila). Added wet garden planting on exposed lakebed.
- **New citation:** Murphy (2014).

### 10. Results 3.5 (Enforcement Conflicts) expanded with BVC governance history
- **What changed:** Added full BVC governance timeline: 1994 RVCs, 1997 Fisheries Act, 2004 COMPASS programme. Added COMPASS self-critique quote on BVC representation failures. Added statutory fine ceiling violations. Added seine-net vs gill-net structural conflict. Added Lake Chiutha violent expulsion precedent.
- **New citations:** Murphy (2014), Njaya (2009), COMPASS (2004), Wilson (unpublished).

### 11. Discussion 4.3 (Population Dynamics) strengthened
- **What changed:** Added EAD's mischaracterisation of migration as a barrier. Added matrilineal kinship scaffold for migration (mkamwini/uxorilocal ties). Added species-level argument against optimal sustainable yield paradigm. Added Sarch and Allison (2000) quote on "opportunistic and unstable" strategies.
- **New citations:** EAD (2000), Murphy (2014), Njaya (2001).

### 12. Discussion 4.4 (Fisheries Governance) strengthened
- **What changed:** Added colonial licensing history from 1950s, 1956 Commission of Enquiry, class formation among commercial fishermen, COMPASS self-assessment quote, Lake Chiutha violent expulsion precedent.
- **New citations:** Chirwa (1996), Njaya (2007), COMPASS (2004), Murphy (2014), Wilson (unpublished).

### 13. References section: 21 new citations added
Allison and Mvula (2002), Breiman (2001), Chirwa (1996), Chiotha (1996), COMPASS (2004), Davidson (2014), EAD (2000), Finlayson et al. (2018), Gorelick et al. (2017), Kalk (1979), Mahdianpari et al. (2017), Mahdianpari et al. (2019), Mitsch and Gosselink (2015), Muro et al. (2016), Murphy (2014), Njaya (2001), Njaya (2007), NSO (2008), Peters (2006), Rebelo et al. (2009), Wilson (unpublished).

---

## Session 2 updates (Wilson citations, in wilson-citations.md ONLY)

All edits below are in `_staging/Manuscript_2025-08_wilson-citations.md` only. Your working copy (`Manuscript_2025-08_integrated.md`) is untouched.

### W1. Study Area, paragraph 1 (Section 1.3, ~line 59)
- **What changed:** Replaced "321 persons km^-2^ (NSO, 2008; Njaya et al., 2014; Peters, 2006)" with "162 persons km^-2^ around the lake itself" per Wilson (2010) census data. Added Ramsar designation sentence citing Wilson (2007). Added poverty statistics (70-81% below poverty line, 0.35 ha average plot) citing Wilson (2010).
- **Citations added:** Wilson (2007), Wilson (2010)

### W2. Study Area, paragraph 2 (Section 1.3, ~line 61)
- **What changed:** Replaced "Since 1900, major recessions... (Murphy, 2014)" with Wilson (2014) as primary source for recession chronology. Added geological context: average max depth 2.95 m, 1760-1850 prolonged dry phase depositing sandy clay substrate. Updated recession date list to match Wilson (2014).
- **Citations added:** Wilson (2014), Crossley et al. (1983)

### W3. Study Area, paragraph 3 (Section 1.3, ~line 63)
- **What changed:** Added sentence on 1979 bumper harvest mechanism (24,310 tonnes, oxidation of organic matter during drying released nutrients upon refilling).
- **Citation added:** Wilson (2014)

### W4. Study Area, paragraph 4 (Section 1.3, ~line 65)
- **What changed:** Added settlement history sentence: Akafula hunter-gatherers, 16th-century Maravi/Mang'anja, mid-19th-century Yao migration from Mozambique, 20th-century Lomwe influx.
- **Citation added:** Wilson (2010)

### W5. Methods, SES Framework (Section 2.1, ~line 87)
- **What changed:** Added sentence after study site boundaries description, noting the fisheries governance structure (6 Fisheries Associations aligned to Traditional Authorities, 53 BVCs aligned to Group Village Headmen).
- **Citation added:** Wilson (2009)

### W6. Discussion 4.3, paragraph 2 (~line 234)
- **What changed:** Added closing sentences to the anticipatory lag paragraph, contextualising current recession cycles within Wilson's 450,000-year geological record of pluvial/inter-pluvial oscillations. Notes that fisher migration is an adaptive response to phenomena far older than the communities.
- **Citation added:** Wilson (2014)

### W7. Discussion 4.4, paragraph 1 (~line 240)
- **What changed:** Added Wilson, Russell, and Dobson (2008) "patchwork of traditional, modern, and post-modern regimes" characterisation as opening frame for the governance discussion.
- **Citation added:** Wilson, Russell, and Dobson (2008)

### W8. Discussion 4.4, paragraph 2 (~line 242)
- **What changed:** Added Wilson (2009) Fisheries Management Plan as documenting the 6 FA / BVC governance structure. Added Hara and Nielsen (2003), from volume co-edited by Wilson, on structural co-management failures across Africa.
- **Citations added:** Wilson (2009), Hara and Nielsen (2003)

### W9. Discussion 4.4, final paragraph (~line 246)
- **What changed:** Replaced "Wilson, unpublished" with "Wilson, 2009" for Lake Chiutha precedent. Added comparative framing across Malombe/Chiutha/Chilwa as "natural experiment in co-management under ecological fluctuation." Added Lopes et al. (1998) on Mozambican co-management transitions (Wilson co-author).
- **Citations added:** Wilson (2009), Wilson, Russell, and Dobson (2008), Lopes et al. (1998)

### W10. Conclusions, paragraph 3 (~line 270)
- **What changed:** Added Wilson, Russell, and Dobson (2008) "patchwork" characterisation to the governance conclusion.
- **Citation added:** Wilson, Russell, and Dobson (2008)

### W11. References section
- **New Wilson-authored citations added (7):**
  - Wilson, J.G.M. (2007) Waterfowl of Lake Chilwa / Ramsar Convention. Lake Chilwa Wetland Project.
  - Wilson, J.G.M. (2009) Lake Chilwa and Mpoto Lagoon Fisheries Management Plan. Unpublished.
  - Wilson, J. (2010) The people of the lake and basin. Unpublished manuscript. Zomba.
  - Wilson, J. (2014) The history of the level of Lake Chilwa. *Society of Malawi Journal*, 67(2), 41-45.
  - Wilson, J.G., Russell, A. and Dobson, T. (2008) Fisheries management in Malawi. *AFS Symposium*, 62.
  - Lopes, S., Poisse, E., Wilson, J. et al. (1998) From no management towards co-management. IFM, 125-150.
  - Hara, M. and Nielsen, J. (2003) Experiences with fisheries co-management in Africa. In: Wilson et al. (eds.) Kluwer, 82-95.
- **Other new citation:** Crossley et al. (1983) SASQUA symposium.
- **Removed:** Wilson (unpublished) replaced by Wilson (2009) where applicable.

---

## Session 3 updates (GEE workflow integration, in wilson-citations.md ONLY)

All edits below integrate preliminary results and technical parameters from 9 DOCX files containing the GEE analytical workflow (conclusions.docx, discussion.docx, environment.docx, index.docx, introduction.docx, methods-landsat.docx, methods-sar.docx, methods-ses.docx, results.docx). Light-touch integration only; consolidation deferred to a later session.

### G1. Remote Sensing Framework, opening paragraph (Section 2.3)
- **What changed:** Added AOI derivation from HydroSHEDS level-6 basin boundaries. Added hydrological base layers: MERIT Hydro (v1.0.3), HydroSHEDS (v1.1), JRC Global Surface Water (v1.4; 454 monthly images).
- **Source:** environment.docx

### G2. SAR Processing (Section 2.3)
- **What changed:** Added IW mode, VV+VH polarisations, descending orbit filter details. Added scene count: 913 Sentinel-1 scenes. Replaced Lee Sigma filter with focal mean filter (7x7 kernel) matching actual GEE implementation. Replaced generic processing steps with four specific steps: (1) focal mean speckle filter, (2) VV/VH ratio and normalised difference bands for water/flooded vegetation discrimination, (3) adaptive VV threshold at -15 dB for water detection, (4) 96 monthly composites (2015-2024).
- **Source:** methods-sar.docx

### G3. Landsat Processing (Section 2.3)
- **What changed:** Added band harmonisation to common six-band schema across L5/L7/L8. Added Collection 2 scale factors (DN x 0.0000275 - 0.2). Added QA_PIXEL cloud masking and 30% cloud cover filter. Added scene count: 1,335 harmonised scenes. Added 39 annual median composites (1984-2024). Replaced "inverted elastic-net / dark spectrum fitting" paragraph with Otsu adaptive thresholding description (per Pekel et al., 2016), including specific MNDWI percentile values (p10=-0.57, p50=-0.51, p90=0.74) from 2020 dry-season composite.
- **Source:** methods-landsat.docx

### G4. Training Samples and Classification (Section 2.3)
- **What changed:** Added training sample generation method (threshold-based MNDWI/NDVI stratification refined by community validation). Added total sample count: 450 (150 water, 200 flooded vegetation, 180 dry vegetation, 120 bare soil, 80 urban/built). Added feature stack: 11 bands (6 spectral + 5 indices). Added 70/30 train/validation split with fixed random seed.
- **Source:** methods-landsat.docx

### G5. References section: 16 missing citations added
These references were cited in the manuscript text but had no corresponding entry in the reference list:
- Eva and Lambin (2000)
- Furse, Morgan, and Kalk (1979)
- Howard-Williams and Gaudet (1985)
- Howard-Williams and Walker (1974)
- Kalk and Schulten-Senden (1977)
- Kambombe et al. (2021)
- Kutser (2012)
- Matthews (2011)
- Ngongondo et al. (2011)
- Nicholson, Klotter, and Chavula (2014)
- Ogashawara, Mishra, and Gitelson (2017)
- Ozesmi and Bauer (2002)
- Palmer, Kutser, and Hunter (2015)
- Philippe and Karume (2019)
- Sulieman and Ahmed (2013)
- Vanden Bossche and Bernacsek (1990)

### G7. Results 3.1: SMA water fraction dynamics added
- **What changed:** Added sub-pixel SMA results paragraph reporting basin-mean open water fraction range (0.04 to 0.21), flooded vegetation fraction range (0.05 to 0.14), and the inverse relationship between the two. Added reference to 39 annual composites and 1,335 scenes in the opening sentence.
- **Source:** methods-landsat.docx (SMA time series data table), results.docx

### G8. Results 3.2: Index comparison data integrated
- **What changed:** Added basin-wide annual mean tracking of recession-refilling cycles across all five indices. Noted MNDWI/NDPI correspondence and WRI divergence during turbid peaks. Explicitly stated NDWI as weakest under turbid/shallow conditions. Added "39-year time series" reference. Strengthened NDPI description to "tracked the water-vegetation boundary most effectively."
- **Source:** results.docx (index comparison narrative), methods-landsat.docx (index time series data)

### G9. Results 3.3: SAR monthly time series integrated
- **What changed:** Added 913 scenes / 96 monthly composites reference. Added new paragraph reporting monthly water fraction oscillation (0.04 in 2015 rising to 0.09-0.11 from 2018 onward), mean VV backscatter range (-7.2 to -12.0 dB), VH offset (~7 dB lower), scene density improvement over time, and January 2020 wet-season composite as contrast reference.
- **Source:** methods-sar.docx (monthly time series data table)

### G11. Results 3.1: SMA figure data integrated
- **What changed:** Added combined fraction stability (0.20-0.25 in non-recession years) and its sharp drop during 1995/2012 as basin-wide desiccation rather than lateral redistribution. Added "Figure X" placeholder.
- **Source:** SMA stacked area chart (rendered output)

### G12. Results 3.2: Correlation matrix and MNDWI time series values integrated
- **What changed:** Added full inter-index correlation matrix values: NDPI-MNDWI r=-1.000, WRI-MNDWI r=0.829, NDWI-AWEIsh r=0.759, NDWI-MNDWI r=0.448, AWEIsh-MNDWI r=0.557. Added interpretation: NDWI and MNDWI capture genuinely different spectral properties (moderate r=0.448 despite both being water indices). Added basin-wide MNDWI value range (-0.25 to -0.50) with recession dips. Added "Figure X" and "Table X" placeholders.
- **Source:** Correlation matrix table, Landsat MNDWI annual median time series plot (rendered outputs)

### G13. Results 3.3: SAR water fraction from figure data corrected
- **What changed:** Corrected wet-season peak water fraction from 0.09-0.11 (January-only values from data table) to 0.25-0.30 (actual seasonal peaks visible in rendered figure). Added dry-season trough at ~0.05. Added five- to six-fold seasonal amplitude. Added regularity description (2016-2024). Added "Figure X" placeholder.
- **Source:** SAR-derived water fraction monthly time series plot (rendered output)

### G14. Notes for future consolidation
- **Lancaster (1979)** is cited in Results 3.1 but has no reference entry. May be an error for Kalk (1979). Needs author verification.
- **Donda, Njaya, and Bene (2012)** is cited in Discussion 4.4 but the reference entry has Njaya as first author: "Njaya, F., Donda, S., & Bene, C. (2012)." Verify whether these are the same or different papers.
- **conclusions.docx** contained only "To be developed" and was not integrated.
- **discussion.docx** contained generic discussion text less specific than the existing manuscript; not integrated.
- **Accuracy figures:** The GEE workflow reports 99.3% overall accuracy / 0.989 kappa, but these derive from threshold-generated training samples rather than independent validation. The manuscript retains the more conservative 81% / 0.77 figures, which should be verified against the final validated results.
