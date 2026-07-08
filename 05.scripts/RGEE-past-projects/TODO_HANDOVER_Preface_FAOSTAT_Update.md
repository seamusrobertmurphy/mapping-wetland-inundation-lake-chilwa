# Handover: Preface FAOSTAT Section Update

**Date:** 2026-04-04
**Status:** Draft edits applied to `index.qmd` but user requests they be moved to a working branch or reverted. All source materials and data are staged.

---

## 1. What Was Done

### FAOSTAT Bulk Data Downloaded

Honduras fire emissions data extracted from the FAOSTAT Emissions from Fires domain bulk download (`Emissions_Land_Use_Fires_E_All_Data.zip` from `bulks-faostat.fao.org`). The clean long-format CSV is at:

```
data/raw/honduras_faostat_fires.csv
```

Contains 1,612 records covering 1990-2024 with these fields:

| Field | Description |
|---|---|
| country | "Honduras" |
| year | 1990-2024 |
| item_code | FAOSTAT item code (6760-6993) |
| item | Vegetation subcategory name |
| element_code | 7225 (CH4), 7230 (N2O), 7245 (biomass), 7246 (burned area), 7273 (CO2) |
| element | Human-readable element name |
| unit | ha, t, or kt |
| value | Numeric value |

### FAOSTAT Item Codes for Honduras

Twelve vegetation subcategories within three fire categories:

| Code | Item | Fire Category |
|---|---|---|
| 6796 | Humid tropical forest | Forest fires |
| 6797 | Other forest | Forest fires |
| 6992 | Forest fires (aggregate) | Forest fires |
| 6760 | Savanna | Savanna fires |
| 6789 | Woody savanna | Savanna fires |
| 6791 | Closed shrubland | Savanna fires |
| 6792 | Open shrubland | Savanna fires |
| 6794 | Grassland | Savanna fires |
| 6790 | Savanna and woody savanna (aggregate) | Savanna fires |
| 6793 | Closed and open shrubland (aggregate) | Savanna fires |
| 6795 | Savanna fires (aggregate) | Savanna fires |
| 6993 | Fires in organic soils | Organic soil fires |

### Honduras 2010 Reference Year Summary

| Fire Category | Burned Area (ha) | Biomass Burned (t DM) | CH4 (kt) | N2O (kt) |
|---|---:|---:|---:|---:|
| Forest fires | 41,131 | 2,212,389 | 5.09 | 0.46 |
|   Humid tropical forest | 37,666 | 2,037,754 | 4.69 | 0.43 |
|   Other forest | 3,465 | 174,635 | 0.40 | 0.04 |
| Savanna fires | 117,256 | 739,478 | 1.70 | 0.16 |
| Organic soil fires | 0 | 0 | 0 | 0 |

### Edits Applied to index.qmd

The following changes were made directly to `index.qmd`. The user may want to revert these and re-apply from a working branch:

1. **Package typo fix** (line 16): `"FAOTSTAT"` changed to `"FAOSTAT"`

2. **FAOSTAT section replaced** (lines 148-370): The old FAOSTAT + QA/QC sections (Database Overview, Three Fire Categories, Temporal Gap-Filling, QA/QC with Tiered Enhancements, Mitigation Tracking) were replaced with a restructured section that introduces each global stratification dataset as a term in Equation 2.27:

   - `### Activity Data: Burned Area (A)` — MCD64A1, gap-filling, tmap map
   - `### Land Cover Classification: Three Fire Categories` — MCD12Q1 reclassification, tmap map
   - `#### Organic Soil Screening` — HWSD/Histosols, callout note updated to "Central and South America"
   - `### Fuel Consumption Values (MB × Cf)` — Table 2.4 (tropical subset), FAO GEZ tmap map
   - `### Combustion Factors (Cf)` — Table 2.6 (tropical subset)
   - `### Emission Factors (Gef)` — Table 2.5 (all five fire types), AR5 GWP-100 values
   - `### QA/QC and Validation` — condensed to two paragraphs (tiered enhancements and mitigation tracking cut)
   - `### Honduras: FAOSTAT Benchmark` — loads `data/raw/honduras_faostat_fires.csv`, creates burned area stacked area plot, CH4/N2O time series, 2010 summary table

3. **Sections cut:**
   - Tiered Enhancements (Tier 1→2→3 progression)
   - Mitigation Action Tracking

4. **Sections preserved unchanged:**
   - Introduction paragraphs and Guide Outline
   - IPCC Tier 1 Fire Emissions (Equation 2.27, vegetation applications, MLP)
   - Uncertainty Considerations
   - Data Sources
   - Footnotes

### tmap Visualization Pattern

All maps use the same pattern established in Chapter 3:

```r
ee_tile_url <- function(ee_image, vis_params) {
  ee_image$getMapId(vis_params)$tile_fetcher$url_format
}

tmap::tmap_mode("view")
tmap::tm_shape(aoi_sf) + tmap::tm_borders(col = "white", lwd = 1.5) +
  tmap::tm_tiles(ee_tile_url(image, vis_params), group = "Label") +
  tmap::tm_basemap("Esri.WorldImagery")
```

GEE is initialized via reticulate (not rgee) matching the chapter pattern:

```r
library(reticulate)
use_python("/opt/local/bin/python", required = TRUE)
ee <- import("ee")
ee$Initialize(project = "murphys-deforisk")
```

### IPCC Tables Included

Three tables from 2019 Refinement Vol 4 Ch 2, filtered to Honduras-relevant tropical rows:

- **Table 2.4** — Fuel consumption (t DM ha⁻¹): primary/secondary tropical forest, shrublands, savanna woodlands/grasslands (early/late season), peatland
- **Table 2.5** — Emission factors (g kg⁻¹): all five fire types (savanna, ag residues, tropical forest, extra-tropical, biofuel)
- **Table 2.6** — Combustion factors: same vegetation types as Table 2.4

All tables use Quarto cross-referenceable captions with `{.striped .hover}` styling.

---

## 2. Scientific Writer Plugin

### Installed

The `scientific-writer` Python package (v2.12.1) was installed via pip. The git repo was cloned to `/tmp/claude-scientific-writer/`.

Nine relevant skills were copied to the project at:

```
.claude/skills/
├── citation-management/
├── document-skills/
├── literature-review/
├── peer-review/
├── research-grants/
├── research-lookup/
├── scientific-critical-thinking/
├── scientific-writing/
└── venue-templates/
```

### Not Done

The plugin was NOT installed via Claude Code `/plugin` commands (Cowork mode cannot execute slash commands). To complete the installation in Claude Code:

```bash
# Option 1: Install from marketplace
/plugin marketplace add https://github.com/K-Dense-AI/claude-scientific-writer
/plugin install claude-scientific-writer

# Option 2: The skills are already in .claude/skills/ and should be
# auto-discovered by Claude Code at session start.
```

The `CLAUDE.scientific-writer.md` template was copied to `.claude/` for reference.

---

## 3. What Remains

1. **Decide whether to keep or revert the index.qmd edits.** If reverting, the new content can be re-applied from this handover on a working branch (`todo/preface-faostat-update`).

2. **Test the tmap + GEE code blocks render.** The code follows established chapter patterns but has not been executed against a live GEE session. Key risks:
   - `FAO/GEZ/2010/V2` asset ID may need verification
   - `ee_to_sf()` for large FeatureCollections may timeout
   - tmap v4 API changes (the chapters use `tm_tiles` which is tmap v4 syntax)

3. **Honduras-specific vegetation table.** The user's original request referenced an Ecuador-style table with country-specific vegetation types (e.g., "Dry Andean Forest") mapped to local fuel loads. The current draft uses IPCC Tier 1 defaults only. To add Honduras-specific names and values (Pinus caribaea, broadleaf humid forest, dry pine-oak), run the GEE overlay of MCD12Q1 × GEZ within the Honduras boundary to identify which strata actually occur, then map them to the hypothetical Tier 2 values already in Chapter 3 (ICF, SIGMOF, La Mosquitia sources).

4. **FAOSTAT R package integration.** The `FAOSTAT` package is loaded in the packages list but the data is currently loaded from the bulk CSV. The package's `getFAO()` function could replace the bulk download for reproducibility, but the 2014-era API may need updating. Test with:
   ```r
   library(FAOSTAT)
   FAOsearch()
   ```

5. **Complete plugin installation** via Claude Code slash commands (see section 2 above).

6. **AR5 GWP-100 values confirmed:** CH₄ = 28, N₂O = 265 (per user decision).
