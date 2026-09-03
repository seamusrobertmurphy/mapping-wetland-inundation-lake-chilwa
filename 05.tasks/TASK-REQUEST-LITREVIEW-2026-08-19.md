# Task request, literature review

Written 2026-08-19. Applies the standard now set in the global `CLAUDE.md`: every manuscript carries
a literature review as a section before the Introduction, written against that manuscript's own
objectives, with numbered findings cited by number in the Introduction, the hypotheses, the Methods,
the Results, the Discussion and the Conclusions. The worked example is
`claude-science-library/publications-academic/beetle-topography-and-wind-study`.

## Manuscript

`01.manuscript/Manuscript_2026-08-03.qmd`

> Mapping wetland inundation, Lake Chilwa

## Current state, measured

Counted on 2026-08-19 with `scientific-library/library/tools/cite-audit.py`, which is used rather
than an ad-hoc grep because a hand-rolled pattern reported zero citations on three manuscripts in
this set that in fact carry dozens. Citation keys here are `Author_Year`, and a pattern assuming
lowercase-then-digits misses every one of them.

```
/Volumes/PortableSSD/Github/mapping-wetland-inundation-lake-chilwa/01.manuscript/Manuscript_2026-08-03.qmd
   citations: 0 total, 0 unique
   bib: references.bib  entries=450  uncited=450
```

**Defect to fix first.** A large majority of the bibliography is uncited, which
usually means reading was done and never connected to the text. The review is the
mechanism for connecting it.

## Objectives the review must be written against

1. To be read from the manuscript. It currently carries no citations at all against a 450-entry bibliography, so the objectives must be restated before the corpus is screened.

## Corpus to screen

Inundation and flood mapping from optical and SAR time series; endorheic lake water balance; Lake Chilwa hydrology and its recession cycles; livelihoods and fisheries dependence on inundation extent.

Seed the search from the central store at `Github/scientific-library`, which held 579 distinct papers
and 380 CrossRef-verified entries on 2026-08-19. Candidates already in the store include:

- Njaya 2011 natural history and fisheries ecology of Lake Chilwa
- Macuiane 2011 seasonal dynamics
- Nicholson 2013 rainfall climatology for Malawi
- Missi 2020 groundwater and surface water
- Adam 2009 wetland vegetation remote sensing review

Anything the review needs and the store lacks is added to the store, not to this repository alone.

## The claim to check hardest

Determine after the objectives are restated. The 450-entry bibliography suggests the reading was done and never connected to the text.

## Procedure

1. Restate the objectives from the manuscript, not from this file.
2. Screen the store on topic, then on whether a study can answer an objective. Record every decision
   and its reason in `02.inputs/literature/review-screening.csv`.
3. Read in full text every study that survives screening. An absence established from titles and
   abstracts is not evidence.
4. Write `02.inputs/literature/review-protocol.md`: search date, stores queried, both screening
   stages, the per-study extraction, and the limits of the review.
5. Add `# Literature review` before `# Introduction`, with an executable chunk that reads the
   screening table so every count is computed at render time and none is typed into prose.
6. Number the findings R1 onward and cite them by number in each later section. A review no later
   section cites has been appended, not integrated.
7. Reconcile the bibliography: every cited key present, every entry CrossRef-verified, uncited
   entries either used or removed.

## Do not

Do not generate a summary of what a paper shows and file it beside a verified one; they become
indistinguishable. Quote verbatim with a page number in a `sources/<key>.md` note in the central
store, and link to it. Do not report a count you did not compute. Do not claim an absence you
checked only by search: on 2026-08-19 a correctly cited claim was reported as fabricated because a
grep missed an `fl` ligature in the PDF text layer.
