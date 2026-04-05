# Project: Mapping Wetland Inundation, Lake Chilwa

Integrated remote sensing and community mapping of an endorheic wetland in southern Malawi. Combines multi-temporal satellite analysis (Sentinel-1, Landsat 1994-2015) with participatory mapping of fishing communities in the Lake Chilwa Basin.

## Writing Style

Write with the economy of good poetry so that every word is chosen for necessity and precision, and none are carried along by inertia. Prefer the concrete over the abstract, the specific over the general. When a single well-chosen word will do, do not reach for three approximate ones.

Favour plain, direct prose using short sentences where the logic is dense, longer ones only where rhythm demands it. Cut qualifications, hedges, and throat-clearing. Never open with preamble. Land the point first.

The standard is not brevity for its own sake, but distillation: language at its most compressed without losing precision or nuance. If the sentence can carry more weight with fewer words, it must.

Avoid AI-register entirely. This means no em-dashes used for dramatic pause or emphasis, no bullet-point decomposition of ideas that belong in continuous prose, no filler transitions. Do not use symbols in place of words. Colons and commas do the work that dashes and slashes pretend to. Parentheses are used sparingly and only where a true aside is needed, not as a substitute for clear sentence structure.

Punctuation should be conventional and minimal: periods, commas, colons, and the occasional semicolon. Nothing more unless the writing specifically calls for it. No asterisks for emphasis, no slashes as shorthand, no brackets outside of citations or technical notation. Write the way a careful human writer edits: reading aloud, cutting what sounds performed, keeping what sounds earned.

## Design Objective

The goal is prose that is simultaneously compressed and musical. Maximum meaning per word, but never at the cost of cadence. A well-compressed sentence naturally has good rhythm; density and flow are inseparable, and neither is sacrificed for the other.

The writer recedes. No inflated claims, no self-congratulatory framing. Present the work plainly and let the reader draw conclusions. Humility is not false modesty but the discipline of letting evidence speak without amplification.

Choose the word that captures meaning most quickly and accurately within the wider argument. The most suitable term is often the simplest one. "Use" not "utilise", "show" not "demonstrate", "because" not "due to the fact that." Reach for a more specific or technical word only when the plain one is genuinely ambiguous or loses a distinction the argument requires. In technical passages (Methods, Results), precision wins when simplicity and precision diverge. In narrative passages (Introduction, Discussion), simplicity wins.

## Scientific Voice

Prefer active voice. Passive is tolerated where it reads naturally, particularly in Methods, but active is the default. Use first person plural ("we") where appropriate.

Hedge minimally. State findings directly. Reserve hedging language for claims where genuine uncertainty exists. Write "the results show" not "the results seem to suggest." Let the data carry the weight; do not soften what the evidence supports.

## Syntax and Grammar

Default to short, declarative sentences. Compound and complex sentences are used sparingly and only when the relationship between ideas demands it. Keep the grammatical core visible.

Always use the Oxford comma. Prefer "that" for restrictive clauses, "which" for non-restrictive, and cut either when the sentence works without it. Do not use contractions in formal prose: "do not" not "don't", "cannot" not "can't." No dangling modifiers: participial phrases must attach to their correct subject.

## Rhetoric and Transitions

State the point, then support it. The reader should know where you stand before they see why. Never build to a reveal.

Use standard connectives but keep them short: one word, not a phrase. "Yet" not "on the other hand." "However" earns its place; "furthermore" and "moreover" almost never do.

## Structure and Flow

Each paragraph opens with its claim. Supporting evidence and elaboration follow. One idea per paragraph, no paragraph exceeding roughly 150 words.

Never use bullet points or numbered lists in manuscript text. Enumerate within sentences: "first... second... third..." or separate items with semicolons. Continuous prose is the only acceptable form for argumentation, description, and reporting.

## Terminology

British English throughout: favour, analyse, centre, colour.

Define specialist terms on first use with a brief explanation, then use freely without further definition. Assume a literate but not specialist reader: someone comfortable with scientific prose but not necessarily trained in remote sensing or wetland ecology.

## Working Rules

Never edit master files directly. All code, manuscript drafts, and scripts must be written to `_staging/` as duplicate copies. The user reviews and approves all changes before they are implemented into the main project tree. This protects active drafts from corruption.

## Project Structure

- `1.manuscript/` — manuscript drafts (DOCX)
- `2.data/` — shapefiles (basin, sub-basins, AOI) and rasters (DEM derivatives)
- `3.literature/` — reference PDFs
- `4.images/` — figures
- `mapping-wetland-inundation-lake-chilwa.qmd` — Quarto analysis document
- `mapping-wetland-inundation-lake-chilwa.Rmd` — R Markdown analysis document
- `references/` — citation files
- `apa.csl` — APA citation style
