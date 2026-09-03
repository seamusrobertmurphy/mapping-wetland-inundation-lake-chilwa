#!/bin/zsh
# Build the Quarto book from the master manuscript, then render it.
#
# Two things about this machine shape the script.
#
# The repository lives on /Volumes/PortableSSD, which is exFAT. macOS cannot
# store extended attributes there, so it writes an AppleDouble shadow beside
# every file it touches: ._name for a file, .__name for a dotfile, plus
# __MACOSX folders from any unzip. Quarto's filesystem calls trip over them,
# and two separate failures were traced to this, a directory walk of
# ._execute-results reported as "Not a directory", and a rename of index.html
# into _book reported as "No such file or directory". So the render happens in
# ~/chilwa-book, which is APFS, and the shadows are swept before and after.
#
# The chapters are generated into 01.manuscript/Archive/book/, so nothing there
# is source, and they are staged as chapters/ for the render. The master at
# 01.manuscript/Manuscript_2026-08-03.qmd is the only file to edit.
#
# Usage:  ./05.scripts/build_book.sh          frozen, nothing executes
#         ./05.scripts/build_book.sh --live   keep the master's own eval settings

set -e
REPO="/Volumes/PortableSSD/Github/mapping-wetland-inundation-lake-chilwa"
STAGE="$HOME/chilwa-book"

sweep() {  # remove every macOS metadata artefact under a directory
  find "$1" \( -name '._*' -o -name '.__*' -o -name '__MACOSX' -o -name '.DS_Store' \) \
       -not -path '*/.git/*' -exec rm -rf {} + 2>/dev/null || true
}

cd "$REPO"
sweep "$REPO"
/usr/local/bin/Rscript 05.scripts/split_manuscript.R "$@"

rm -rf "$STAGE"
mkdir -p "$STAGE/chapters" "$STAGE/03.outputs/PNG" "$STAGE/04.references" "$STAGE/05.scripts"

cp 05.scripts/book/_quarto.yml "$STAGE/_quarto.yml"
cp 05.scripts/book/styles.css "$STAGE/styles.css"
cp 01.manuscript/Archive/book/index.qmd "$STAGE/index.qmd"
cp 01.manuscript/Archive/book/0[1-5]-*.qmd "$STAGE/chapters/"
cp 04.references/references.bib 04.references/apa.csl "$STAGE/04.references/"
cp 05.scripts/_common.R        "$STAGE/05.scripts/"

# Only the figures the chapters actually cite, read from the chapters themselves
# so a renamed figure is reported rather than silently dropped.
grep -oh '\.\./03\.outputs/PNG/[A-Za-z0-9_.-]*' "$STAGE"/chapters/*.qmd 2>/dev/null \
  | sed 's|\.\./||' | sort -u | while read -r rel; do
      if [ -f "$REPO/$rel" ]; then cp "$REPO/$rel" "$STAGE/03.outputs/PNG/"
      else echo "  WARNING: cited figure missing, $rel"; fi
    done

# The margin-header figure is named in _quarto.yml rather than in a chapter, so
# the citation scan above does not see it.
cp 03.outputs/PNG/fig01_study_area.png "$STAGE/03.outputs/PNG/"

sweep "$STAGE"
echo "staged $(du -sh "$STAGE" | cut -f1) in $STAGE"

cd "$STAGE"
quarto render
sweep "$STAGE"

echo
echo "Built: $STAGE/_book"
echo "Publish with:  cd $STAGE && quarto publish quarto-pub"
