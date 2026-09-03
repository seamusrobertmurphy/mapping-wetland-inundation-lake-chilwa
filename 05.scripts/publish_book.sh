#!/bin/zsh
# Publish the built book to GitHub Pages.
#
# `quarto publish gh-pages` cannot be used here, because it renders inside the
# repository and this volume is exFAT, where Quarto's filesystem calls fail on
# the AppleDouble shadow files macOS writes. So 05.scripts/build_book.sh renders
# on the internal disk and this script pushes that output to the gh-pages branch
# with plain git, which has no such trouble.
#
# The worktree also lives on the internal disk, so nothing copies onto exFAT.
#
# Usage:  ./05.scripts/build_book.sh && ./05.scripts/publish_book.sh

set -e
REPO="/Volumes/PortableSSD/Github/mapping-wetland-inundation-lake-chilwa"
BOOK="$HOME/chilwa-book/_book"
TREE="$HOME/chilwa-ghpages"
BRANCH="gh-pages"

[ -d "$BOOK" ] || { echo "No build at $BOOK. Run 05.scripts/build_book.sh first."; exit 1; }

cd "$REPO"
git worktree remove --force "$TREE" 2>/dev/null || true
rm -rf "$TREE"

if git show-ref --verify --quiet "refs/heads/$BRANCH"; then
  git worktree add "$TREE" "$BRANCH" -q
else
  # First publish. An orphan branch keeps the site's history separate from the
  # manuscript's, so the built pages never clutter the source log.
  git worktree add --detach "$TREE" -q
  git -C "$TREE" checkout --orphan "$BRANCH" -q
  git -C "$TREE" rm -rf . -q 2>/dev/null || true
fi

find "$TREE" -mindepth 1 -not -path "$TREE/.git*" -delete 2>/dev/null || true
cp -R "$BOOK"/. "$TREE"/
# Without this, GitHub runs Jekyll and drops every directory beginning with an
# underscore, which is most of what Quarto writes.
touch "$TREE/.nojekyll"
find "$TREE" \( -name '._*' -o -name '.__*' -o -name '.DS_Store' \) -delete 2>/dev/null || true

cd "$TREE"
git add -A
if git diff --cached --quiet; then
  echo "No change since the last publish."
else
  git commit -q -m "Publish book, $(date -u '+%Y-%m-%d %H:%M UTC')"
fi
git push -q -u origin "$BRANCH"
cd "$REPO"
git worktree remove --force "$TREE"

echo "Pushed to $BRANCH."
echo "Site: https://seamusrobertmurphy.github.io/mapping-wetland-inundation-lake-chilwa/"
