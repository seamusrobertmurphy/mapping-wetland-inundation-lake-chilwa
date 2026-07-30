#!/usr/bin/env python3
"""Guard against the index drifting from the repository.

Checks that every path named in SPEC.md, DECISIONS.md, and INDEX.md exists, and
that no governing file states D8 as a method in use. Run from anywhere:

    python3 05.scripts/check_repo_index.py
"""

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GOVERNING = ["SPEC.md", "DECISIONS.md", "INDEX.md", "CLAUDE.md"]

# Backtick-quoted tokens that look like repo paths.
PATH_RE = re.compile(r"`([0-9A-Za-z_][0-9A-Za-z_./-]*\.[A-Za-z0-9]{1,6}|[0-9A-Za-z_][0-9A-Za-z_./-]*/)`")

# D8 is retired. These phrasings retire it; anything else naming D8 is suspect.
RETIREMENT_RE = re.compile(
    r"retired|not (?:appear|be used)|never|rather than|instead of|comparator|"
    r"forbid|excluded|must not|no longer|defect|open item|D8 pointer|disclosed",
    re.I,
)

# A line may legitimately name a path in order to say it is absent.
ABSENCE_RE = re.compile(
    r"there is no|does not exist|no longer exists|never (?:existed|supplied)|"
    r"absent|missing|has never been|slated for deletion",
    re.I,
)


def repo_dir_names():
    """Directory names anywhere in the repo, so context-relative shorthand
    such as `PDF/` under 03.outputs resolves without a full path."""
    names = set()
    for d in ROOT.rglob("*"):
        if d.is_dir() and ".git" not in d.parts:
            names.add(d.name)
    return names


def check_paths():
    problems = []
    dir_names = repo_dir_names()
    for name in GOVERNING:
        f = ROOT / name
        if not f.exists():
            problems.append(f"{name}: governing file missing")
            continue
        for lineno, line in enumerate(f.read_text(encoding="utf-8").splitlines(), 1):
            for token in PATH_RE.findall(line):
                if token.startswith(("http", "@")) or " " in token:
                    continue
                # Globs and wildcards are descriptive, not literal paths.
                if any(c in token for c in "*?"):
                    continue
                # Only repo-relative paths are checkable. A bare filename in
                # running text is a reference, not a location. A trailing
                # slash always means a directory, so check those too.
                if not token.endswith("/") and "/" not in token:
                    continue
                if (ROOT / token).exists():
                    continue
                # Context-relative directory shorthand, e.g. `PDF/` under 03.outputs.
                if token.endswith("/") and token.rstrip("/") in dir_names:
                    continue
                if ABSENCE_RE.search(line):
                    continue
                problems.append(f"{name}:{lineno}: path does not exist: {token}")
    return problems


def check_d8():
    problems = []
    for name in GOVERNING:
        f = ROOT / name
        if not f.exists():
            continue
        for lineno, line in enumerate(f.read_text(encoding="utf-8").splitlines(), 1):
            if "D8" in line and not RETIREMENT_RE.search(line):
                problems.append(f"{name}:{lineno}: D8 named without retirement context: {line.strip()[:90]}")
    return problems


def main():
    problems = check_paths() + check_d8()
    if problems:
        print(f"FAIL: {len(problems)} problem(s)\n")
        for p in problems:
            print("  " + p)
        return 1
    print("OK: governing files match the repository, and D8 appears only as retired.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
