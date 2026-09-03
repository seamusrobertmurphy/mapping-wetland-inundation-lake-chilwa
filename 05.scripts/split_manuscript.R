#!/usr/bin/env Rscript
# Generate the Quarto book chapters from the master manuscript.
#
# The master at 01.manuscript/Manuscript_2026-08-03.qmd stays the single
# canonical file. The chapters this writes are build artefacts, regenerated on
# every render, exactly as the DOCX is, so they cannot drift away from it.
# Nothing here edits the master.
#
# The master carries five level-two headings and nothing above them. Each
# becomes a chapter, and every heading is promoted one level so a chapter title
# sits at level one, as a book expects. Headings inside code chunks are left
# alone, since a comment beginning with a hash is not a heading.
#
# Run from the repository root:  Rscript 05.scripts/split_manuscript.R

master  <- "01.manuscript/Manuscript_2026-08-03.qmd"
out_dir <- "01.manuscript/Archive/book"   # generated pages, git-ignored

stopifnot(file.exists(master))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
lines <- readLines(master, warn = FALSE)

# ---- split the YAML header from the body -----------------------------------
fences   <- which(trimws(lines) == "---")
yaml_txt <- lines[(fences[1] + 1):(fences[2] - 1)]
body     <- lines[(fences[2] + 1):length(lines)]

yaml_value <- function(key) {
  # Pull one scalar or folded block out of the header without a YAML parser,
  # which would mangle the folded abstract and the smart quotes in the title.
  i <- grep(paste0("^", key, "\\s*:"), yaml_txt)
  if (length(i) == 0) return("")
  first <- sub(paste0("^", key, "\\s*:\\s*>?\\s*"), "", yaml_txt[i[1]])
  if (nzchar(trimws(first))) return(trimws(first))
  j <- i[1] + 1
  acc <- character()
  while (j <= length(yaml_txt) && grepl("^\\s+\\S", yaml_txt[j])) {
    acc <- c(acc, trimws(yaml_txt[j])); j <- j + 1
  }
  paste(acc, collapse = " ")
}

title    <- gsub('^"|"$', "", yaml_value("title"))
subtitle <- gsub('^"|"$', "", yaml_value("subtitle"))
abstract <- yaml_value("abstract")

# ---- mark which body lines sit inside a code chunk --------------------------
in_chunk <- logical(length(body))
open <- FALSE
for (i in seq_along(body)) {
  if (grepl("^```\\{", body[i])) open <- TRUE
  in_chunk[i] <- open
  if (open && trimws(body[i]) == "```") open <- FALSE
}

# ---- find the chapter boundaries -------------------------------------------
is_h2  <- grepl("^## [^#]", body) & !in_chunk
starts <- which(is_h2)
stopifnot(length(starts) == 5)
names_ <- trimws(sub("^##\\s*", "", body[starts]))
files  <- sprintf("%02d-%s.qmd", seq_along(names_), tolower(names_))
ends   <- c(starts[-1] - 1, length(body))

# The `here` package treats a _quarto.yml as a project root, so the book config
# in 01.manuscript/ would otherwise pull here() into that folder and every
# here("03.outputs", ...) call in the manuscript would resolve one level too
# deep. here::i_am() states where this file sits relative to the true root and
# pins it back to the repository.
setup_chunk <- function(chapter_file) c(
  "```{r}",
  "#| label: setup",
  "#| include: false",
  paste0("here::i_am(\"chapters/", chapter_file, "\")"),
  "source(here::here('05.scripts', '_common.R'))",
  "```",
  ""
)

# Chunks in the master carry their own `#| eval: true`, which beats the
# project-level default. For a frozen render those have to be turned off in the
# generated copy; the master is never touched.
frozen <- !("--live" %in% commandArgs(trailingOnly = TRUE))

for (k in seq_along(starts)) {
  block <- body[starts[k]:ends[k]]
  inner <- in_chunk[starts[k]:ends[k]]
  # Promote every heading by one level, outside code chunks only.
  h <- grepl("^#{2,6} ", block) & !inner
  block[h] <- sub("^#", "", block[h])
  # The chapter's own heading is now the YAML title, so the promoted copy would
  # print the name twice on the page and nest it under itself in the sidebar.
  block <- block[-1]
  has_code <- any(grepl("^```\\{", block))
  if (frozen) block[block == "#| eval: true"] <- "#| eval: false"
  out <- c(
    "---",
    paste0("title: \"", names_[k], "\""),
    "---",
    "",
    if (has_code) setup_chunk(files[k]),
    block
  )
  writeLines(out, file.path(out_dir, files[k]))
  message(sprintf("  %-22s %5d lines%s", files[k], length(out),
                  if (has_code) "  (sources _common.R)" else ""))
}

# ---- the landing page, also generated so it cannot go stale -----------------
writeLines(c(
  "---",
  "title: \"Abstract\"",
  "---",
  "",
  # No title or subtitle here. Quarto prints both from the book metadata at the
  # top of this page, and repeating them produced a second copy in the body and
  # a third in the right-hand table of contents.
  abstract,
  "",
  "",
  "## About this book",
  "",
  paste("Every page here is generated from the master manuscript at",
        paste0("`", master, "`"), "by `05.scripts/split_manuscript.R`.",
        "Edit the master, never these pages."),
  ""
), file.path(out_dir, "index.qmd"))
message("  index.qmd              generated")
message("\nchapters written to ", out_dir, "/  (git-ignored, rebuilt every time)")
message(if (frozen)
  "mode: FROZEN, every chunk eval: false. Re-run with --live once rgee is up."
  else "mode: LIVE, the master's own eval settings are kept.")
