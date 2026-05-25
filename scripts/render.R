#!/usr/bin/env Rscript
# Render an Rmd from the command line, honoring its YAML `knit:` hook.
# Usage (from project root):  Rscript scripts/render.R paper.Rmd

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Usage: Rscript scripts/render.R <file.Rmd>")
input <- args[1]

fm <- rmarkdown::yaml_front_matter(input)
if (!is.null(fm$knit)) {
  eval(parse(text = fm$knit))(input)
} else {
  rmarkdown::render(input)
}

# Clean up xelatex/knitr build artefacts left alongside the input file.
# Artefact names follow the YAML `output_file` stem (if set), otherwise the input stem.
input_dir <- dirname(input)
stems <- tools::file_path_sans_ext(basename(input))
if (!is.null(fm$output_file)) {
  stems <- c(stems, tools::file_path_sans_ext(basename(fm$output_file)))
}
# Also pull output_file from a `knit:` hook string, e.g. output_file = "paper.pdf"
if (!is.null(fm$knit)) {
  m <- regmatches(fm$knit, regexec('output_file\\s*=\\s*["\']([^"\']+)["\']', fm$knit))[[1]]
  if (length(m) >= 2) stems <- c(stems, tools::file_path_sans_ext(basename(m[2])))
}
stems <- unique(stems)
for (stem in stems) {
  path_stem <- file.path(input_dir, stem)
  for (ext in c("log", "tex", "knit.md", "aux", "out", "toc", "nav", "snm", "vrb")) {
    file.remove(Sys.glob(paste0(path_stem, ".", ext)))
  }
  unlink(paste0(path_stem, "_files"), recursive = TRUE)
  unlink(paste0(path_stem, "_cache"), recursive = TRUE)
}
