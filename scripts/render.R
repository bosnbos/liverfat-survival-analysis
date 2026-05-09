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

# Clean up xelatex/knitr build artefacts left in the working directory.
stem <- tools::file_path_sans_ext(input)
for (ext in c("log", "tex", "knit.md", "aux", "out", "toc", "nav", "snm", "vrb")) {
  file.remove(Sys.glob(paste0(stem, ".", ext)))
}
unlink(paste0(stem, "_files"), recursive = TRUE)
unlink(paste0(stem, "_cache"), recursive = TRUE)
