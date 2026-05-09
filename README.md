# Liver Fat and Chemotherapy Toxicity

Survival analysis of the association between liver radiodensity (Hounsfield Units, measured on diagnostic CT) and time to dose-limiting toxicity in colon-cancer patients receiving capecitabine (CAP), oxaliplatin (OXA), or both (CAPOX). Data from the Dutch multi-centre COLON cohort. Hospital-level heterogeneity is addressed via stratified Cox proportional-hazards models.

## Structure

```
liverfat_survival_paper.Rmd     – formal PDF paper (xelatex)
liverfat_survival_slides.Rmd    – Beamer presentation (metropolis, 16:9)
R/                              – supplementary analysis scripts
data/                           – dataset (not committed)
docs/                           – study documentation, codebook, feedback
reports/                        – knitted output PDFs
archive/                        – superseded Rmd versions
```

## Canonical Rmd files

- **`liverfat_survival_paper.Rmd`** — written report. Knits to PDF via xelatex.
- **`liverfat_survival_slides.Rmd`** — slide deck. Knits to `reports/liverfat_survival_slides.pdf`.

## Archive

Older iterations are kept in `archive/` for traceability:

- `liverfat_survival_paper_v1.Rmd` — first paper draft
- `liverfat_survival_slides_v1.Rmd` — first slide deck (Madrid theme)
- `liverfat_survival_slides_v2.Rmd` — second slide iteration (Madrid, dual-author)
- `liverfat_survival_report.Rmd` — abandoned HTML walkthrough variant

## Running the analysis

1. Place `Dataset_liverfat.RData` in `data/`
2. Open `liverfat_survival_analysis.Rproj` in RStudio
3. Knit `liverfat_survival_paper.Rmd` or `liverfat_survival_slides.Rmd`

Required packages: `survival`, `survminer`, `dplyr`, `knitr`, `broom`, `ggplot2`, `forestmodel`, `kableExtra`, `tibble`. Missing packages are auto-installed by the Rmd setup chunk.
