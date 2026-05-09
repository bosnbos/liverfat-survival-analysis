# Liver Fat and Chemotherapy Toxicity

Survival analysis of the association between liver radiodensity (Hounsfield Units, measured on diagnostic CT) and time to dose-limiting toxicity in colon-cancer patients receiving capecitabine (CAP), oxaliplatin (OXA), or both (CAPOX). Data from the Dutch multi-centre COLON cohort. Hospital-level heterogeneity is addressed via stratified Cox proportional-hazards models.

## Structure

```
R/main_analysis.R   – CANONICAL analysis (data prep + Cox fits + PH tests)
analysis.Rmd        – analyst's notebook (raw output, sensitivities, feedback log)
paper.Rmd           – formal PDF paper (xelatex)
slides.Rmd          – Beamer presentation (metropolis, 16:9)
data/               – dataset (not committed)
docs/               – study documentation, codebook, feedback
reports/            – knitted output PDFs
archive/            – superseded Rmd versions
```

## How the pieces fit together

`R/main_analysis.R` is the single source of truth for the canonical models. The three Rmds all `source()` it from their setup chunks, so there is no duplicated data prep or model fitting:

```
        R/main_analysis.R   ← edit here to change the analysis
               │
   ┌───────────┼───────────┐
   ▼           ▼           ▼
analysis.Rmd  paper.Rmd  slides.Rmd
(notebook)    (PDF)      (Beamer)
```

**Rule:** keep `R/main_analysis.R` deterministic. Experimental code, alternative specifications, and sensitivity checks belong in `analysis.Rmd` — promote them only after they're decided.

## Canonical files

- **`R/main_analysis.R`** — the analysis itself. Edit here.
- **`analysis.Rmd`** — analyst's notebook: raw model output, sensitivity analyses, feedback responses, decision log. Knits to HTML.
- **`paper.Rmd`** — written report. Knits to PDF via xelatex.
- **`slides.Rmd`** — slide deck. Knits to `reports/slides.pdf`.

## Archive

Older iterations are kept in `archive/` for traceability:

- `paper_v1.Rmd` — first paper draft
- `slides_v1.Rmd` — first slide deck (Madrid theme)
- `slides_v2.Rmd` — second slide iteration (Madrid, dual-author)
- `report.Rmd` — abandoned HTML walkthrough variant

## Running the analysis

1. Place `Dataset_liverfat.RData` in `data/`
2. Open `liverfat_survival_analysis.Rproj` in RStudio
3. Knit any of: `analysis.Rmd`, `paper.Rmd`, `slides.Rmd`

Each Rmd sources `R/main_analysis.R` automatically — no manual setup needed.

Required packages: `survival`, `survminer`, `dplyr`, `knitr`, `broom`, `ggplot2`, `forestmodel`, `kableExtra`, `tibble`. Missing packages are auto-installed by the Rmd setup chunks.
