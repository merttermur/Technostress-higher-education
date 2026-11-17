# Technostress Inhibitors, Technostress Creators and Work Performance in Higher Education

Analysis code for the study of technostress inhibitors, technostress creators and
self-reported work performance among university faculty members in Türkiye
(N = 408), submitted to *Scientific Reports*.

## Model

- Partial least squares structural equation modelling (PLS-SEM) with the R package `plspm`
- 9 reflective constructs (Mode A), 39 indicators, 23 structural paths
- Technostress inhibitors: literacy facilitation (LF), technical support provision (TSP),
  involvement facilitation (IF)
- Technostress creators: techno-overload (TO), techno-invasion (TI), techno-complexity (TC),
  techno-insecurity (TIS), techno-uncertainty (TU)
- Outcome: work performance (WP)
- 5,000 bootstrap resamples, 5,000 permutations, random seed = 42

## Files

| File | Content |
|---|---|
| `01_manuscript_tables.R` | Tables 1–7 of the manuscript: participant characteristics, reliability and convergent validity, Fornell–Larcker criterion, HTMT, structural model with bootstrap confidence intervals, university-type multigroup comparison with MICOM Step 2, cross-loadings |
| `02_peer_review_analyses.R` | Additional analyses: missing-data audit, HTMT inference (90% and 95% bootstrap intervals), item diagnostics, direct/indirect/total associations, inner VIF and f², sensitivity analysis without techno-invasion, MICOM Steps 1–3, first-component (Harman) diagnostic, ordinal EFA and CFA, exploratory age and gender comparisons, item-level sensitivity analysis without TO6 and TIS1 |

Run `01_manuscript_tables.R` first and `02_peer_review_analyses.R` afterwards in the same R session.

## How to run

### Google Colab
1. Open a new notebook and choose **Runtime → Change runtime type → R**.
2. Upload the data file `Anket.csv` to `/content/`.
3. Paste and run `01_manuscript_tables.R` in the first cell and `02_peer_review_analyses.R` in the second cell.

### Local R
Place `Anket.csv` in the working directory and run:

```r
source("01_manuscript_tables.R")
source("02_peer_review_analyses.R")
```

All results are printed to the console.

## Requirements

- R (≥ 4.3)
- `plspm` (installed from CRAN, or from GitHub `gastonstat/plspm` if CRAN installation fails)
- `psych`, `GPArotation`, `lavaan` (installed automatically by the second script if missing)

## Data

The survey data are not included in this repository. Data are available from the
corresponding author on reasonable request.

## License

MIT License (see `LICENSE`).
