
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Installation

``` r
devtools::install_github('mcbg/AptamerTools')
```

# AptamerTools

Contrived example since the data is already normalized using ANML.

``` r
library(magrittr)
library(data.table)
#> Warning: pakke 'data.table' blev bygget under R version 4.2.3
library(knitr)
#> Warning: pakke 'knitr' blev bygget under R version 4.2.3
library(AptamerTools)

temp <- tempfile()
url <- 'https://github.com/SomaLogic/SomaLogic-Data/raw/master/example_data.adat'
download.file(url, destfile = temp)

ds <- SomaDataIO::read_adat(temp)

# extract metadata
meta <- ds %>% attributes() %>% .$Col.Meta %>% data.table()
meta[, variable := gsub('-', '.', SeqId) %>% paste0('seq.', .)]

# convert to data.table
setDT(ds)
```

Somalogic provides data in a wide format. Here we transform to a long
format.

``` r
seq_cols <- names(ds) %>% { .[grepl('^seq', .)] }

ds.long <- ds %>% 
  melt(measure.vars = seq_cols) %>% 
  merge(meta, by = 'variable')
```

# Normalization example

The following code derives the scale factors of each step in the
pipeline and the normalized data. The variable `value_HIPCA` contains
the fully normalized measurement. Note that the ANML step isn’t
recommended by the author.

``` r
# derive references
ds.long[SampleType %in% c('Sample', 'QC'), ref.sample := median(value), .(variable)]
ds.long[SampleType %in% c('Sample'), mad := mad(value), .(variable)]
ds.long[, ref.hyb := median(value), .(variable)]
ds.long[SampleType == 'Calibrator', ref.cal := median(value), .(variable)]

# apply normalization
normalization_pipeline(
  ds.long,
  'value',
  steps = list(hyb_norm, intraplate_norm, plate_norm, calibration_norm, anml_norm),
  step_names = c('H', 'HI', 'HIP', 'HIPC', 'HIPCA'),
  refs = list('ref.hyb', derive_intraplate_ref, 'ref.cal', 'ref.cal', 'ref.sample')
)
```

``` r
ds.long[, .(SampleId, variable, value_HIPCA)] %>% 
  head() %>% 
  kable()
```

| SampleId | variable     | value_HIPCA |
|:---------|:-------------|------------:|
| 1        | seq.10000.28 |    513.1384 |
| 2        | seq.10000.28 |    515.2335 |
| 3        | seq.10000.28 |    459.2683 |
| 4        | seq.10000.28 |    467.9756 |
| 5        | seq.10000.28 |    445.0118 |
| 6        | seq.10000.28 |    462.6189 |

# Analysis

In the paper we recommend either robust regression or non-parametric
methods to analyse SomaScan data.

## Robust regression using MM-estimation

Uses the package `RobStatTM`.

``` r
library(RobStatTM)
#> Warning: pakke 'RobStatTM' blev bygget under R version 4.2.3
#> 
#> Vedhæfter pakke: 'RobStatTM'
#> Det følgende objekt er maskeret fra 'package:datasets':
#> 
#>     stackloss

analysis.MM.estimation <- ds.long[, {
    ml <- lmrobdetMM(log2(value_HIPCA) ~ Sex, data = .SD)
    coef <- ml %>% summary() %>% .$coefficients
    pvalue <- coef %>% .[2, 4]
    logFC <- coef %>% .[2, 1]
    converged <- ml$converged
    data.table(pvalue, logFC, FC = 2^logFC, converged)
}, .(variable)]

head(analysis.MM.estimation) %>% kable()
```

| variable     |    pvalue |      logFC |        FC | converged |
|:-------------|----------:|-----------:|----------:|:----------|
| seq.10000.28 | 0.7062293 |  0.0094066 | 1.0065415 | TRUE      |
| seq.10001.7  | 0.4637004 | -0.0664417 | 0.9549905 | TRUE      |
| seq.10003.15 | 0.4601715 |  0.0298400 | 1.0208989 | TRUE      |
| seq.10006.25 | 0.4195581 |  0.0265707 | 1.0185880 | TRUE      |
| seq.10008.43 | 0.6066472 | -0.0238581 | 0.9835988 | TRUE      |
| seq.10011.65 | 0.9753342 | -0.0014553 | 0.9989918 | TRUE      |

## Robust regression using M-estimation

Uses the package `RobStatTM`.

``` r
# settings for regression
control <- lmrobM.control(bb = 0.5, efficiency = 0.95, family = 'bisquare', mscale_maxit = 1000, max.it = 1000)

analysis.M.estimation <- ds.long[, {
    ml <- lmrobM(log2(value_HIPCA) ~ Sex, data = .SD, control = control)
    coef <- ml %>% summary() %>% .$coefficients
    pvalue <- coef %>% .[2, 4]
    logFC <- coef %>% .[2, 1]
    converged <- ml$converged
    data.table(pvalue, logFC, FC = 2^logFC, converged)
}, .(variable)]

head(analysis.M.estimation) %>% kable()
```

| variable     |    pvalue |      logFC |        FC | converged |
|:-------------|----------:|-----------:|----------:|:----------|
| seq.10000.28 | 0.6595544 |  0.0111812 | 1.0077803 | TRUE      |
| seq.10001.7  | 0.4470061 | -0.0714593 | 0.9516749 | TRUE      |
| seq.10003.15 | 0.3776596 |  0.0357996 | 1.0251248 | TRUE      |
| seq.10006.25 | 0.5217764 |  0.0216700 | 1.0151339 | TRUE      |
| seq.10008.43 | 0.5599208 | -0.0260886 | 0.9820793 | TRUE      |
| seq.10011.65 | 0.6553905 | -0.0203214 | 0.9860130 | TRUE      |

## Kruskal-Wallis (non-parametric)

``` r
# settings for regression
control <- lmrobM.control(bb = 0.5, efficiency = 0.95, family = 'bisquare', mscale_maxit = 1000, max.it = 1000)

analysis.kruskal <- ds.long[, {
    ml <- kruskal.test(value_HIPCA ~ Sex, data = .SD)
    pvalue <- ml$p.value
    data.table(pvalue)
}, .(variable)]

head(analysis.kruskal) %>% kable()
```

| variable     |    pvalue |
|:-------------|----------:|
| seq.10000.28 | 0.2871925 |
| seq.10001.7  | 0.4368176 |
| seq.10003.15 | 0.6750879 |
| seq.10006.25 | 0.3888388 |
| seq.10008.43 | 0.4497916 |
| seq.10011.65 | 0.6147377 |
