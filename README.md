
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Sparse Correspondence Analysis

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/sparseCorrespondenceAnalysis)](https://CRAN.R-project.org/package=sparseCorrespondenceAnalysis)
<!-- badges: end -->

The goal of sparseCorrespondenceAnalysis is to illustrate Correspondence
Analysis and its sparsification to a data-set of the cause of deaths in
the United States in 2019.

## Installation

You can install the development version of
`sparseCorrespondenceAnalysis` from [GitHub](https://github.com/) with:

``` r
devtools::install_github("vguillemot/sparseCorrespondenceAnalysis")
```

## Running Sparse Correspondence Analysis

First load the package and the “Cause of death” data set.

``` r
library(sparseCorrespondenceAnalysis)
#> Le chargement a nécessité le package : PMA
#> Le chargement a nécessité le package : ggplot2
#> Le chargement a nécessité le package : ggrepel
data("death.2019")
```

Then apply the `sCAwithPMD` to the data:

``` r
sca.res <- sCAwithPMD(
  DATA = death.2019, # Contingency table
  dimensions = 2L, # the number of dimensions
  doublecentering = TRUE, # center the data
  s1 = rep(0.5 * sqrt(nrow(death.2019)), 2), # Asking for a medium amount of sparsity
  s2 = rep(0.5 * sqrt(ncol(death.2019)), 2)
)
```

``` r
sca.fi.map.12 <- createFactorMap(X = sca.res$fi,
                       col.background = NULL,
                       col.axes = "#42376B", 
                       width.axes = 0.5,
                       title = "SCA: row factor scores",
                       alpha.axes = 0.5,
                       alpha.points = 0.5,
                       pch = 16,
                       axis1 = 1,
                       axis2 = 2,
                       constraints = NULL, text.cex = 4)

sca.fi.plot.12 <- sca.fi.map.12$zeMap_background + sca.fi.map.12$zeMap_dots + sca.fi.map.12$zeMap_text + geom_path(color = "darkorchid4") + theme(axis.title = element_text(color = "#42376B"), axis.text = element_text(color = "#42376B"), title = element_text(color = "#42376B"), panel.border = element_rect(size = 1.5, color = "#42376B", fill = NA)) + labs(x = "Dimension 1", y = "Dimension 2")

sca.fi.plot.12
```

<img src="man/figures/README-row factor map-1.svg" width="100%" />

``` r
sca.fj.map.12 <- createFactorMap(X = sca.res$fj,
                       col.background = NULL,
                       col.axes = "#42376B", 
                       width.axes = 0.5,
                       title = "SCA: row factor scores",
                       alpha.axes = 0.5,
                       alpha.points = 0.5,
                       pch = 16,
                       axis1 = 1,
                       axis2 = 2,
                       constraints = NULL, 
                       text.cex = 4)

sca.fj.plot.12 <- sca.fj.map.12$zeMap_background + sca.fj.map.12$zeMap_dots + sca.fj.map.12$zeMap_text + theme(axis.title = element_text(color = "#42376B"), axis.text = element_text(color = "#42376B"), title = element_text(color = "#42376B"), panel.border = element_rect(size = 1.5, color = "#42376B", fill = NA)) + labs(x = "Dimension 1", y = "Dimension 2")

sca.fj.plot.12
#> Warning: ggrepel: 8 unlabeled data points (too many overlaps). Consider
#> increasing max.overlaps
```

<img src="man/figures/README-column factor map-1.svg" width="100%" />
