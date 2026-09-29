
<!-- README.md is generated from README.Rmd. Please edit that file -->

# NMFregress

## Goal:

Convert text into a term document matrix and decompose it using
separable NMF so that you can assess effect sizes of covariates on the
topic allocation.

## Installation

You can install the development version of NMFregress from
[GitHub](https://github.com/) with:

``` r
devtools::install_github("iamdavecampbell/NMFregress", build_vignettes = TRUE)
```

## Example

- List the vignette name:

``` r
library(NMFregress)
vignette(package="NMFregress")
```

- See the vignette: **romeo_and_juliet**

``` r
vignette("romeo_juliet", package = "NMFregress")
#> starting httpd help server ... done
```

## Published article:

G.Phelan and D. A.Campbell, “Testing Hypotheses of Covariate Effects on
Topics of Discourse,” Statistical Analysis and Data Mining: An ASA Data
Science Journal 19, no. 2 (2026): e70066,
<https://doi.org/10.1002/sam.70066>.
\[<https://onlinelibrary.wiley.com/share/DWKBHB6HB58AKPJWIYM4?target=10.1002/sam.70066>\]
