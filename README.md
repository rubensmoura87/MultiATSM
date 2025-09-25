
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/MultiATSM)](https://cran.r-project.org/package=MultiATSM)
[![R-CMD-check](https://github.com/mothur/phylotypr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mothur/phylotypr/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/mothur/phylotypr/branch/main/graph/badge.svg)](https://app.codecov.io/gh/mothur/phylotypr?branch=main)
[![CRAN
downloads](https://cranlogs.r-pkg.org/badges/MultiATSM)](https://CRAN.R-project.org/package=MultiATSM)

<!-- badges: end -->

## Overview

MultiATSM is an R package for estimating, analyzing, and forecasting
multi-country macro-finance affine term structure models (ATSMs). The
package implements a variety of models, including JPS, GVAR, and JLL
frameworks (see complete references below), and provides tools for bias
correction, bootstrap analysis, and graphical/numerical outputs.

## Features

- Estimation of multi-country ATSMs (JPS, 2014, JLL, 2015, GVAR,
  2023-2024, and variants)
- Flexible input structure for macroeconomic and yield data
- Graphical and numerical outputs for model diagnostics and
  interpretation
- Bias correction for a VAR dynamics (BRW, 2012)
- Bootstrap analysis for confidence intervals
- Out-of-sample forecasting

## Installation

You can install the development version of MultiATSM from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("rubensmoura87/MultiATSM")
```

You can also get the official release version from CRAN

``` r
install.packages("MultiATSM")
```

## Usage and documentation

- See the package manual and
  [vignettes](https://https://rubensmoura87.github.io/articles/MultiATSM/)
  for detailed function documentation and examples.
- Main functions: `InputsForOpt`, `Bias_Correc_VAR`, `Optimization`,
  `NumOutputs`, `Bootstrap`, `ForecastYields`, etc.

## More information about `{MultiATSM}`

You can learn more about the underlying models in the working paper
published in the [LIDAM-LFIN
Series.](https://dial.uclouvain.be/pr/boreal/object/boreal%3A259119/datastream/PDF_01/view#:~:text=Abstract,sample%20forecasting%20of%20bond%20yields.).

## References

- Bauer, M. D., Rudebusch, G. D., & Wu, J. C. (2012). Correcting
  estimation bias in dynamic term structure models. *Journal of Business
  & Economic Statistics*, 30(3), 454-467.

- Candelon, B., & Moura, R. (2023). Sovereign yield curves and the
  COVID-19 in emerging markets. *Economic Modelling*, 127, 106453.

- Candelon, B., & Moura, R. (2024). A Multicountry Model of the Term
  Structures of Interest Rates with a GVAR. *Journal of Financial
  Econometrics*, 22(5), 1558-1587.

- Joslin, S., Priebsch, M., & Singleton, K. J. (2014). Risk premiums in
  dynamic term structure models with unspanned macro risks. *The Journal
  of Finance*, 69(3), 1197-1233.

- Joslin, S., Singleton, K. J., & Zhu, H. (2011). A new perspective on
  Gaussian dynamic term structure models. *The Review of Financial
  Studies*, 24(3), 926-970.

- Jotikasthira, C., Le, A., & Lundblad, C. (2015). Why do term
  structures in different currencies co-move?. *Journal of Financial
  Economics*, 115(1), 58-83.

- See also the references in the package documentation.
