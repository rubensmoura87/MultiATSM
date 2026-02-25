# Estimates an unbiased VAR(1) using stochastic approximation (Bauer, Rudebusch and Wu, 2012)

Estimates an unbiased VAR(1) using stochastic approximation (Bauer,
Rudebusch and Wu, 2012)

## Usage

``` r
Bias_Correc_VAR(
  ModelType,
  BRWinputs,
  RiskFactors,
  Economies,
  FactorLabels,
  GVARinputs = NULL,
  JLLinputs = NULL,
  verbose = TRUE
)
```

## Arguments

- ModelType:

  character. Model type to be estimated. Permissible choices: "JPS
  original", "JPS global", "GVAR single", "JPS multi", "GVAR multi",
  "JLL original", "JLL No DomUnit", "JLL joint Sigma".

- BRWinputs:

  list. Contains the necessary inputs for the BRW model estimation:

  1.  `Cent_Measure`: "Mean" or "Median" (unbiased estimation type)

  2.  `gamma`: Numeric. Adjustment parameter between 0 and 1. Default is
      0.5.

  3.  `N_iter`: Integer. Number of iterations for the stochastic
      approximation algorithm after burn-in. Default is 5000.

  4.  `N_burn`: Integer. Number of burn-in iterations. Default is 15

  5.  `B`: Integer. Number of bootstrap samples per iteration for
      calculating the noisy measure of the biased estimator's mean or
      median. Default is 50.

  6.  `check`: Logical. Indicates whether to perform a closeness check.
      Default is TRUE.

  7.  `B_check`: Integer. Number of bootstrap samples for the closeness
      check. Default is 100000.

  8.  `Eigen_rest`: Numeric. Restriction on the largest eigenvalue under
      the P-measure. Default is 1.

- RiskFactors:

  numeric matrix (`Td x K`). Time series of risk factors.

- Economies:

  character vector. Names of the `C` economies included in the system.

- FactorLabels:

  list. Labels for all variables in the model.

- GVARinputs:

  list. Inputs for GVAR model estimation (see
  [`GVAR`](https://rubensmoura87.github.io/MultiATSM/reference/GVAR.md)).
  Default is NULL.

- JLLinputs:

  list. Inputs for JLL model estimation (see
  [`JLL`](https://rubensmoura87.github.io/MultiATSM/reference/JLL.md)).
  Default is NULL.

- verbose:

  logical. Flag controlling function messaging. Default TRUE.

## Value

Bias-corrected VAR parameters based on the framework of Bauer, Rudebusch
and Wu (2012). The list contains:

1.  `KOZ_BC`: estimated intercept (K x 1);

2.  `K1Z_BC`: estimated feedback matrix (K x K);

3.  `SSZ_BC`: estimated variance-covariance matrix (K x K);

4.  `dist`: root mean square distance (scalar);

## General Notation

- `Td` denotes the model time series dimension.

- `C` number of countries in the system.

- `K` denotes the total number of risk factors.

## References

Bauer, Rudebusch and, Wu (2012). "Correcting Estimation Bias in Dynamic
Term Structure Models"  
This function offers an independent R implementation that is informed by
the conceptual framework outlined in Bauer, Rudebusch and Wu (2012), but
adapted to the present modeling context. Related Matlab routines are
available on Cynthia Wu's website
(https://sites.google.com/view/jingcynthiawu/).

## Examples

``` r
# \donttest{
data(RiskFacFull)
Factors <- t(RiskFacFull[1:7, ])

BRWinputs <- list(
  Cent_Measure = "Mean", gamma = 0.4, N_iter = 1000, N_burn = 100,
  B = 10, check = 1, B_check = 5000
)

Economies <- "China"
N <- 3
ModelType <- "JPS original"
FactorLabels <- NULL

BRWpara <- Bias_Correc_VAR(ModelType, BRWinputs, Factors, Economies, FactorLabels, verbose = FALSE)
# }
```
