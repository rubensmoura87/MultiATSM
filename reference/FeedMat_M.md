# Computes an average or median feedback matrix across several bootstrap iterations

Computes an average or median feedback matrix across several bootstrap
iterations

## Usage

``` r
FeedMat_M(
  K1Z_j,
  N_boot,
  RiskFactors,
  GVARinputs,
  JLLinputs,
  FactorLabels,
  Economies,
  ModelType,
  Use_Mean
)
```

## Arguments

- K1Z_j:

  Feedback matrix at the j_th iteration

- N_boot:

  Number of bootstrap samples per iteration.

- RiskFactors:

  A numeric matrix (T x F) representing the time series of risk factors.

- GVARinputs:

  List. Inputs for GVAR model estimation.

- JLLinputs:

  List. Inputs for JLL model estimation.

- FactorLabels:

  A list of character vectors with labels for all variables in the model

- Economies:

  A character vector containing the names of the economies included in
  the system.

- ModelType:

  A character vector indicating the model type to be estimated.

- Use_Mean:

  Choose between the mean or the median across the bootstrap iterations.
