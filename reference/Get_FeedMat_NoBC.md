# Estimate feedback matrix from several models (No bias-corrected version)

Estimate feedback matrix from several models (No bias-corrected version)

## Usage

``` r
Get_FeedMat_NoBC(
  RiskFactors,
  ModelType,
  Economies,
  GVARinputs,
  JLLinputs,
  FactorLabels,
  CheckInpts = FALSE
)
```

## Arguments

- RiskFactors:

  A numeric matrix (T x F) representing the time series of risk factors.

- Economies:

  A character vector containing the names of the economies included in
  the system.

- GVARinputs:

  List. Inputs for GVAR model estimation.

- JLLinputs:

  List. Inputs for JLL model estimation.

- FactorLabels:

  A list of character vectors with labels for all variables in the
  model.

- CheckInpts:

  Check for the consistency of the provided inputs. Default is FALSE.
