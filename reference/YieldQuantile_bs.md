# Compute quantiles for model bond yield-related outputs

Compute quantiles for model bond yield-related outputs

## Usage

``` r
YieldQuantile_bs(
  DrawSet,
  LabIRF,
  ndraws,
  quants,
  Horiz,
  FacDim,
  YieDim,
  LabelsYies,
  ModelType,
  Ortho = FALSE
)
```

## Arguments

- DrawSet:

  Draw-specific set

- LabIRF:

  vector containing the labels "IRF" and "GIRF"

- ndraws:

  number of draws

- quants:

  quantile of the confidence bounds

- Horiz:

  horizon of numerical outputs

- FacDim:

  dimension of the risk factor set

- YieDim:

  dimension of the bond yield set

- LabelsYies:

  labels of the factor set

- ModelType:

  desired model type

- Ortho:

  Orthogonolized version for the JLL models. Default is FALSE.
