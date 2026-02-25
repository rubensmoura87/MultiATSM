# Compute the confidence bounds for the model bond yield-related outputs

Compute the confidence bounds for the model bond yield-related outputs

## Usage

``` r
YieldBounds_FEVDandGFEVD(
  ModelBootstrap,
  quants,
  ModelType,
  ndraws,
  Horiz,
  FacDim,
  YieDim,
  LabFEVD,
  Economies
)
```

## Arguments

- ModelBootstrap:

  numerical output set from the bootstrap analysis

- quants:

  quantile of the confidence bounds

- ModelType:

  desired model type

- ndraws:

  number of draws

- Horiz:

  horizon of numerical outputs

- FacDim:

  dimension of the risk factor set

- YieDim:

  dimension of the bond yield set

- LabFEVD:

  vector containing the labels "FEVD" and "GFEVD"

- Economies:

  Economies that are part of the economic system
