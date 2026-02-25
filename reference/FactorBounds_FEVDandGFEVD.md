# Compute the confidence bounds for the model bond P-dynamics-related outputs

Compute the confidence bounds for the model bond P-dynamics-related
outputs

## Usage

``` r
FactorBounds_FEVDandGFEVD(
  ModelBootstrap,
  quants,
  ModelType,
  ndraws,
  Horiz,
  FacDim,
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

- LabFEVD:

  vector containing the labels "FEVD" and "GFEVD"

- Economies:

  Economies that are part of the economic system
