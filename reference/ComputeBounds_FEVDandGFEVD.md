# Compute the confidence bounds around the P-dynamics and bond yields for FEVD and GFEVD

Compute the confidence bounds around the P-dynamics and bond yields for
FEVD and GFEVD

## Usage

``` r
ComputeBounds_FEVDandGFEVD(
  ModelBootstrap,
  quants,
  FacDim,
  YieDim,
  ModelType,
  Economies,
  ndraws,
  Horiz
)
```

## Arguments

- ModelBootstrap:

  numerical output set from the bootstrap analysis

- quants:

  quantile of the confidence bounds

- FacDim:

  dimension of the risk factor set

- YieDim:

  dimension of the bond yield set

- ModelType:

  desired model type

- Economies:

  Economies that are part of the economic system

- ndraws:

  number of draws

- Horiz:

  horizon of numerical outputs
