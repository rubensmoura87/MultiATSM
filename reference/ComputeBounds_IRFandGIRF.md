# Compute the confidence bounds from the model's numerical outputs

Compute the confidence bounds from the model's numerical outputs

## Usage

``` r
ComputeBounds_IRFandGIRF(
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

  Desired model type

- Economies:

  Economies that are part of the economic system

- ndraws:

  number of draws selected

- Horiz:

  horizon of numerical outputs
