# Compute the confidence bounds for the model P-dynamics

Compute the confidence bounds for the model P-dynamics

## Usage

``` r
FactorBounds_IRFandGIRF(
  ModelBootstrap,
  quants,
  ModelType,
  ndraws,
  Horiz,
  FacDim,
  LabIRF,
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

- LabIRF:

  vector containing the labels "IRF" and "GIRF"

- Economies:

  Economies that are part of the economic system
