# Compute quantiles for model P-dynamics

Compute quantiles for model P-dynamics

## Usage

``` r
FacQuantile_bs(
  DrawSet,
  LabIRF,
  ndraws,
  quants,
  Horiz,
  FacDim,
  DimLabelsFac,
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

- DimLabelsFac:

  labels of the factor set

- ModelType:

  desired model type

- Ortho:

  Orthogonolized version for the JLL models. Default is FALSE.
