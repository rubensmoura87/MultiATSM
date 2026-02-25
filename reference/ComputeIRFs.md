# Compute IRFs of all models

Compute IRFs of all models

## Usage

``` r
ComputeIRFs(
  SIGMA,
  K1Z,
  BLoad,
  FactorLabels,
  FacDim,
  MatLength,
  IRFhoriz,
  YieldsLabel,
  ModelType,
  Economy = NULL,
  PI = NULL,
  Mode = FALSE
)
```

## Arguments

- SIGMA:

  Variance-covariance matrix

- K1Z:

  Loading As

- BLoad:

  Loading Bs

- FactorLabels:

  List containing the label of factors

- FacDim:

  Dimension of the P-dynamics

- MatLength:

  Length of the maturity vector

- IRFhoriz:

  Horizon of the analysis

- YieldsLabel:

  Label of bond yields

- ModelType:

  Desired model type

- Economy:

  specific economy under study

- PI:

  matrix PI for JLL-based models

- Mode:

  allows for the orthogonalized version in the case of JLL-based models
