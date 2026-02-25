# Compute FEVDs for all models

Compute FEVDs for all models

## Usage

``` r
ComputeFEVDs(
  SIGMA,
  K1Z,
  G0,
  BLoad,
  FactorLabels,
  FacDim,
  MatLength,
  FEVDhoriz,
  YieldsLabel,
  ModelType,
  Economy = NULL,
  CholFac_JLL = NULL,
  PI = NULL,
  Mode = FALSE
)
```

## Arguments

- SIGMA:

  Variance-covariance matrix

- K1Z:

  Loading As

- G0:

  contemporaneous terms

- BLoad:

  Loading Bs

- FactorLabels:

  List containing the label of factors

- FacDim:

  Dimension of the P-dynamics

- MatLength:

  Length of the maturity vector

- FEVDhoriz:

  Horizon of the analysis

- YieldsLabel:

  Label of bond yields

- ModelType:

  Desired model type

- Economy:

  specific economy under study

- CholFac_JLL:

  Cholesky factorization term from JLL models

- PI:

  matrix PI for JLL-based models

- Mode:

  allows for the orthogonalized version in the case of JLL-based models

## References

- This function is a modified and extended version of the "fevd"
  function from Smith, L.V. and A. Galesi (2014). GVAR Toolbox 2.0,
  available at https://sites.google.com/site/gvarmodelling/gvar-toolbox.

- Pesaran and Shin, 1998. "Generalized impulse response analysis in
  linear multivariate models" (Economics Letters)
