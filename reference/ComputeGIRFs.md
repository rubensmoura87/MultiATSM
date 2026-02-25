# Compute GIRFs for all models

Compute GIRFs for all models

## Usage

``` r
ComputeGIRFs(
  Sigma.y,
  F1,
  BLoad,
  G0.y,
  FactorLabels,
  FacDim,
  MatLength,
  GIRFhoriz,
  YieldsLabel,
  ModelType,
  Economy = NULL,
  PI = NULL,
  Mode = FALSE
)
```

## Arguments

- Sigma.y:

  Variance-covariance matrix

- F1:

  Feedback matrix

- BLoad:

  Loading Bs

- G0.y:

  matrix of contemporaneous terms

- FactorLabels:

  List containing the labels of the factors

- FacDim:

  Dimension of the P-dynamics

- MatLength:

  Length of the maturity vector

- GIRFhoriz:

  Horizon of the analysis

- YieldsLabel:

  Label o yields

- ModelType:

  desired Model type

- Economy:

  Economy under study

- PI:

  matrix PI for JLL-based models

- Mode:

  allows for the orthogonalized version in the case of JLL-based models

  \#' @references

  - This function is partially based on the version of the "irf"
    function from Smith, L.V. and A. Galesi (2014). GVAR Toolbox 2.0,
    available at
    https://sites.google.com/site/gvarmodelling/gvar-toolbox.

  - Pesaran and Shin, 1998. "Generalized impulse response analysis in
    linear multivariate models" (Economics Letters)
