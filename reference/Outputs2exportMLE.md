# Prepares inputs to export

Prepares inputs to export

## Usage

``` r
Outputs2exportMLE(
  Label_Multi_Models,
  Economies,
  RiskFactors,
  Yields,
  mat,
  ModelInputsGen,
  ModelInputsSpec,
  PdynPara,
  ModelType
)
```

## Arguments

- Label_Multi_Models:

  string-vector containing the names of the multicountry setups

- Economies:

  string-vector containing the names of the economies which are part of
  the economic system

- RiskFactors:

  time series of risk factors (F x T). Could be stored in a list
  depending on the model

- Yields:

  matrix (CJ x T) or a list containing those matrices, where C is the
  number of countries, J - the number of maturities and T - time series
  length  

- mat:

  vector of maturities (in years) used in the estimation

- ModelInputsGen:

  List of generic inputs

- ModelInputsSpec:

  List of specific inputs

- PdynPara:

  Model parameters estimated in the P-dynamics the

- ModelType:

  string-vector containing the label of the model to be estimated
