# Gathers the general inputs for model estimation

Gathers the general inputs for model estimation

## Usage

``` r
GeneralMLEInputs(
  Yields,
  RiskFactors,
  FactorLabels,
  mat,
  DataFrequency,
  Label_Multi_Models,
  Economies,
  ModelType,
  UnitYields
)
```

## Arguments

- Yields:

  matrix (CJ x T) or a list containing those matrices, where C is the
  number of countries, J - the number of maturities and T - time series
  length  

- RiskFactors:

  time series of risk factors (F x T). Could be stored in a list
  depending on the model

- FactorLabels:

  string-list based which contains the labels of all variables present
  in the model

- mat:

  vector of maturities (in years) used in the estimation

- DataFrequency:

  single element character-based vector. Available options are: "Daily
  All Days",  
  "Daily Business Days", "Weekly", "Monthly", "Quarterly", "Annually"

- Label_Multi_Models:

  string-vector containing the names of the multicountry setups

- Economies:

  string-vector containing the names of the economies which are part of
  the economic system

- ModelType:

  string-vector containing the label of the model to be estimated
