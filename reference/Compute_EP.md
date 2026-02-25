# Compute the expected component for all models

Compute the expected component for all models

## Usage

``` r
Compute_EP(
  ModelPara,
  ModelType,
  UnitYields,
  matAdjUnit,
  N,
  rhoList,
  Economies,
  FactorLabels,
  WishFP = FALSE,
  matMIN = FALSE,
  matMAX = FALSE
)
```

## Arguments

- ModelPara:

  list of model parameter estimates

- ModelType:

  Desired model type

- UnitYields:

  Available options: "Month" and "Year"

- matAdjUnit:

  Adjusted vector of matutities

- N:

  number of country-specific spanned factors

- rhoList:

  List of risk-neutral parameters

- Economies:

  string-vector containing the names of the economies which are part of
  the economic system

- FactorLabels:

  List of factor labels

- WishFP:

  If users wants to compute the forward premia. Default is FALSE.

- matMIN:

  For the forward premia, the shortest maturity of the remium of
  interest

- matMAX:

  For the forward premia, the longest maturity of the remium of interest
