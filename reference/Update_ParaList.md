# Update parameters in the optimization process

Update parameters in the optimization process

## Usage

``` r
Update_ParaList(
  x,
  ModelType,
  FactorLabels,
  Economies,
  JLLinputs = NULL,
  GVARinputs = NULL,
  ListInputSet
)
```

## Arguments

- x:

  vector containing all the vectorized auxiliary parameters

- ModelType:

  string-vector containing the label of the model to be estimated

- FactorLabels:

  string-list based which contains the labels of all the variables
  present in the model

- Economies:

  string-vector containing the names of the economies which are part of
  the economic system

- JLLinputs:

  Set of necessary inputs used in the estimation of the JLL-based models

- GVARinputs:

  Set of necessary inputs used in the estimation of the GVAR-based
  models

- ListInputSet:

  variable inputs used in the optimization (see "Optimization" function)
