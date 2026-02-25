# Use function ML to generate the outputs from a ATSM

Use function ML to generate the outputs from a ATSM

## Usage

``` r
FunctionML_vec(
  x,
  ML_fun,
  ListInputSet,
  ModelType,
  FactorLabels,
  Economies,
  JLLinputs,
  GVARinputs,
  WithEstimation
)
```

## Arguments

- x:

  vector containing all the vectorized auxiliary parameters

- ML_fun:

  vector-valued objective function (function)

- ListInputSet:

  variable inputs used in the optimization (see inputs from
  "optimization" function)

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
  (see "JLL" function)

- GVARinputs:

  Set of necessary inputs used in the estimation of the GVAR-based
  models (see "GVAR" function)

- WithEstimation:

  if TRUE, returns only the values of the likelihood function, else
  generates the entire set of outputs
