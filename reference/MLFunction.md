# Set up the vector-valued objective function (Point estimate)

Set up the vector-valued objective function (Point estimate)

## Usage

``` r
MLFunction(
  MLEinputs,
  Economies,
  DataFrequency,
  FactorLabels,
  ModelType,
  BS_outType = FALSE
)
```

## Arguments

- MLEinputs:

  Set of inputs that are necessary to the log-likelihood function

- Economies:

  string-vector containing the names of the economies which are part of
  the economic system

- DataFrequency:

  character-based vector: "Daily All Days", "Daily Business Days",
  "Weekly", "Monthly", "Quarterly", "Annually"

- FactorLabels:

  string-list based which contains the labels of all the variables
  present in the model

- ModelType:

  string-vector containing the label of the model to be estimated

- BS_outType:

  Generates simplified output list in the bootstrap setting. Default is
  set to FALSE.

## Value

objective function
