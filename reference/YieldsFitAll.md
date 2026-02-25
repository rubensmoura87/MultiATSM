# Fit yields for all maturities of interest

Fit yields for all maturities of interest

## Usage

``` r
YieldsFitAll(MatInt, ModelPara, FactorLabels, ModelType, Economies, YLab)
```

## Arguments

- MatInt:

  numerical vector containing the fit maturities of interest

- ModelPara:

  List of model parameter estimates (See the "Optimization" function)

- FactorLabels:

  a string-list based which contains all the labels of all the variables
  present in the model

- ModelType:

  a string-vector containing the label of the model to be estimated

- Economies:

  a string-vector containing the names of the economies which are part
  of the economic system

- YLab:

  Label of yields ("Months" or "Yields")
