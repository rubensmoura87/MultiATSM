# FEVDs and GFEVDs after bootstrap for all models

FEVDs and GFEVDs after bootstrap for all models

## Usage

``` r
FEVDandGFEVD_BS(ModelType, ModelParaBoot, FEVDhoriz, FactorLabels, Economies)
```

## Arguments

- ModelType:

  string-vector containing the label of the model to be estimated

- ModelParaBoot:

  list of model parameter estimates (see the "Optimization" function)
  after a bootstrap draw

- FEVDhoriz:

  single numerical vector containing the desired horizon of analysis for
  the FEVDs

- FactorLabels:

  string-list based which contains all the labels of all the variables
  present in the model

- Economies:

  string-vector containing the names of the economies which are part of
  the economic system
