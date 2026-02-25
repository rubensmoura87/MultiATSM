# IRFs and GIRFs after bootstrap for all models

IRFs and GIRFs after bootstrap for all models

## Usage

``` r
IRFandGIRF_BS(ModelType, ModelParaBoot, IRFhoriz, FactorLabels, Economies)
```

## Arguments

- ModelType:

  string-vector containing the label of the model to be estimated

- ModelParaBoot:

  list of model parameter estimates (see the "Optimization" function)
  after a bootstrap draw

- IRFhoriz:

  single numerical vector containing the desired horizon of analysis for
  the IRFs

- FactorLabels:

  string-list based which contains all the labels of all the variables
  present in the model

- Economies:

  string-vector containing the names of the economies which are part of
  the economic system
