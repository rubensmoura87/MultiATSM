# FEVDs and GFEVDs for all models

FEVDs and GFEVDs for all models

## Usage

``` r
FEVDandGFEVD(ModelType, ModelPara, FEVDhoriz, FactorLabels, Economies)
```

## Arguments

- ModelType:

  string-vector containing the label of the model to be estimated

- ModelPara:

  list of model parameter estimates (see the "Optimization" function)

- FEVDhoriz:

  single numerical vector containing the desired horizon of analysis for
  the FEVDs and GFEVDs

- FactorLabels:

  string-list based which contains all the labels of all the variables
  present in the model

- Economies:

  string-vector containing the names of the economies which are part of
  the economic system

## Details

Structural shocks are identified via Cholesky decomposition
