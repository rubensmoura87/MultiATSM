# IRFs and GIRFs for all models

IRFs and GIRFs for all models

## Usage

``` r
IRFandGIRF(ModelType, ModelPara, IRFhoriz, FactorLabels, Economies)
```

## Arguments

- ModelType:

  string-vector containing the label of the model to be estimated

- ModelPara:

  list of model parameter estimates (See the "Optimization" function)

- IRFhoriz:

  single numerical vector containing the desired horizon of analysis for
  the IRFs

- FactorLabels:

  string-list based which contains the labels of all the variables
  present in the model

- Economies:

  string-vector containing the names of the economies which are part of
  the economic system

## Details

The Structural shocks from the IRFs are identified via Cholesky
decomposition
