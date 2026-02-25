# Compute the B matrix of loadings

Compute the B matrix of loadings

## Usage

``` r
Get_BFull(ModelParaPE, FactorLabels, mat, Economies, ModelType)
```

## Arguments

- ModelParaPE:

  List of point estimates of the model parameter

- FactorLabels:

  String-list based which contains the labels of all the variables
  present in the model

- mat:

  Vector of bond yield maturities

- Economies:

  String-vector containing the names of the economies which are part of
  the economic system

- ModelType:

  A character vector indicating the model type to be estimated
