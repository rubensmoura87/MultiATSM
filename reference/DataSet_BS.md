# Prepare the factor set for GVAR models (Bootstrap version)

Prepare the factor set for GVAR models (Bootstrap version)

## Usage

``` r
DataSet_BS(ModelType, RiskFactors, Wgvar, Economies, FactorLabels)
```

## Arguments

- ModelType:

  A character vector containing the label of the model to be estimated.

- RiskFactors:

  A matrix of the complete set of risk factors (K x T).

- Wgvar:

  A transition matrix from GVAR models (C x C).

- Economies:

  A character vector containing the names of the economies included in
  the system.

- FactorLabels:

  A list of character vectors with labels for all variables in the
  model.

## Value

A list containing the factor set for GVAR models.
