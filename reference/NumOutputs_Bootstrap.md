# Numerical outputs (IRFs, GIRFs, FEVD, and GFEVD) for bootstrap

Numerical outputs (IRFs, GIRFs, FEVD, and GFEVD) for bootstrap

## Usage

``` r
NumOutputs_Bootstrap(
  ModelType,
  ModelParaBoot,
  InputsForOutputs,
  FactorLabels,
  Economies
)
```

## Arguments

- ModelType:

  A character vector indicating the model type to be estimated.

- ModelParaBoot:

  A list of model parameter estimates (see the "Optimization" function)
  after a bootstrap draw

- InputsForOutputs:

  A list containing the necessary inputs for generating IRFs, GIRFs,
  FEVDs, GFEVDs and Term Premia.

- FactorLabels:

  A list of character vectors with labels for all variables in the
  model.

- Economies:

  A character vector containing the names of the economies included in
  the system.
