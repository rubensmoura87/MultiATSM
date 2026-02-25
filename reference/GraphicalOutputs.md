# Generate the graphical outputs for the selected models (Point estimate)

Generate the graphical outputs for the selected models (Point estimate)

## Usage

``` r
GraphicalOutputs(
  ModelType,
  ModelPara,
  NumOut,
  InputsForOutputs,
  Economies,
  FactorLabels,
  Folder2save,
  verbose
)
```

## Arguments

- ModelType:

  A character vector indicating the model type to be estimated.

- ModelPara:

  List of model parameter estimates (See the
  [`Optimization`](https://rubensmoura87.github.io/MultiATSM/reference/Optimization.md)
  function)

- NumOut:

  list of computed outputs containing the model fit, IRFs, FEVDs, GIRFs,
  GFEVDs and Term Premia

- InputsForOutputs:

  list containing the desired inputs for the construction of the desired
  output

- Economies:

  A character vector containing the names of the economies included in
  the system.

- FactorLabels:

  A list of character vectors with labels for all variables in the
  model.

- Folder2save:

  Folder path where the outputs will be stored.

- verbose:

  Logical flag controlling function messaging.
