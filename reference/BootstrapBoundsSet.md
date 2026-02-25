# Builds the confidence bounds and graphs (Bootstrap set)

Builds the confidence bounds and graphs (Bootstrap set)

## Usage

``` r
BootstrapBoundsSet(
  ModelType,
  ModelBootstrap,
  NumOutPE,
  InputsForOutputs,
  Economies,
  Folder2save,
  verbose
)
```

## Arguments

- ModelType:

  string-vector containing the label of the model to be estimated

- ModelBootstrap:

  list containing the complete set of model parameters after the
  bootstrap estimation procedure

- NumOutPE:

  point estimate from the numerical outputs (see the outputs of the
  "NumOutputs" function)

- InputsForOutputs:

  list containing the desired inputs for the construction of IRFs,
  GIRFs, FEVDs, and GFEVDs

- Economies:

  string-vector containing the names of the economies which are part of
  the economic system

- Folder2save:

  Folder path where the outputs will be stored.

- verbose:

  Logical flag controlling function messaging.
