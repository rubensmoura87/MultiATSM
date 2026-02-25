# Creates the confidence bounds and the graphs of IRFs and GIRFs after bootstrap

Creates the confidence bounds and the graphs of IRFs and GIRFs after
bootstrap

## Usage

``` r
IRFandGIRFbs(
  ModelType,
  ModelBootstrap,
  NumOutPE,
  InputsForOutputs,
  Economies,
  PathsGraphs,
  Folder2save,
  verbose
)
```

## Arguments

- ModelType:

  string-vector containing the label of the model to be estimated

- ModelBootstrap:

  list containing the complete set of model parameters after bootstrap
  estimation procedure

- NumOutPE:

  list of model parameter point estimates

- InputsForOutputs:

  list containing the desired inputs for the construction of the outputs
  of interest

- Economies:

  string-vector containing the names of the economies which are part of
  the economic system

- PathsGraphs:

  path of the folder in which the graphs will be saved

- Folder2save:

  Folder path where the outputs will be stored.

- verbose:

  Logical flag controlling function messaging.
