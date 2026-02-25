# Model fit graphs for all models

Model fit graphs for all models

## Usage

``` r
Fitgraphs(
  ModelType,
  WishFitgraphs,
  ModelPara,
  NumOut,
  Economies,
  PathsGraphs,
  Folder2save,
  verbose
)
```

## Arguments

- ModelType:

  character. Estimated model type.

- WishFitgraphs:

  logical. Set TRUE to generate fit graphs, FALSE otherwise.

- ModelPara:

  list. Model parameter estimates (see
  [`Optimization`](https://rubensmoura87.github.io/MultiATSM/reference/Optimization.md)).

- NumOut:

  list. Outputs containing model fit, IRFs, FEVDs, GIRFs, GFEVDs and
  Term premia.

- Economies:

  character vector. Names of the economies included in the system.

- PathsGraphs:

  character. Path of the folder in which the graphs will be saved.

- Folder2save:

  character. Desired folder path to save outputs.

- verbose:

  logical. Flag controlling function messaging.

## Available Methods

\- \`autoplot(object, type = "Fit")\`

## Examples

``` r
data("ParaSetEx")
data("NumOutEx")
ModelType <- "JPS original"
Economy <- "Brazil"
Fitgraphs(ModelType,
  WishFitgraphs = TRUE, ParaSetEx, NumOutEx, Economy, PathsGraphs = NULL,
  Folder2save = NULL, verbose = FALSE
)
#> $Brazil

#> 
```
