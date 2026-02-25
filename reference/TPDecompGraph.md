# Term Premia decomposition graphs for all models

Term Premia decomposition graphs for all models

## Usage

``` r
TPDecompGraph(
  ModelType,
  NumOut,
  ModelPara,
  WishRPgraphs,
  UnitYields,
  Economies,
  PathsGraphs,
  Folder2Save,
  verbose
)
```

## Arguments

- ModelType:

  character. Estimated model type. Permissible choices: "JPS original",
  "JPS global", "GVAR single", "JPS multi", "GVAR multi", "JLL
  original", "JLL No DomUnit", "JLL joint Sigma".

- NumOut:

  list. Computed outputs containing model fit, IRFs, FEVDs, GIRFs,
  GFEVDs and risk premia.

- ModelPara:

  list. Model parameter estimates (see
  [`Optimization`](https://rubensmoura87.github.io/MultiATSM/reference/Optimization.md)).

- WishRPgraphs:

  logical. Set TRUE to generate term premia graphs, FALSE otherwise.

- UnitYields:

  character. "Month" if yields are in months, "Year" if in years.

- Economies:

  character vector. Names of the `C` economies included in the system.

- PathsGraphs:

  character. Path of the folder in which the graphs will be saved.

- Folder2Save:

  character. Folder path where the outputs will be stored.

- verbose:

  logical. Flag controlling function messaging.

## Available Methods

\- \`autoplot(object, type = "TermPremia")\`

## Examples

``` r
data("ParaSetEx")
data("NumOutEx")
ModelType <- "JPS original"
Economy <- "Brazil"
UnitYields <- "Month"
TPDecompGraph(ModelType, NumOutEx, ParaSetEx,
  WishRPgraphs = FALSE, UnitYields, Economy,
  PathsGraphs = NULL, Folder2Save = NULL, verbose = FALSE
)
```
