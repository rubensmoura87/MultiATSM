# IRF and GIRF graphs for all models

IRF and GIRF graphs for all models

## Usage

``` r
IRFandGIRFgraphs(
  ModelType,
  NumOut,
  WishPdynamicsgraphs,
  WishYieldsgraphs,
  IRFhoriz,
  PathsGraphs,
  OutputType,
  Economies,
  Folder2save,
  verbose
)
```

## Arguments

- ModelType:

  character. Estimated model type.Permissible choices: "JPS original",
  "JPS global", "GVAR single", "JPS multi", "GVAR multi", "JLL
  original", "JLL No DomUnit", "JLL joint Sigma".

- NumOut:

  list. Computed outputs containing model fit, IRFs, FEVDs, GIRFs,
  GFEVDs and term premia.

- WishPdynamicsgraphs:

  logical. Set TRUE to generate risk factor graphs, FALSE otherwise.

- WishYieldsgraphs:

  logical. Set TRUE to generate bond yield graphs, FALSE otherwise.

- IRFhoriz:

  integer. Desired horizon of analysis for the IRFs.

- PathsGraphs:

  character. Path of the folder in which the graphs will be saved.

- OutputType:

  character. Available options: "IRF", "GIRF", "IRF Ortho", "GIRF
  Ortho".

- Economies:

  character vector. Names of the `C` economies included in the system.

- Folder2save:

  character. Folder path where the outputs will be stored.

- verbose:

  logical. Flag controlling function messaging.

## Available Methods

\- \`autoplot(object, type = "IRF_Factor")\`, \`autoplot(object, type =
"IRF_Yields")\`, \`autoplot(object, type = "GIRF_Yields")\`,
\`autoplot(object, type = "GIRF_Yields")\`. For JLL-based models:
\`autoplot(object, type = "IRF_Factor-\_Ortho")\`,  
\`autoplot(object, type = "IRF_Yields_Ortho")\`, \`autoplot(object, type
= "GIRF_Yields_Ortho")\`, \`autoplot(object, type =
"GIRF_Yields_Ortho")\`.

## Examples

``` r
data("NumOutEx")
ModelType <- "JPS original"
Economy <- "Brazil"
IRFhoriz <- 20
irf_Out <- IRFandGIRFgraphs(ModelType, NumOutEx,
  WishPdynamicsgraphs = FALSE, WishYieldsgraphs = TRUE, IRFhoriz,
  PathsGraphs = NULL, OutputType = "GIRF", Economy, Folder2save = NULL,
  verbose = FALSE
)
#> Warning: number of rows of result is not a multiple of vector length (arg 1)
#> Warning: number of rows of result is not a multiple of vector length (arg 1)
#> Warning: number of rows of result is not a multiple of vector length (arg 1)
#> Warning: number of rows of result is not a multiple of vector length (arg 1)
#> Warning: number of rows of result is not a multiple of vector length (arg 1)
#> Warning: number of rows of result is not a multiple of vector length (arg 1)
```
