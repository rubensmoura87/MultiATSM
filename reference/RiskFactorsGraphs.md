# Spanned and unspanned factors plot

Spanned and unspanned factors plot

## Usage

``` r
RiskFactorsGraphs(
  ModelType,
  WishRFgraphs,
  ModelOutputs,
  Economies,
  FactorLabels,
  Folder2save,
  verbose
)
```

## Arguments

- ModelType:

  character. Estimated model type. Permissible choices: "JPS original",
  "JPS global", "GVAR single", "JPS multi", "GVAR multi", "JLL
  original", "JLL No DomUnit", "JLL joint Sigma".

- WishRFgraphs:

  logical. Set TRUE to generate graphs, FALSE otherwise.

- ModelOutputs:

  list. Model parameter estimates (see
  [`Optimization`](https://rubensmoura87.github.io/MultiATSM/reference/Optimization.md)).

- Economies:

  character vector. Names of the `C` economies included in the system.

- FactorLabels:

  list. Labels for all variables in the model.

- Folder2save:

  character. Folder path where the outputs will be stored.

- verbose:

  logical. Flag controlling function messaging.

## Available Methods

\- \`autoplot(object, type = "RiskFactors")\`

## Examples

``` r
data("ParaSetEx")
# Adapt factor labels according to the example
ModelType <- "JPS original"
Economy <- "Brazil"
FacLab <- LabFac(N = 1, DomVar = "Eco_Act", GlobalVar = "Gl_Eco_Act", Economy, ModelType)

RiskFactorsGraphs(ModelType,
  WishRFgraphs = FALSE, ParaSetEx, Economy, FacLab,
  Folder2save = NULL, verbose = FALSE
)
```
