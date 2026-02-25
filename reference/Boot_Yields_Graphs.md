# Build P-dynamic graphs after the bootstrap implementation

Build P-dynamic graphs after the bootstrap implementation

## Usage

``` r
Boot_Yields_Graphs(
  NumOutBounds,
  NumOutPE,
  ModelType,
  FacDim,
  YielDim,
  Horiz,
  Economies,
  PathsGraphs,
  OutInt,
  Folder2save,
  WishYieldGraphs,
  WishYieldGraphsOrtho = NULL
)
```

## Arguments

- NumOutBounds:

  numerical output set from the bootstrap analysis

- NumOutPE:

  numerical output set from the point estimate analysis

- ModelType:

  desired model type

- FacDim:

  dimension of the risk factor set

- Horiz:

  horizon of numerical outputs

- Economies:

  Economies that are part of the economic system

- PathsGraphs:

  Path to save the desired graphs

- OutInt:

  Available option are "IRF", "FEVD", "GIRF" or "GFEVD"

- Folder2save:

  Folder path where the outputs will be stored.

- WishYieldGraphs:

  Binary variable reflecting the graphs of interest

- WishYieldGraphsOrtho:

  Binary variable reflecting the graphs of interest (orthogonalized
  version). Default is NULL
