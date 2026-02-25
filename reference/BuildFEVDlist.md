# Build the list of IRF and GIRF for both factors and bond yields

Build the list of IRF and GIRF for both factors and bond yields

## Usage

``` r
BuildFEVDlist(
  NumOut,
  Economies,
  ModelType,
  FEVDhoriz,
  FacDim,
  YieldsDim,
  OutputType
)
```

## Arguments

- NumOut:

  list of computed outputs containing the model fit, IRFs, FEVDs, GIRFs,
  GFEVDs and Term Premia

- Economies:

  Economies of the economic system

- ModelType:

  Desired model type

- FEVDhoriz:

  time-horizon of the FEVD and GFEVD

- FacDim:

  dimension of the risk factor vector

- YieldsDim:

  dimension of the model set of yields

- OutputType:

  available option are 'FEVD' and 'GFEVD'
