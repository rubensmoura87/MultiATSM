# Build the list of IRF and GIRF for both factors and bond yields

Build the list of IRF and GIRF for both factors and bond yields

## Usage

``` r
BuildIRFlist(NumOut, Economies, ModelType, IRFhoriz, FacDim, OutputType)
```

## Arguments

- NumOut:

  list of computed outputs containing the model fit, IRFs, FEVDs, GIRFs,
  GFEVDs and Term Premia

- Economies:

  Economies of the economic system

- ModelType:

  Desired model type

- IRFhoriz:

  time-horizon of the IRF and GIRF

- FacDim:

  dimension of the risk factor vector

- OutputType:

  available option are 'IRF' and 'GIRF'
