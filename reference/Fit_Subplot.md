# Build subplot for fitted yields

Build subplot for fitted yields

## Usage

``` r
Fit_Subplot(
  YieldData,
  ModelFit,
  ModelImplied,
  MatLength,
  mat,
  Economies,
  ModelType,
  PathsGraphs
)
```

## Arguments

- YieldData:

  Time series of bond yields

- ModelFit:

  Time series of fitted bond yields

- ModelImplied:

  Time series of model-implied bond yields

- MatLength:

  number of country-specific maturities

- mat:

  vector of maturities

- Economies:

  Economies of the economic system

- ModelType:

  Desired estimated model

- PathsGraphs:

  Path to save the graphs
