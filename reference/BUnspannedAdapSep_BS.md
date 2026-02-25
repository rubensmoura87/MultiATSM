# Obtain the full form of B unspanned for "sep Q" models within the bootstrap setting

Obtain the full form of B unspanned for "sep Q" models within the
bootstrap setting

## Usage

``` r
BUnspannedAdapSep_BS(G, M, ModelParaBoot, Economies, Economy, ModelType, tt)
```

## Arguments

- G:

  number of global unspanned factors

- M:

  number of country-specific domestic unspanned factors

- ModelParaBoot:

  list of model parameter estimates (see the "Optimization" function)
  after a bootstrap draw

- Economies:

  string-vector containing the names of the economies which are part of
  the economic system

- Economy:

  string-vector containing the names of the economy under study

- ModelType:

  string-vector containing the label of the model to be estimated

- tt:

  number of the bootstrap draw
