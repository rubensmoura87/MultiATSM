# Transform B_spanned into B_unspanned for sepQ models

Transform B_spanned into B_unspanned for sepQ models

## Usage

``` r
BUnspannedAdapSep(G, M, ModelPara_Short, Economies, Economy, ModelType)
```

## Arguments

- G:

  number of global unspanned factors

- M:

  number of domestic unspanned factors per country

- ModelPara_Short:

  list of model parameter estimates (See the "Optimization" function)

- Economies:

  complet set of economies of the economic system

- Economy:

  specific economy under study

- ModelType:

  a string-vector containing the label of the model to be estimated
