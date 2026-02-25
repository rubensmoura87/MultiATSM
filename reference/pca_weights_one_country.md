# Computes the PCA weights for a single country

Computes the PCA weights for a single country

## Usage

``` r
pca_weights_one_country(Yields, Economy)
```

## Arguments

- Yields:

  matrix (`J x Td`). Bond yields for a single country.

- Economy:

  character. Name of the economy.

## Value

matrix (`J x J`). Eigenvectors of the variance-covariance matrix of
yields.

## General Notation

- `Td`: model time series dimension

- `J`: number of bond yields per country used in estimation

## Examples

``` r
data(Yields)
Economy <- "Mexico"
pca_weights <- pca_weights_one_country(Yields, Economy)
```
