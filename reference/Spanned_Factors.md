# Computes the country-specific spanned factors

Computes the country-specific spanned factors

## Usage

``` r
Spanned_Factors(Yields, Economies, N)
```

## Arguments

- Yields:

  matrix (`J x Td`). Bond yields for all countries.

- Economies:

  character vector. Names of the `C` economies included in the system.

- N:

  integer. Desired number of country-specific spanned factors (maximum
  allowed is `N = J`).

## Value

matrix. Contains the `N` spanned factors for all countries in the system
(`CJ x Td`).

## General Notation

- `Td`: model time series dimension

- `C`: number of countries in the system

- `N`: number of country-specific spanned factors

- `J`: number of bond yields per country used in estimation

## Examples

``` r
data(Yields)
Economies <- c("China", "Brazil", "Mexico", "Uruguay")
N <- 3
SpaFact_TS <- Spanned_Factors(Yields, Economies, N)
```
