# Generates the labels for risk factors used in the model

Generates the labels for risk factors used in the model

## Usage

``` r
LabFac(N, DomVar, GlobalVar, Economies, ModelType)
```

## Arguments

- N:

  positive integer. Number of country-specific spanned factors. Must be
  between 1 and 8.

- DomVar:

  character vector. Names of the domestic variables.

- GlobalVar:

  character vector. Names of the global variables.

- Economies:

  character vector. Names of the economies included in the system.

- ModelType:

  character. Model type to be estimated. Permissible choices: "JPS
  original", "JPS global", "GVAR single", "JPS multi", "GVAR multi",
  "JLL original", "JLL No DomUnit", "JLL joint Sigma".

## Value

List containing the risk factor labels for spanned, domestic, star, and
global variables, as well as tables for each country and all countries.

## Examples

``` r
N <- 2
DomVar <- c("inflation", "Output gap")
GlobalVar <- "Commodity Prices"
Economies <- c("U.S.", "Canada", "Germany", "Japan")
ModelType <- "JPS original"
VarLabels <- LabFac(N, DomVar, GlobalVar, Economies, ModelType)
```
