# Gather data of several countries in a list. Particularly useful for GVAR-based setups (Compute "GVARFactors")

Gather data of several countries in a list. Particularly useful for
GVAR-based setups (Compute "GVARFactors")

## Usage

``` r
DatabasePrep(
  t_First,
  t_Last,
  Economies,
  N,
  FactorLabels,
  ModelType,
  Macro_FullData,
  Yields_FullData,
  Wgvar = NULL
)
```

## Arguments

- t_First:

  character. Start date of the sample period in the format yyyy-mm-dd.

- t_Last:

  character. End date of the sample period in the format yyyy-mm-dd.

- Economies:

  character vector. Names of the `C` economies included in the system.

- N:

  positive integer. Number of country-specific spanned factors per
  country.

- FactorLabels:

  list. Labels for all variables present in the model, as returned by
  [`LabFac`](https://rubensmoura87.github.io/MultiATSM/reference/LabFac.md).

- ModelType:

  character. Model type to be estimated. Permissible choices: "JPS
  original", "JPS global", "GVAR single", "JPS multi", "GVAR multi",
  "JLL original", "JLL No DomUnit", "JLL joint Sigma".

- Macro_FullData:

  list. Full set of macroeconomic data, as returned by
  [`Load_Excel_Data`](https://rubensmoura87.github.io/MultiATSM/reference/Load_Excel_Data.md).

- Yields_FullData:

  list. Full set of bond yield data, as returned by
  [`Load_Excel_Data`](https://rubensmoura87.github.io/MultiATSM/reference/Load_Excel_Data.md).

- Wgvar:

  GVAR transition matrix. For GVAR models, either a matrix (`C x C`) for
  fixed weights, or a named list of matrices for time-varying weights.
  Default is NULL. Required for GVAR models.

## Value

List containing the risk factor set for all countries and global
factors. Particularly useful for GVAR-based models.

## General Notation

- `C`: number of countries in the system.

- `N`: number of country-specific spanned factors.

## Examples

``` r
# Load data from excel
macro_data <- Load_Excel_Data(system.file("extdata", "MacroData.xlsx", package = "MultiATSM"))
yields_data <- Load_Excel_Data(system.file("extdata", "YieldsData.xlsx", package = "MultiATSM"))
trade_data <- Load_Excel_Data(system.file("extdata", "TradeData.xlsx", package = "MultiATSM"))
#> New names:
#> • `` -> `...1`
#> New names:
#> • `` -> `...1`
#> New names:
#> • `` -> `...1`
#> New names:
#> • `` -> `...1`
#> New names:
#> • `` -> `...1`

# Adjust trade data
trade_data <- lapply(trade_data, function(df) {
  countries <- df[[1]]
  df <- as.data.frame(df[-1])
  rownames(df) <- countries
  df
})

# Define features of interest
ModelType <- "GVAR multi"
Economies <- c("China", "Uruguay", "Russia")
GlobalVar <- c("GBC", "CPI_OECD")
DomVar <- c("Eco_Act", "Inflation")
N <- 3
t0 <- "2006-09-01"
tF <- "2019-01-01"


# Compute some inputs
FactorLabels <- LabFac(N, DomVar, GlobalVar, Economies, ModelType)
Wgvar <- Transition_Matrix(
  t_First = "2006", t_Last = "2019", Economies,
  type = "Sample Mean", trade_data
)

# Compute GVARFactors
GVARFactors <- DatabasePrep(
  t0, tF, Economies, N, FactorLabels, ModelType, macro_data,
  yields_data, Wgvar
)
```
