# Retrieves data from Excel and builds the database used in the model estimation

Retrieves data from Excel and builds the database used in the model
estimation

## Usage

``` r
DataForEstimation(
  t0,
  tF,
  Economies,
  N,
  FactorLabels,
  ModelType,
  DataFrequency,
  Macro_FullData,
  Yields_FullData,
  DataConnect = NULL,
  W_type = NULL,
  t_First_Wgvar = NULL,
  t_Last_Wgvar = NULL
)
```

## Arguments

- t0:

  character. Start date of the sample period in the format yyyy-mm-dd.

- tF:

  character. End date of the sample period in the format yyyy-mm-dd.

- Economies:

  character vector. Names of the `C` economies included in the system.

- N:

  positive integer. Number of country-specific spanned factors.

- FactorLabels:

  list. Labels for all variables present in the model, as returned by
  [`LabFac`](https://rubensmoura87.github.io/MultiATSM/reference/LabFac.md).

- ModelType:

  character. Model type to be estimated. Permissible choices: "JPS
  original", "JPS global", "GVAR single", "JPS multi", "GVAR multi",
  "JLL original", "JLL No DomUnit", "JLL joint Sigma".

- DataFrequency:

  character. Data frequency. Permissible choices: "Daily All Days",
  "Daily Business Days", "Weekly", "Monthly", "Quarterly", "Annually".

- Macro_FullData:

  list. Full set of macroeconomic data.

- Yields_FullData:

  list. Full set of bond yield data.

- DataConnect:

  list. Data for computing bilateral connectedness measures. Default is
  NULL. Required for GVAR-based models.

- W_type:

  character. Weight matrix type. Permissible choices: "Full Sample" (all
  years), "Sample Mean" (average over sample), or a specific year (e.g.
  "1998", "2005"). Default is NULL.

- t_First_Wgvar:

  character. First year for weight matrix computation. Default is NULL.

- t_Last_Wgvar:

  character. Last year for weight matrix computation. Default is NULL.

## Value

A list containing:

1.  Yields: matrix (`J x Td` or `CJ x Td`) of bond yields for all
    countries.

2.  RiskFactors: matrix (`K x Td`) of risk factors for all countries.

3.  GVARFactors: list of variables used in VARX estimation (see
    `GVARFactors` data file). NULL if not GVAR-based.

## General Notation

- `Td`: model time series dimension.

- `C`: number of countries in the system.

- `N`: number of country-specific spanned factors.

- `K`: total number of risk factors.

- `J`: number of bond yields per country used in estimation.

## See also

[`Load_Excel_Data`](https://rubensmoura87.github.io/MultiATSM/reference/Load_Excel_Data.md)

## Examples

``` r
DomVar <- c("Eco_Act", "Inflation")
GlobalVar <- c("GBC", "CPI_OECD")
t0 <- "2006-09-01"
tF <- "2019-01-01"
Economies <- c("China", "Brazil", "Mexico", "Uruguay", "Russia")
N <- 2
ModelType <- "JPS original"
FactorLabels <- LabFac(N, DomVar, GlobalVar, Economies, ModelType)
DataFrequency <- "Monthly"
MacroData <- Load_Excel_Data(system.file("extdata", "MacroData.xlsx", package = "MultiATSM"))
YieldData <- Load_Excel_Data(system.file("extdata", "YieldsData.xlsx", package = "MultiATSM"))
DataModel <- DataForEstimation(
  t0, tF, Economies, N, FactorLabels, ModelType, DataFrequency,
  MacroData, YieldData
)
```
