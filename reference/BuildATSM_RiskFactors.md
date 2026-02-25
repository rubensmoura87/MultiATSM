# Builds the time series of the risk factors that are used in the estimation of the ATSM

Builds the time series of the risk factors that are used in the
estimation of the ATSM

## Usage

``` r
BuildATSM_RiskFactors(
  InitialSampleDate,
  FinalSampleDate,
  Yields,
  GlobalMacroFactors,
  DomMacroFactors,
  Economies,
  FactorLabels,
  ModelType,
  BS_Adj = FALSE
)
```

## Arguments

- InitialSampleDate:

  Sample starting date

- FinalSampleDate:

  Sample last date

- Yields:

  matrix (J x T), where J - the number of maturities and T - time series
  length

- GlobalMacroFactors:

  time series of the global macroeconomic risk factors (G x T)

- DomMacroFactors:

  time series of the country-specific macroeconomic risk factors for all
  C countries (CM x T)

- Economies:

  string-vector containing the names of the economies which are part of
  the economic system

- FactorLabels:

  string-list based which contains the labels of all variables present
  in the model

- ModelType:

  string-vector containing the label of the model to be estimated

- BS_Adj:

  adjustment of global series for sepQ model in the Bootstrap setting.
  Default is set to FALSE.

## Value

Generates the complete set of risk factors that are used in the
estimation of the ATSM
