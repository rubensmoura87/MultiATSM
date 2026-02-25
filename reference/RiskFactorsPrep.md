# Builds the complete set of time series of the risk factors (spanned and unspanned)

Builds the complete set of time series of the risk factors (spanned and
unspanned)

## Usage

``` r
RiskFactorsPrep(
  FactorSet,
  Economies,
  FactorLabels,
  Initial_Date,
  Final_Date,
  DataFrequency
)
```

## Arguments

- FactorSet:

  Factor set list (see e.g. "GVARFactors" data file)

- Economies:

  A character vector containing the names of the economies included in
  the system.

- FactorLabels:

  A list of character vectors with labels for all variables in the
  model.

- Initial_Date:

  Start date of the sample period in the format yyyy-mm-dd

- Final_Date:

  End date of the sample period in the format yyyy-mm-dd

- DataFrequency:

  A character vector specifying the data frequency. Available options:
  "Daily All Days", "Daily Business Days", "Weekly", "Monthly",
  "Quarterly", "Annually".

## Value

Risk factors used in the estimation of the desired ATSM
