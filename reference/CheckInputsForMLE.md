# Check consistence of inputs

Check consistence of inputs

## Usage

``` r
CheckInputsForMLE(
  t0,
  tF,
  Economies,
  DomesticMacroFac,
  GlobalMacroFac,
  UnitYields,
  DataFreq,
  Label_Single_Models,
  Label_Multi_Models,
  FactorLabels,
  GVARlist,
  ModelType
)
```

## Arguments

- t0:

  Sample starting date

- tF:

  Sample last date

- Economies:

  string-vector containing the names of the economies of the system.

- DomesticMacroFac:

  time series of the country-specific macroeconomic risk factors for all
  C countries (CM x T)

- GlobalMacroFac:

  time series of the global macroeconomic risk factors (G x T)

- UnitYields:

  \(i\) "Month": if maturity of yields are expressed in months or (ii)
  "Year": if maturity of yields are expressed in years

- DataFreq:

  single element character-based vector. Available options are: "Daily
  All Days",  
  "Daily Business Days", "Weekly", "Monthly", "Quarterly", "Annually"

- Label_Single_Models:

  string-vector containing the names of the single country setups

- Label_Multi_Models:

  string-vector containing the names of the multicountry setups

- FactorLabels:

  string-list based which contains the labels of all variables present
  in the model

- GVARlist:

  list of necessary inputs for the estimation of GVAR-based models (see
  "GVAR" function)

- ModelType:

  string-vector containing the label of the model to be estimated
