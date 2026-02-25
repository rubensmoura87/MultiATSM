# Build the time-series of the risk factors in each bootstrap draw

Build the time-series of the risk factors in each bootstrap draw

## Usage

``` r
BuildRiskFactors_BS(
  ModelParaPE,
  residPdynOriginal,
  residYieOriginal,
  InputsForOutputs,
  Economies,
  ModelType,
  FactorLabels,
  GVARlist,
  JLLlist,
  WishBRW,
  BRWlist,
  nlag = 1,
  verbose
)
```

## Arguments

- ModelParaPE:

  List of point estimates of the model parameter

- residPdynOriginal:

  Time-series of the residuals from the P-dynamics equation (T x F)

- residYieOriginal:

  Time-series of the residuals from the observational equation (T x J or
  T x CJ)

- InputsForOutputs:

  List containing the desired inputs for the construction of the
  numerical outputs.

- Economies:

  String-vector containing the names of the economies which are part of
  the economic system

- ModelType:

  Desired model to be estimated

- FactorLabels:

  String-list based which contains the labels of all the variables
  present in the model

- GVARlist:

  List of necessary inputs for the estimation of GVAR-based models

- JLLlist:

  List of necessary inputs for the estimation of JLL-based models

- WishBRW:

  Whether the user wishes to estimate the physical parameter model with
  the Bias correction model from BRW (2012) (see "Bias_Correc_VAR"
  function). Default is set to 0.

- BRWlist:

  List of necessary inputs for performing the bias-corrected estimation
  (see "Bias_Correc_VAR" function)

- nlag:

  Number of lags in the P-dynamics. Default is set to 1.

- verbose:

  Logical flag controlling function messaging.
