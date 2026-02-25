# Compute the Feedback matrix of each bootstrap draw

Compute the Feedback matrix of each bootstrap draw

## Usage

``` r
FeedbackMat_BS(
  ModelType,
  RiskFactors_TS,
  FactorLabels,
  Economies,
  GVARlist,
  JLLlist,
  WishBRW,
  BRWlist,
  verbose
)
```

## Arguments

- ModelType:

  String-vector containing the label of the model to be estimated

- RiskFactors_TS:

  Time-series of risk factors of the bootstrap (F x T)

- FactorLabels:

  String-list based which contains the labels of all the variables
  present in the model

- Economies:

  String-vector containing the names of the economies which are part of
  the economic system

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
  (see
  [`Bias_Correc_VAR`](https://rubensmoura87.github.io/MultiATSM/reference/Bias_Correc_VAR.md)
  function)

- verbose:

  Logical flag controlling function messaging.
