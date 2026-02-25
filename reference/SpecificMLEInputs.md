# Concatenate the model-specific inputs in a list

Concatenate the model-specific inputs in a list

## Usage

``` r
SpecificMLEInputs(
  ModelType,
  Economies,
  RiskFactors,
  FactorLabels,
  GVARlist = NULL,
  JLLlist = NULL,
  WishBRW = 0,
  BRWlist = NULL
)
```

## Arguments

- ModelType:

  string-vector containing the label of the model to be estimated

- Economies:

  string-vector containing the names of the economies of the system

- RiskFactors:

  time series of risk factors (F x T)

- FactorLabels:

  string-list based which contains the labels of all the variables
  present in the model

- GVARlist:

  A list of required inputs to estimate the GVAR-based setups:

  1.  VARXtype string-vector containing the VARX feature (see "GVAR"
      function) (GVAR-based models)

  2.  t_First_Wgvar Sample starting date (year) (GVAR-based models)

  3.  t_Last_Wgvar Sample last date (year) (GVAR-based models)

  4.  W_type Criterion used in the computation of the star variables
      (see "Transition_Matrix" function) (GVAR-based models)

- JLLlist:

  A list of required inputs to estimate the JLL-based setups:

  1.  DomUnit name of the economy which is assigned as the dominant unit
      (JLL-based models)

  2.  WishSigmas equal to "1" if one wishes the variance-covariance
      matrices and the Cholesky factorizations (JLL-based models)

  3.  SigmaNonOrtho NULL or some F x F matrix from the
      non-orthogonalized dynamics (JLL-based models)

- WishBRW:

  Whether the user wishes to estimate the physical parameter model with
  the Bias correction model from BRW (2012) (see "Bias_Correc_VAR"
  function).  
  Default is set to 0.

- BRWlist:

  A list of required inputs to estimate the bias corrected setups of the
  type of BRW:

  1.  BiasCorrection binary variable. it takes value equal to 1 if the
      user whishes the estimates to be bias-corrected and 0, otherwise.
      (BRW model)

  2.  flag_mean flag whether mean- (TRUE) or median- (FALSE) unbiased
      estimation is desired

  3.  gamma adjustment parameter (BRW model)

  4.  N_iter number of iterations (BRW model)

  5.  N_burn number of burn-in iterations (BRW model)

  6.  B number of bootstrap samples (BRW model)

  7.  checkBRW flag whether the user wishes to perform the closeness
      check (BRW model)

  8.  B_check number of bootstrap samples for closeness check
