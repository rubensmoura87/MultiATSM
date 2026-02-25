# Compute the variance-covariance matrix after the bias correction procedure

Compute the variance-covariance matrix after the bias correction
procedure

## Usage

``` r
Get_SSZ_BC(K1Z_BC, RiskFactors, GVARinputs, JLLinputs, FactorLabels, ModelType)
```

## Arguments

- K1Z_BC:

  Feedback matrix resulting from the bias-correction procedure

- RiskFactors:

  time series of the risk factors (T x F)

- GVARinputs:

  inputs used in the estimation of the GVAR-based models (see "GVAR"
  function). Default is set to NULL

- JLLinputs:

  inputs used in the estimation of the JLL-based models (see "JLL"
  function). Default is set to NULL

- FactorLabels:

  string-list based which contains the labels of all variables present
  in the model

- ModelType:

  string-vector containing the label of the model to be estimated
