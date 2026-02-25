# Estimate numerically the variance-covariance matrix from the GVAR-based models

Estimate numerically the variance-covariance matrix from the GVAR-based
models

## Usage

``` r
EstimationSigma_GVARrest(SigmaUnres, res, IdxVarRest)
```

## Arguments

- SigmaUnres:

  Unrestricted variance-covariance matrix (K x K)

- res:

  residuals from the VAR of a GVAR model (K x T)

- IdxVarRest:

  index of the variable that is selected as strictly exogenous

## Value

restricted version of the variance-covariance matrix a GVAR model (K x
K)
