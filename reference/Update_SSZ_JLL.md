# Update the variance-covariance matrix from the "JLL joint Sigma" model. Necessary for optimization

Update the variance-covariance matrix from the "JLL joint Sigma" model.
Necessary for optimization

## Usage

``` r
Update_SSZ_JLL(SSZ, Z, N, JLLinputs)
```

## Arguments

- SSZ:

  Variance-covariance matrix from JLL model

- Z:

  complete set of spanned and unspanned factors (F x T)

- N:

  number of country-specific spanned factors

- JLLinputs:

  List of inputs from JLL models
