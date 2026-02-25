# Autoplot method for ATSMNumOutputs objects

Autoplot method for ATSMNumOutputs objects

## Usage

``` r
# S3 method for class 'ATSMNumOutputs'
autoplot(x, type, ...)
```

## Arguments

- x:

  An object of class 'ATSMNumOutputs'

- type:

  Plot type: "RiskFactors", "Fit", "TermPremia", or one of "IRF",
  "FEVD", "GIRF", "GFEVD" (each must be suffixed with "\_Factors" or
  "\_Yields"). For JLL-based models, an additional "\_Ortho" suffix
  produces orthogonalized outputs.

- ...:

  Additional arguments (not used)
