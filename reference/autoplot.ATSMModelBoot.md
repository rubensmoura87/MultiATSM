# Autoplot method for ATSMModelBoot objects

Autoplot method for ATSMModelBoot objects

## Usage

``` r
# S3 method for class 'ATSMModelBoot'
autoplot(x, NumOutPE, type, ...)
```

## Arguments

- x:

  An object of class 'ATSMModelBoot'

- NumOutPE:

  An object of class 'ATSMNumOutputs': point estimates of the numerical
  outputs

- type:

  Plot type: one of "IRF", "FEVD", "GIRF", "GFEVD" (each must be
  suffixed with "\_Factors" or "\_Yields"). For JLL-based models, an
  additional "\_Ortho" suffix produces orthogonalized outputs. All
  inputs must end by "\_Boot" as a reference to the bootstrap procedure.

- ...:

  Additional arguments (not used)
