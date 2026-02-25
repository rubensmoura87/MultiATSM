# Get the intercept, feedback matrix and the variance-covariance matrix from GVAR without global factors

Get the intercept, feedback matrix and the variance-covariance matrix
from GVAR without global factors

## Usage

``` r
Get_G0G1Sigma(ParaVARX, GVARinputs, Ai0, Ai1, Wi)
```

## Arguments

- ParaVARX:

  Set of VARX model parameters

- GVARinputs:

  List of inputs for GVAR-based models

- Ai0:

  list containing the country-specific intercepts

- Ai1:

  list containing the country-specific feedback matrices

- Wi:

  list containing the country-specific link matrices
