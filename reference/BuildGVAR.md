# Build the GVAR(1) from the country-specific VARX(1,1,1)

Build the GVAR(1) from the country-specific VARX(1,1,1)

## Usage

``` r
BuildGVAR(ParaVARX, GlobalPara, GVARinputs, DomLabels, GlobalLabels, N)
```

## Arguments

- ParaVARX:

  Set of VARX model parameters

- GlobalPara:

  Set of marginal model parameters

- GVARinputs:

  List of inputs for GVAR-based models

- DomLabels:

  string-based vector containing label of the domestic risk factors

- GlobalLabels:

  string-based vector containing label of the global risk factors

- N:

  number of country-specific spanned factors (scalar)
