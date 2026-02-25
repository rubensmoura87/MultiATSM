# Map auxiliary (unconstrained) parameters a to constrained parameters b

Map auxiliary (unconstrained) parameters a to constrained parameters b

## Usage

``` r
GetTruePara(
  ParaValue,
  Const_Type_Full,
  FactorLabels,
  Economies,
  JLLinputs = NULL,
  GVARinputs = NULL
)
```

## Arguments

- ParaValue:

  unconstrained auxiliary parameter

- Const_Type_Full:

  One of the following options:

  - 'Jordan'

  - 'Jordan; stationary'

  - 'Jordan MultiCountry'

  - 'Jordan MultiCountry; stationary'

  - 'psd';

  - 'BlockDiag'

  - 'JLLstructure'

- FactorLabels:

  string-list based which contains the labels of all the variables
  present in the model

- Economies:

  string-vector containing the names of the economies which are part of
  the economic system

- JLLinputs:

  Inputs used in the estimation of the JLL-based models

- GVARinputs:

  Inputs used in the estimation of the GVAR-based models
