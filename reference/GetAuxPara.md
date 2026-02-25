# Compute the auxiliary parameters a.

Compute the auxiliary parameters a.

## Usage

``` r
GetAuxPara(
  ParaValue,
  Const_Type_Full,
  Economies,
  FactorLabels,
  JLLinputs = NULL
)
```

## Arguments

- ParaValue:

  Parameter value

- Const_Type_Full:

  character-based vector that describes the constraints. Constraints
  are:

  - 'Jordan' for single-country models;

  - 'Jordan; stationary' for single-country models;

  - 'Jordan MultiCountry' for multicountry models;

  - 'Jordan MultiCountry; stationary' for multicountry models;

  - 'psd' for JPS-based models;

  - 'BlockDiag' for GVAR-based models;

  - 'JLLstructure' for JLL-based models;

- Economies:

  string-vector containing the names of the economies which are part of
  the economic system

- FactorLabels:

  string-list based which contains the labels of all the variables
  present in the model

- JLLinputs:

  list of necessary inputs for the estimation of JLL-based models (see
  "JLL" function)
