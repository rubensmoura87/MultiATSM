# Generates the desired bootstrap graphs

Generates the desired bootstrap graphs

## Usage

``` r
Boot_DataGraphFact_perShock(
  NumOutBounds,
  NumOutPE,
  IdxShock,
  nmResponse,
  Lab_Int,
  ModelType,
  FacDim,
  Horiz,
  Economies = NULL,
  Ortho = FALSE
)
```

## Arguments

- NumOutBounds:

  numerical output set from the bootstrap analysis

- NumOutPE:

  numerical output set from the point estimate analysis

- IdxShock:

  index associated with the shock variable

- nmResponse:

  Label of the response variable

- Lab_Int:

  Output types "IRF", "GIRF" and "IRF Ortho"

- ModelType:

  desired model type

- FacDim:

  dimension from the P-dynamics

- Horiz:

  horizon of analysis

- Economies:

  name of economies forming the economic system

- Ortho:

  Option for orthogonal outputs, for JLL models. Default is FALSE.
