# Build Confidence intervals for yield-related outputs

Build Confidence intervals for yield-related outputs

## Usage

``` r
BuildCI_Yields(
  NumOutBounds,
  NumOutPE,
  Lab_Int,
  ModelType,
  Economies,
  IdxResp,
  IdxShock,
  Ortho = FALSE
)
```

## Arguments

- NumOutBounds:

  numerical output set from the bootstrap analysis

- NumOutPE:

  numerical output set from the point estimate analysis

- Lab_Int:

  Label of interest. available options are "IRF" and "FEVD"

- ModelType:

  desired model type

- Economies:

  name of the economies forming the economic system

- IdxResp:

  index associated with the response variable

- IdxShock:

  index associated with the shock variable

- Ortho:

  Option for orthogonal outputs, for JLL models. Default is FALSE.
