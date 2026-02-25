# Constructs the model numerical outputs (model fit, IRFs, GIRFs, FEVDs, GFEVDs, and term premia)

Constructs the model numerical outputs (model fit, IRFs, GIRFs, FEVDs,
GFEVDs, and term premia)

## Usage

``` r
NumOutputs(
  ModelType,
  ModelPara,
  InputsForOutputs,
  FactorLabels,
  Economies,
  Folder2save = NULL,
  verbose = TRUE
)
```

## Arguments

- ModelType:

  character. Model type to be estimated. Permissible choices: "JPS
  original", "JPS global", "GVAR single", "JPS multi", "GVAR multi",
  "JLL original", "JLL No DomUnit", "JLL joint Sigma".

- ModelPara:

  list. Point estimates of the model parameters. See outputs from
  [`Optimization`](https://rubensmoura87.github.io/MultiATSM/reference/Optimization.md)

- InputsForOutputs:

  list. Inputs for generating IRFs, GIRFs, FEVDs, GFEVDs, and Term
  Premia.

- FactorLabels:

  list. Labels for all variables present in the model, as returned by
  [`LabFac`](https://rubensmoura87.github.io/MultiATSM/reference/LabFac.md).

- Economies:

  character vector. Names of the `C` economies included in the system.

- Folder2save:

  Folder path where the outputs will be stored. Default option saves the
  outputs in a temporary directory.

- verbose:

  Logical flag controlling function messaging. Default is TRUE.

## Value

An object of class 'ATSMNumOutputs' containing the following keys
elements:

1.  Model parameter estimates

2.  Model fit of bond yields

3.  IRFs

4.  FEVDs

5.  GIRFs

6.  GFEVDs

7.  Bond yield decomposition

## Details

Both IRFs and FEVDs are computed using the Cholesky decomposition
method. The risk factors are ordered as follows: (i) global unspanned
factors, and (ii) domestic unspanned and spanned factors for each
country. The order of countries follows the sequence defined in the
`Economies` vector.

## Available methods

\- \`autoplot(object, type)\`

## References

Pesaran, H. Hashem, and Shin, Yongcheol. "Generalized impulse response
analysis in linear multivariate models." Economics letters 58.1 (1998):
17-29.

## Examples

``` r
data("ParaSetEx")
data("InpForOutEx")
# Adjust inputs according to the loaded features
ModelType <- "JPS original"
Economy <- "Brazil"
FacLab <- LabFac(N = 1, DomVar = "Eco_Act", GlobalVar = "Gl_Eco_Act", Economy, ModelType)

NumOut <- NumOutputs(ModelType, ParaSetEx, InpForOutEx, FacLab, Economy,
  Folder2save = NULL, verbose = FALSE
)
```
