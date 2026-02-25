# Perform the minimization of ML function

Perform the minimization of ML function

## Usage

``` r
Optimization_PE(
  ML_fun,
  ListInputSet,
  FactorLabels,
  Economies,
  ModelType,
  JLLinputs = NULL,
  GVARinputs = NULL,
  tol = 1e-04,
  EstType,
  TimeCount = TRUE,
  verbose
)
```

## Arguments

- ML_fun:

  vector-valued objective ML function

- ListInputSet:

  list containing :

  1.  a starting value for K1XQ and/or SSZ

  2.  a variable label among the following:

      - 'Jordan' or 'Jordan; stationary' for single countries setups or
        'Jordan MultiCountry' or 'Jordan MultiCountry; stationary'. All
        cases related to the computation of a K1XQ parameter;

      - 'psd': PSD matrix, used for JPS-based models. It relates to the
        SSZ parameter;

      - 'BlockDiag': block diagonal matrix, used for JPS-based models.
        It relates to the SSZ parameter.

      - 'JLLstructure': to impose the zero-restrictions on the SSZ term
        along the lines of the JLL models

- FactorLabels:

  list. Labels for all variables present in the model, as returned by
  [`LabFac`](https://rubensmoura87.github.io/MultiATSM/reference/LabFac.md).

- Economies:

  character vector. Names of the `C` economies included in the system.

- ModelType:

  character. Model type to be estimated. Permissible choices: "JPS
  original", "JPS global", "GVAR single", "JPS multi", "GVAR multi",
  "JLL original", "JLL No DomUnit", "JLL joint Sigma".

- JLLinputs:

  List. Inputs for JLL model estimation (see
  [`JLL`](https://rubensmoura87.github.io/MultiATSM/reference/JLL.md)).
  Default is NULL.

- GVARinputs:

  List. Inputs for GVAR model estimation (see
  [`GVAR`](https://rubensmoura87.github.io/MultiATSM/reference/GVAR.md)).
  Default is NULL.

- tol:

  convergence tolerance (scalar). Default value is 1e-4.

- EstType:

  Available options are"BFGS" and/or "Nelder-Mead".

- TimeCount:

  computes the required time for estimation of the model. Default is
  TRUE.

- verbose:

  Logical flag controlling function messaging.
