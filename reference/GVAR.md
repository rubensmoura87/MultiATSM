# Estimates a GVAR(1) and VARX(1,1,1) models

Estimates a GVAR(1) and VARX(1,1,1) models

## Usage

``` r
GVAR(GVARinputs, N, CheckInputs = FALSE)
```

## Arguments

- GVARinputs:

  list. Inputs for GVAR model estimation:

  1.  `Economies`: character vector. Contains the `C` names of the
      economies included in the system.

  2.  `GVARFactors`: list. All variables used in the estimation of the
      VARX model  
      (see e.g. `GVARFactors` file for details);

  3.  `VARXtype`: Permissible:

      - `'unconstrained'`: model is estimated without constraints (each
        equation is estimated individually by ordinary least square);

      - `'constrained: Spanned Factors'`: The model is estimated with
        the restriction that foreign pricing factors do NOT affect (i)
        domestic economic variables and (ii) domestic pricing factors
        (estimation via restricted least squares).

      - `'constrained : [factor_name]'`: The model is estimated with the
        restriction that the specified risk factor is influenced only by
        its own lagged values and the lagged values of its corresponding
        star variables. (estimation via restricted least squares.)

  4.  `Wgvar`: The GVAR transition matrix (`C x C`) used in the model
      solution.  
      (See the output from the
      [`Transition_Matrix`](https://rubensmoura87.github.io/MultiATSM/reference/Transition_Matrix.md)
      function.).

- N:

  positive integer. Number of country-specific spanned factors.

- CheckInputs:

  logical. Whether to perform a prior consistency check on the inputs
  provided in `GVARinputs`. Default is FALSE.

## Value

list. Contains:

1.  parameters of the country-specific VARX(1,1,1):

    - intercept (M + N x 1)

    - phi_1 (M + N x M + N)

    - phi_1\* (M + N x M + N)

    - phi_g (M + N x M + N)

    - Sigma (M + N x G)

2.  parameters of the GVAR:

    - F0 (K x K)

    - F1 (K x K)

    - Sigma_y (K x K)

## General Notation

- `C`: number of countries in the system

- `G`: number of global unspanned factors

- `M`: number of country-specific unspanned factors

- `N`: number of country-specific spanned factors

- `K`: total number of risk factors (K = C x (N + M) + G)

## References

Chudik, A. and Pesaran, M. H. (2016). "Theory and Practice of GVAR
modelling" (Journal of Economic Surveys)

## Examples

``` r
data(GVARFactors)

GVARinputs <- list(
  Economies = c("China", "Brazil", "Mexico", "Uruguay"),
  GVARFactors = GVARFactors, VARXtype = "unconstrained"
)

GVARinputs$Wgvar <- matrix(c(
  0, 0.83, 0.86, 0.38,
  0.65, 0, 0.13, 0.55,
  0.32, 0.12, 0, 0.07,
  0.03, 0.05, 0.01, 0
), nrow = 4, ncol = 4)
N <- 3

GVARPara <- GVAR(GVARinputs, N)
```
