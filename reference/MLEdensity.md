# Compute the maximum likelihood function of all models

Compute the maximum likelihood function of all models

## Usage

``` r
MLEdensity(
  K1XQ,
  r0,
  SSZ,
  K0Z,
  K1Z,
  se,
  Gy.0,
  mat,
  Y,
  Z,
  P,
  Wpca,
  We,
  WpcaFull,
  dt,
  Economies,
  FactorLabels,
  ModelType,
  GVARinputs = NULL,
  JLLinputs = NULL,
  BS_outputs = FALSE,
  ExportListOut = TRUE
)
```

## Arguments

- K1XQ:

  risk-neutral feedback matrix (N x N or CN x CN)

- r0:

  long-run interest rate (scalar or vector with length C)

- SSZ:

  variance-covariance matrix (F x F)

- K0Z:

  intercept from the P-dynamics (F x 1)

- K1Z:

  feedback matrix from the P-dynamics (F x F)

- se:

  Variance of the portfolio of yields observed with error (scalar).
  Default is set to NULL.

- Gy.0:

  matrix of contemporaneous terms from the P-dynamics (F x F)

- mat:

  vector of maturities (in years) of yields used in estimation (J x 1)

- Y:

  matrix of yields used in estimation (J x T or CJ x T)

- Z:

  complete set of spanned and unspanned factors (F x T)

- P:

  complete set of spanned factors (N x T or CN x T)

- Wpca:

  matrix of weights of the portfolios observed without errors (N x J or
  CN x J)

- We:

  matrix of weights of the portfolios observed with errors ((J-N) x J or
  C(J-N) x CJ)

- WpcaFull:

  composite matrix of weights the portfolios observed with and without
  errors

- dt:

  time interval unit of the model (scalar). For instance, if data is (i)
  monthly, dt \<- 12; (ii) quarterly, dt \<- 4; (iii) yearly, dt \<- 1.

- Economies:

  string-vector containing the names of the economies which are part of
  the economic system

- FactorLabels:

  string-list based which contains the labels of all the variables
  present in the model

- ModelType:

  string-vector containing the label of the model to be estimated

- GVARinputs:

  if the model chosen is the "GVAR single" or "GVAR multi", the
  "GVARinputs" should be specified (see "GVAR" function)

- JLLinputs:

  if the model chosen is JLL-based. "JLLinputs" should contain (i)
  DomUnit, (ii) WishSigmas, (iii) SigmaNonOrtho, (iv) JLLModelType (See
  JLL function)

- BS_outputs:

  Generates simplified output list in the bootstrap setting. Default is
  set to FALSE.

- ExportListOut:

  export the complete ATSM outputs. Default is TRUE.

## References

- Candelon, C. and Moura, R. (2024). “A Multicountry Model of the Term
  Structures of Interest Rates with a GVAR.” Journal of Financial
  Econometrics 22 (5): 1558–87.

- Jotikasthira, C; Le, A. and Lundblad, C (2015). “Why Do Term
  Structures in Different Currencies Co-Move?” Journal of Financial
  Economics 115: 58–83.

- Joslin, S,; Priebsch, M. and Singleton, K. (2014). “Risk Premiums in
  Dynamic Term Structure Models with Unspanned Macro Risks.” Journal of
  Finance 69 (3): 1197–1233.

- Joslin, S., Singleton, K. and Zhu, H. (2011). "A new perspective on
  Gaussian dynamic term structure models". The Review of Financial
  Studies.

- Le, A. and Singleton, K. (2018). "A Small Package of Matlab Routines
  for the Estimation of Some Term Structure Models." Euro Area Business
  Cycle Network Training School - Term Structure Modelling.
