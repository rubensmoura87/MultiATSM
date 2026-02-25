# Compute long-run risk neutral mean (r0) for the various models

Compute long-run risk neutral mean (r0) for the various models

## Usage

``` r
Get_r0(Y, P, N, mat, dt, B_list, Wpca, We, Economies, ModelType)
```

## Arguments

- Y:

  matrix of yields used in estimation (J x T or CJ x T)

- P:

  complete set of spanned factors (N x T or CN x T)

- N:

  number of country-specific spanned factors

- mat:

  vector of maturities (in years) of yields used in estimation (J x 1)

- dt:

  time interval unit of the model (scalar). For instance, if data is (i)
  monthly, dt \<- 12; (ii) quarterly, dt \<- 4; (iii) yearly, dt \<- 1.

- B_list:

  list containing the B loadings

- Wpca:

  matrix of weights of the portfolios observed without errors (N x J or
  CN x J)

- We:

  matrix of weights of the portfolios observed with errors ((J-N) x J or
  C(J-N) x CJ)

- Economies:

  string-vector containing the names of the economies which are part of
  the economic system

- ModelType:

  string-vector containing the label of the model to be estimated
