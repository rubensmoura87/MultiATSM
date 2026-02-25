# Compute the A loadings

Compute the A loadings

## Usage

``` r
Get_As(LoadBs, Wpca, r0, dt, Economies, ModelType)
```

## Arguments

- LoadBs:

  list containing the B loadings

- Wpca:

  matrix of weights of the portfolios observed without errors (N x J or
  CN x J)

- r0:

  long-run interest rate (scalar or vector with length C)

- dt:

  time interval unit of the model (scalar). For instance, if data is (i)
  monthly, dt \<- 12; (ii) quarterly, dt \<- 4; (iii) yearly, dt \<- 1

- Economies:

  string-vector containing the names of the economies which are part of
  the economic system

- ModelType:

  string-vector containing the label of the model to be estimated
