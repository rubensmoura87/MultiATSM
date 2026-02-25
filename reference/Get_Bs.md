# Build the B loadings

Build the B loadings

## Usage

``` r
Get_Bs(mat, dt, K1XQ, SSZ, Wpca, FactorLabels, Economy, ModelType)
```

## Arguments

- mat:

  vector of maturities (in years) of yields used in estimation (J x 1)

- dt:

  time interval unit of the model (scalar). For instance, if data is (i)
  monthly, dt \<- 12; (ii) quarterly, dt \<- 4; (iii) yearly, dt \<- 1

- K1XQ:

  risk-neutral feedback matrix (N x N or CN x CN)

- SSZ:

  variance-covariance matrix (F x F)

- Wpca:

  matrix of weights of the portfolios observed without errors (N x J or
  CN x J)

- FactorLabels:

  string-list based which contains the labels of all the variables
  present in the model

- Economy:

  string-vector containing the names of the economies which are part of
  the economic system

- ModelType:

  string-vector containing the label of the model to be estimated
