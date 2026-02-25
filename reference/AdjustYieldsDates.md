# Makes sure that the time series of yields and risk factors have coincident sample spans

Makes sure that the time series of yields and risk factors have
coincident sample spans

## Usage

``` r
AdjustYieldsDates(Yields, PdynamicsFactors, Economies)
```

## Arguments

- Yields:

  time series of bond yields (CJ x T1 )

- PdynamicsFactors:

  time series of risk factors (F x T2)

- Economies:

  string-vector containing the names of the economies of the system.
