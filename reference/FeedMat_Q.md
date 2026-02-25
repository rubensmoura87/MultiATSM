# Get an estimate for the risk-neutral (Q) feedback matrix

Get an estimate for the risk-neutral (Q) feedback matrix

## Usage

``` r
FeedMat_Q(
  Yields,
  Spa_Fac,
  Economies,
  UnitYields,
  time_step,
  check_inputs = TRUE
)
```

## Arguments

- Yields:

  matrix of bond yields (J x T)

- Spa_Fac:

  matrix of spanned factors (N x T)

- Economies:

  string-vector containing the names of the economies which are part of
  the economic system

- UnitYields:

  \(i\) "Month": if maturity of yields are expressed in months or (ii)
  "Year": if maturity of yields are expressed in years

- time_step:

  time unit of the model (scalar). For instance, if data is (i) monthly,
  dt \<- 12; (ii) quarterly, dt \<- 4; (iii) yearly, dt \<- 1.

- check_inputs:

  Perform input validation. Default is TRUE.

## References

Le, A., & Singleton, K. J. (2018). Small Package of Matlab Routines for
Estimation of Some Term Structure Models. EABCN Training School.  
This function offers an independent R implementation that is informed by
the conceptual framework outlined in Le and Singleton (2018), but
adapted to the present modeling context.
