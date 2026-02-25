# Main Jacobian approximation

Main Jacobian approximation

## Usage

``` r
Jac_approx(f, x, nlevels = 5L)
```

## Arguments

- f:

  function to differentiate

- x:

  evaluation point

- nlevels:

  depth of extrapolation. Default is 5.

## References

Le, A., & Singleton, K. J. (2018). Small Package of Matlab Routines for
Estimation of Some Term Structure Models. EABCN Training School.  
This function offers an independent R implementation that is informed by
the conceptual framework outlined in Le and Singleton (2018), but
adapted to the present modeling context.
