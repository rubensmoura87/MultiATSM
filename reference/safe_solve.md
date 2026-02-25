# Safe matrix inversion with conditioning check

Safe matrix inversion with conditioning check

## Usage

``` r
safe_solve(A, B = NULL, reg = 1e-12, rcond_tol = 1e-12, verbose = FALSE)
```

## Arguments

- A:

  matrix

- B:

  matrix

- reg:

  Numeric scalar (default: `1e-12`). Small ridge term added to the
  diagonal of `A` when the system is ill-conditioned.

- rcond_tol:

  Numeric, default = 1e-12. Tolerance threshold for the reciprocal
  condition number.

- verbose:

  Logical flag controlling function messaging.
