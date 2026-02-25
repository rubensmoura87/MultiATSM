# Richardson extrapolation

Richardson extrapolation

## Usage

``` r
richardson_diff(f, x, i, h, fx0, nlevels = 5L)
```

## Arguments

- f:

  function to differentiate

- x:

  evaluation point

- i:

  coordinate index

- h:

  base step size

- fx0:

  f(x), precomputed

- nlevels:

  depth of extrapolation. Default is 5
