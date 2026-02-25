# Adjust delta for numerical differentiation

Adjust delta for numerical differentiation

## Usage

``` r
adjust_delta(f, x, delta, i, fx0, direction)
```

## Arguments

- f:

  function which contains vector (J x T) valued function handle

- x:

  parameter values

- delta:

  initial delta value

- i:

  index of the parameter being adjusted

- fx0:

  initial function value

- direction:

  direction of adjustment (1 for positive, -1 for negative)

## Value

adjusted delta value
