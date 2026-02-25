# Shrinking the largest eigenvalue

Shrinking the largest eigenvalue

## Usage

``` r
shrink_FeedMat_BC(K1Z_BC, K1Z_NoBC, ev_restr)
```

## Arguments

- K1Z_BC:

  VAR (1) bias-corrected feedback matrix from Bauer, Rudebusch and, Wu
  (2012)

- K1Z_NoBC:

  VAR (1) with no bias-corrected feedback matrix from the selected ATSM

- ev_restr:

  maximum eigenvalue desired in the feedback matrix after the adjustment
