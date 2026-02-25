# Create the variable labels used in the estimation

Create the variable labels used in the estimation

## Usage

``` r
ParaLabelsOpt(ModelType, WishStationarityQ, MLEinputs, BS_outputs = FALSE)
```

## Arguments

- ModelType:

  a string-vector containing the label of the model to be estimated

- WishStationarityQ:

  User must set TRUE is she wishes to impose the largest eigenvalue
  under the Q to be strictly smaller than 1. Otherwise set FALSE

- MLEinputs:

  Set of inputs that are necessary to the log-likelihood function

- BS_outputs:

  Generates simplified output list in the bootstrap setting. Default is
  set to FALSE.
