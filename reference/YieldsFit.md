# Computes two measures of model fit for bond yields (all models)

Computes two measures of model fit for bond yields (all models)

## Usage

``` r
YieldsFit(ModelType, ModelPara, FactorLabels, Economies)
```

## Arguments

- ModelType:

  a string-vector containing the label of the model to be estimated

- ModelPara:

  List of model parameter estimates (See the "Optimization" function)

- FactorLabels:

  a string-list based which contains the labels of all the variables
  present in the model

- Economies:

  a string-vector containing the names of the economies which are part
  of the economic system

## Details

"Model-implied yields" is the measure of fit based exclusively on the
risk-neutral parameters, whereas the "Model-Fit" takes into account both
the risk-neutral and the physical parameters.

## References

See, for instance, Jotiskhatira, Le and Lundblad (2015). "Why do
interest rates in different currencies co-move?" (Journal of Financial
Economics)
