# Data: Risk Factors for the GVAR - Candelon and Moura (2024, JFEC)

Risk factors data used in the GVAR-ATSM from Candelon and Moura (2024,
JFEC)

## Usage

``` r
data("GVARFactors")
```

## Format

List of risk factors organized for GVAR estimation. It includes global
unspanned factors (economic activity, inflation) and domestic
factors—both unspanned (economic activity, inflation) and spanned
(level, slope, curvature) with their starred counterparts. The dataset
covers Brazil, China, Mexico, and Uruguay at a monthly frequency from
June 2004 to January 2020.

## Source

- Global unspanned factors:

  See `data("GlobalMacro")` for a detailed data description.

- Domestic unspanned factors:

  See `data("DomMacro")` for a detailed data description.

- Domestic spanned factors:

  First three principal components of each country set of bond yields.
  See `data("Yields")` for a detailed data description.

- Domestic star factors:

  Weighted average of foreign factors. See
  [`Transition_Matrix`](https://rubensmoura87.github.io/MultiATSM/reference/Transition_Matrix.md)
  for the computation of weights.

## References

Candelon, B. and Moura, R. (2024) "A Multicountry Model of the Term
Structures of Interest Rates with a GVAR". (Journal of Financial
Econometrics)
