# Loads data sets from several papers

Loads data sets from several papers

## Usage

``` r
LoadData(DataPaper)
```

## Arguments

- DataPaper:

  Available options are `BR_2017` (Bauer and Rudebusch, 2017) ,
  `CM_2023` (Candelon and Moura, 2023), `CM_2024` (Candelon and Moura,
  2024)

## Value

Complete set of data from several papers.

## References

1.  Bauer and Rudebusch (2017). "Resolving the Spanning Puzzle in
    Macro-Finance Term Structure Models" (Review of Finance)

2.  Candelon and Moura (2023). "Sovereign yield curves and the COVID-19
    in emerging markets" (Economic Modelling)

3.  Candelon and Moura (2024). "A Multicountry Model of the Term
    Structures of Interest Rates with a GVAR" (Journal of Financial
    Econometrics)

## Examples

``` r
# Example 1:
LoadData("BR_2017")

# Example 2:
LoadData("CM_2023")

# Example 3:
LoadData("CM_2024")
```
