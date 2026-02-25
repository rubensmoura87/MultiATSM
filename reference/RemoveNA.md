# Exclude series that contain NAs

Exclude series that contain NAs

## Usage

``` r
RemoveNA(YieldsData, MacroData, Economies)
```

## Arguments

- YieldsData:

  List of country-specific bond yields

- MacroData:

  List of country-specific and global economic variables

- Economies:

  string-vector containing the names of the economies which are part of
  the economic system

## Value

return the time series data that were not initially composed by NAs.
