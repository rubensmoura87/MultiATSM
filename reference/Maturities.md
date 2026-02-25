# Create a vector of numerical maturities in years

Create a vector of numerical maturities in years

## Usage

``` r
Maturities(DataYields, Economies, UnitYields)
```

## Arguments

- DataYields:

  matrix containing all yields of the system (JxT,if the model is
  single-country-based or CJxT if the model is multy-country-based)

- Economies:

  vector containing names of all the economies of the system

- UnitYields:

  \(i\) "Month": if maturity of yields are expressed in months or (ii)
  "Year": if maturity of yields are expressed in years

## Value

Vector containing all observed maturities expressed in years
