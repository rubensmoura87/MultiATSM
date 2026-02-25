# Data: Risk Factors for the GVAR - Candelon and Moura (2023)

Domestic risk factors data used in the GVAR models - Candelon and Moura
(2023)

## Usage

``` r
data("DomMacro_covid")
```

## Format

A matrix of country-specific risk factors (inflation, output growth,
CDS, and COVID-19 reproduction rate) for Brazil, India, Mexico, and
Russia. The data have weekly frequency and span the period from March
22, 2020, to September 26, 2021.

## Source

- Inflation:

  Monthly CPI (from OECD) interpolated to daily data (spline), converted
  to weekly year-over-year changes, and detrended
  \<https://www.oecd.org/en/data/indicators/inflation-cpi.html\>

- Output growth:

  Detrended weekly estimate of GDP year-over-year growth derived from
  the OECD Weekly Tracker index
  \<https://web-archive.oecd.org/sections/weekly-tracker-of-gdp-growth/index.htm\>

- CDS:

  5-year maturity CDS. Simulated data constructed using Bloomberg bond
  yield series.

- COVID-19 reproduction rate:

  detrended R rate from the Our World in Data database
  \<https://ourworldindata.org/coronavirus\>

## References

Candelon, B. and Moura, R. (2023) "Sovereign yield curves and the
COVID-19 in emerging markets". (Economic Modelling)
