# Preliminary checks for inputs provided for the performing out-of-sample forecasting

Preliminary checks for inputs provided for the performing out-of-sample
forecasting

## Usage

``` r
ChecksOOS(t0Forecast, t0Sample, nForecasts, ForecastType, TimeLength)
```

## Arguments

- t0Forecast:

  Index of the last set of observations in the information set at the
  first forecasting round

- t0Sample:

  Index of the first set of observations in the information set at the
  first forecasting round

- nForecasts:

  Number of forecasting sets generated

- ForecastType:

  Forecast type. Available options are "Rolling" and "Expanding".

- TimeLength:

  Time-series dimension of the model
