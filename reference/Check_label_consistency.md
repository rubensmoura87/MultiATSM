# Check consistency of labels (economies, domestic and global variables)

Check consistency of labels (economies, domestic and global variables)

## Usage

``` r
Check_label_consistency(
  Economies,
  DomesticMacroFac,
  GlobalMacroFac,
  DomVarLab,
  GloVarLab
)
```

## Arguments

- Economies:

  string-vector containing the names of the economies which are part of
  the economic system

- DomesticMacroFac:

  time series of the country-specific domestic risk factors (C(M+N) x T)

- GlobalMacroFac:

  time series of the global risk factors (G x T)

- DomVarLab:

  string-vector containing the names of the desired domestic risk
  factors

- GloVarLab:

  string-vector containing the names of the desired global risk factors
