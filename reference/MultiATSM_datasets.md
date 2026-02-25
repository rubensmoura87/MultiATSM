# Overview of Datasets Included in the MultiATSM Package

The package includes several pre-processed datasets used for estimation
and replication examples:

## Details

- GlobalMacro:

  Global macro-financial risk factors, namely global economic activity
  and global inflation.

- GlobalMacro_covid:

  Global macro-financial risk factors, namely the output growth rate
  from the U.S. and China and the S&P 500 index.

- DomMacro:

  Domestic macroeconomic risk factors, namely economic activity and
  inflation.

- DomMacro_covid:

  Domestic macroeconomic risk factors, namely otput growth, inflation,
  CDS and the COVID-19 reproduction rate

- TradeFlows:

  Bilateral trade flow series used in GVAR examples as a proxy measure
  of cross-country conectdness.

- TradeFlows_covid:

  Bilateral trade flow series used in GVAR examples as a proxy measure
  of cross-country conectdness.

- Yields:

  Monthly series of bond yields by maturity for multiple economies.

- Yields_covid:

  Weekly series of sovereign bond yields by maturity for multiple
  economies.

- RiskFacFull:

  Full set of risk factors (global and domestic) data used throughout
  the package

- GVARFactors:

  List of risk factors used in the estimation of GVAR models.

- BR_jps_out:

  Replications of the JPS outputs by Bauer and Rudebusch (2017)

- InpForOutEx:

  List of inputs for an illustrative JPS model with Brazilian data

- ParaSetEx:

  List of set of parameterafter optimization for an illustrative JPS
  model with Brazilian data

- NumOutEx:

  List of numerical outputs for an illustrative JPS model with Brazilian
  data

- Out_Example:

  re-loaded examaple of a complete list of several model outputs. Used
  in the package vignette.

Each dataset is documented separately using \`?GlobalMacro\`,
\`?DomMacro\`, \`?TradeFlows\`, \`?Yields\`, etc. Datasets ending with
the suffix `_covid` are based on those used in Candelon and Moura (2023)
and cover Brazil, India, Mexico, and Russia. The remaining datasets
correspond to Candelon and Moura (2024) and include Brazil, China,
Mexico, and Uruguay.

## References

1.  Candelon, B. and Moura, R. (2023) "Sovereign yield curves and the
    COVID-19 in emerging markets". (Economic Modelling)

2.  Candelon, B. and Moura, R. (2024) "A Multicountry Model of the Term
    Structures of Interest Rates with a GVAR". (Journal of Financial
    Econometrics)
