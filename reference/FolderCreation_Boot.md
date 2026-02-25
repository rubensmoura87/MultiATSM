# Creates folder to store graphs generated from the bootstrap analysis

Creates folder to store graphs generated from the bootstrap analysis

## Usage

``` r
FolderCreation_Boot(
  ModelType,
  LabIRF,
  Economies,
  OutType,
  Folder2save,
  Ortho = FALSE
)
```

## Arguments

- ModelType:

  Desired model type

- LabIRF:

  Output types "IRF", "GIRF" and "IRF Ortho"

- Economies:

  economies of the economic system

- OutType:

  Available option "Factors" or "Yields

- Folder2save:

  Folder path where the outputs will be stored.

- Ortho:

  Option for orthogonal outputs, for JLL models. Default is FALSE.
