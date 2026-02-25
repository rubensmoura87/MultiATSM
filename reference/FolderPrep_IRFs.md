# Create folders for storing IRFs and GIRFs

Create folders for storing IRFs and GIRFs

## Usage

``` r
FolderPrep_IRFs(
  OutputType,
  WishPdynamicsgraphs,
  WishYieldsgraphs,
  Economies,
  ModelType,
  Folder2Save
)
```

## Arguments

- OutputType:

  available options are "IRF", "GIRF", "IRF Ortho" and "GIRF Ortho"

- WishPdynamicsgraphs:

  binary variable specifing whether the user whishes IRFs and/or GIRFs
  of risk factors

- WishYieldsgraphs:

  binary variable specifing whether the user whishes IRFs and/or GIRFs
  of bond yields

- Economies:

  Set of economies that are part of the economic system

- ModelType:

  Desired modem type

- Folder2Save:

  Folder path where the outputs will be stored.
