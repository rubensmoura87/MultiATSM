# Read data from Excel files and return a named list of data frames

Read data from Excel files and return a named list of data frames

## Usage

``` r
Load_Excel_Data(ExcelFilePath)
```

## Arguments

- ExcelFilePath:

  character. Path to the Excel file (.xlsx) to load. Must be a valid
  file path. The file can contain multiple sheets; each sheet will be
  loaded as a separate data frame in the output list.

## Value

Named list of data frames, one for each sheet in the Excel file. The
names of the list elements correspond to the sheet names.

## Details

Uses the readxl package to read all sheets from the specified Excel
file. Each sheet is returned as a data frame. The output is a named
list, with names matching the sheet names in the Excel file.

## Examples

``` r
if (!requireNamespace("readxl", quietly = TRUE)) {
  stop(
    "Please install package \"readxl\" to use this feature.",
    call. = FALSE
  )

  Load_Excel_Data(system.file("extdata", "MacroData.xlsx", package = "MultiATSM"))
  Load_Excel_Data(system.file("extdata", "YieldsData.xlsx", package = "MultiATSM"))
}
```
