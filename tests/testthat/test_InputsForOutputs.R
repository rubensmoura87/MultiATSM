library(testthat)
library(MultiATSM)

# 1) Test basic output structure
test_that("InputsForOutputs returns correct structure and options", {
  ModelType <- "JPS original"
  Horiz <- 100
  DesiredOutputGraphs <- c("Fit", "GIRF", "GFEVD")
  OutputLabel <- "Test"
  WishStationarityQ <- 1
  DataFrequency <- "Monthly"
  WishGraphYields <- 1
  WishGraphRiskFactors <- 0
  res <- InputsForOutputs(
    ModelType, Horiz, DesiredOutputGraphs, OutputLabel,
    WishStationarityQ, DataFrequency, WishGraphYields, WishGraphRiskFactors
  )
  expect_type(res, "list")
  expect_equal(res$`Label Outputs`, OutputLabel)
  expect_equal(res$StationaryQ, WishStationarityQ)
  expect_equal(res$DataFreq, DataFrequency)
  expect_true(ModelType %in% names(res))
  expect_true("Fit" %in% names(res[[ModelType]]))
})


# 2) Test JLL model
test_that("InputsForOutputs returns correct structure for a JLL model", {
  ModelType <- "JLL original"
  Horiz <- 10
  DesiredOutputGraphs <- c("IRF", "FEVD")
  OutputLabel <- "Test"
  WishStationarityQ <- 1
  DataFrequency <- "Monthly"
  WishGraphYields <- 1
  WishGraphRiskFactors <- 0
  WishOrthoJLLgraphs <- 1
  res <- InputsForOutputs(
    ModelType, Horiz, DesiredOutputGraphs, OutputLabel,
    WishStationarityQ, DataFrequency, WishGraphYields, WishGraphRiskFactors
  )
  expect_type(res, "list")
  expect_equal(res$`Label Outputs`, OutputLabel)
  expect_equal(res$StationaryQ, WishStationarityQ)
  expect_equal(res$DataFreq, DataFrequency)
  expect_true(ModelType %in% names(res))
  expect_true("Fit" %in% names(res[[ModelType]]))
})



# 3) Test bootstrap options
test_that("InputsForOutputs handles bootstrap options", {
  ModelType <- "JPS original"
  Horiz <- 100
  DesiredOutputGraphs <- c("Fit")
  OutputLabel <- "Test"
  WishStationarityQ <- 1
  DataFrequency <- "Monthly"
  WishBootstrap <- 1
  ListBoot <- list(methodBS = "bs", BlockLength = 5, ndraws = 10, pctg = c(95, 99))
  res <- InputsForOutputs(ModelType, Horiz, DesiredOutputGraphs, OutputLabel,
    WishStationarityQ, DataFrequency,
    WishBootstrap = WishBootstrap, ListBoot = ListBoot
  )
  expect_equal(res[[ModelType]]$Bootstrap$WishBoot, 1)
  expect_equal(res[[ModelType]]$Bootstrap$methodBS, "bs")
})


# 4) Test forecasting options
test_that("InputsForOutputs handles forecasting options", {
  ModelType <- "JPS original"
  Horiz <- 100
  DesiredOutputGraphs <- c("Fit")
  OutputLabel <- "Test"
  WishStationarityQ <- 1
  DataFrequency <- "Monthly"
  WishForecast <- 1
  ListForecast <- list(ForHoriz = 10, t0Sample = 1, t0Forecast = 2, ForType = "Rolling")
  res <- InputsForOutputs(ModelType, Horiz, DesiredOutputGraphs, OutputLabel,
    WishStationarityQ, DataFrequency,
    WishForecast = WishForecast, ListForecast = ListForecast
  )
  expect_equal(res[[ModelType]]$Forecasting$WishForecast, 1)
  expect_equal(res[[ModelType]]$Forecasting$ForHoriz, 10)
})
