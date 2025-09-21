library(testthat)
library(MultiATSM)

# Load inputs
data(CM_Factors)
RF_TS <- RiskFactors[1:12, ]
N <- 3

JLLinputs <- list(Economies = c("China", "Brazil"), DomUnit = "China",
                  WishSigmas = 1, SigmaNonOrtho = NULL, JLLModelType = "JLL original")


# 1) Test output structure
test_that("JLL returns correct output structure", {
  K <- nrow(RF_TS)
  res <- JLL(RF_TS, N, JLLinputs)
  expect_type(res, "list")
  expect_true(all(c("a_W", "a_DU_CS", "b", "c", "PIb", "PIac", "PI", "Ye", "k0_e", "k1_e", "k0", "k1", "Sigmas")
                  %in% names(res)))
  expect_equal(dim(res$k0), c(K, 1))
  expect_equal(dim(res$k1), c(K, K))
  expect_equal(dim(res$PI), c(K, K))
  expect_type(res$Sigmas, "list")
})


# 2) Verify the absence of a sigma list
test_that("JLL issues an error for invalid output", {
  JLLinputs$WishSigmas <- 0
  res <- JLL(RF_TS, N, JLLinputs)
  expect_false(is.list(res$Sigmas))
})


# 3) Test error for No DomUnit model
test_that("JLL issues an error for invalid output (No DomUnit model)", {
  JLLinputs$WishSigmas <- 0
  JLLinputs$JLLModelType <- "JLL No DomUnit"
  expect_error(JLL(RF_TS, N, JLLinputs, CheckInputs = TRUE))
})


# 4) Test error for no DomUnit assigned
test_that("JLL issues an error for invalid output (No DomUnit assigned)", {
  JLLinputs$WishSigmas <- 0
  JLLinputs$DomUnit <- "None"
  expect_error(JLL(RF_TS, N, JLLinputs, CheckInputs = TRUE))
})
