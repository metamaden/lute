#!/usr/bin/env R

# Author: Sean Maden
#
# Test RMSE function behavior.
#
#
#

test_that(
  "Difference equal to expectation", {
  differenceExpect <- 1e-5
  value1 <- 1e-4
  value2 <- 9e-5
  value1-value2==differenceExpect
})

test_that(
  "Square equal to expectation", {
    numDigits <- 10
    squareRootExpect <- 1e-10
    value1 <- 1e-4
    value2 <- 9e-5
    expect_equal(
      format((value1-value2)^2, scientific = TRUE, digits = numDigits),
      format(squareRootExpect, scientific = TRUE, digits = numDigits))
  })

test_that(
  "Mean squared error equal to expectation", {
    numDigits <- 10
    meanSquaredErrorExpect <- 1e-10
    value1 <- 1e-4
    value2 <- 9e-5
    expect_equal(
      format(mean((value1-value2)^2), scientific = TRUE, digits = numDigits),
      format((value1-value2)^2, scientific = TRUE, digits = numDigits))
  })

test_that(
  "Root mean squared error equal to expectation", {
    numDigits <- 10
    meanSquaredErrorExpect <- 1e-5
    value1 <- 1e-4
    value2 <- 9e-5
    expect_equal(
      format(sqrt(mean((value1-value2)^2)), scientific = TRUE, digits = numDigits),
      format(meanSquaredErrorExpect, scientific = TRUE, digits = numDigits))
  })

test_that(
  "Root mean squared error equal to expectation (vector)", {
    numDigits <- 10
    meanSquaredErrorExpect <- 1e-5
    value1 <- c(1e-4, 1e-4)
    value2 <- c(9e-5, 9e-5)
    expect_equal(
      format(sqrt(mean((value1-value2)^2)), scientific = TRUE, digits = numDigits),
      format(meanSquaredErrorExpect, scientific = TRUE, digits = numDigits))
  })

test_that(
  "Root mean squared error (lute::rmse) equal to expectation (vector)", {
    numDigits <- 10
    meanSquaredErrorExpect <- 1e-5
    value1 <- seq(1e-10,2e-10,1e-11)
    value2 <- rev(value1)
    expect_equal(
      format(sqrt(mean((value1-value2)^2)), 
             scientific = TRUE, digits = numDigits),
      format(rmseTest(value1, value2), 
             scientific = TRUE, digits = numDigits))
    
    
    expect_true(
      tryCatch(
        rmseTest(as.character(value1), value2),
        error = function(e){TRUE},
        finally = FALSE
      ) 
    )
    
    
  })

test_that(
  "Root mean squared error (lute::rmse) equal to expectation (vector)", {
    numDigits <- 10
    value1 <- seq(1e-10,2e-10,1e-11)
    value2 <- rev(value1)
    expect_equal(
      format(sqrt(mean((value1-value2)^2)), 
             scientific = TRUE, digits = numDigits),
      format(rmse(value1, value2, "mean"), 
             scientific = TRUE, digits = numDigits))
    expect_equal(
      format(sqrt(median((value1-value2)^2)), 
             scientific = TRUE, digits = numDigits),
      format(rmse(value1, value2, "median"), 
             scientific = TRUE, digits = numDigits))
    expect_equal(
      round(sqrt(median((value1-value2)^2)), 10),
      round(rmse(value1, value2, "median"), 10))
    expect_equal(
      round(sqrt(mean((value1-value2)^2)), 10),
      round(rmse(value1, value2, "mean"), 10))
})

