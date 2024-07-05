#!/usr/bin/env R

# Author: Sean Maden
#
# Functions for lute RMSE calculations.
#
#
#

#' rmseTest
#'
#' Tests the rmse function for rounding imprecision.
#'
#' @param firstVector First numeric vector.
#' @param secondVector Second numeric vector.
#' @description Takes 2 vectors of numerics
#' @returns Single numeric value
#' @examples
#' proportionsVectorPred <- seq(1e-10,2e-10,1e-11)
#' proportionsVectorTrue <- rev(proportionsVectorPred)
#' rmseTest(proportionsVectorTrue, proportionsVectorPred)
#' @details Function to test RMSE values (`./unitTests/test_rmse.R`).
#' @export
rmseTest <- function(firstVector, secondVector){
  sqrt(mean((firstVector-secondVector)^2))
}

#' rmse
#'
#' Calculates the root mean squared error (RMSE) for specified true and 
#' predicted cell type proportions.
#'
#' @param proportionsTrue cell type proportions taken as true
#' @param proportionsPred cell type proportions taken as false
#' @param summaryType Toggle summary type (either "mean" or "median")
#' @description Takes 2 vectors of numerics
#' @returns single numeric
#' @examples
#' proportionsVectorPred <- seq(1e-10,2e-10,1e-11)
#' proportionsVectorTrue <- rev(proportionsVectorPred)
#' rmse(proportionsVectorTrue, proportionsVectorPred)
#' 
#' @details
#' Function does not distinguish between true and predicted status, variable 
#' labels provided for convenience.
#' @export
rmse <- function(proportionsTrue, proportionsPred, summaryType = "mean"){
  if(summaryType == "mean"){
    sqrt(mean((proportionsTrue-proportionsPred)^2))
  } else if(summaryType == "median"){
    sqrt(median((proportionsTrue-proportionsPred)^2))
  } else{
    stop('Unrecognized summaryType')
  }
}
