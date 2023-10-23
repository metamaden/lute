#!/usr/bin/env R

### Author: Sean Maden
###
### Get metadata about the lute package.
###

#' luteSupportedDeconvolutionAlgorithms
#'
#' View details about supported deconvolution algorithms.
#'
#' @returns Table of supported deconvolution algorithms.
#'
#' @importFrom utils read.csv
#' 
#' @examples
#' luteSupportedDeconvolutionAlgorithms()
#' 
#' @export
luteSupportedDeconvolutionAlgorithms <- function(){
  csvName <- "lute-deconvolution_transfer-learning-table.csv"
  readPath <- file.path(
    system.file("csv", package="lute"), csvName)
  csvResult <- read.csv(readPath)
  return(csvResult)
}