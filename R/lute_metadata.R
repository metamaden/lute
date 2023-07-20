#!/usr/bin/env R

# Author: Sean Maden
#
# Get metadata about the lute package.
#

#' lute_supported_deconvolution_algorithms
#'
#' View details about supported deconvolution algorithms.
#'
#' @returns Table of supported deconvolution algorithms.
#'
#' @importFrom utils read.csv
#' 
#' @export
lute_supported_deconvolution_algorithms <- function(){
  csv.name <- "lute-deconvolution_transfer-learning-table.csv"
  path <- system.file("csv", package = "lute")
  path <- file.path(path, csv.name)
  csv <- read.csv(path)
  return(csv)
}