#!/usr/bin/env R

# Functions to get randomized or synthetic datasets.
#
#
#

#' make_lpv
#'
#' Get complementary proportions for k types. The first type k1 is the vector of
#' proportions for the first type. The remaining types up to ktotal are based on
#' the reverse of k1. Types k > 1 are assumed to have equal proportions 
#' complementary to k1.
#' 
#' For k1 = c(0, 0.5, 1), ktotal = 2 will generate an additional type with 
#' proportions c(1, 0.5, 0).
#' 
#' For the same k1 above, ktotal = 3, will generate 2 types with the same 
#' proportions as c(0.5, 0.25, 0).
#'
#' @param ktotal Total types to simulate.
#' @param k1 Vector of first type proportions.
#' @returns lpv, a list of proportions vectors for simulation iterations.
#' @examples
#' make_lpv(k1 = c(0, 0.5, 1))
#' @export
make_lpv <- function(ktotal = 2, k1 = seq(0, 1, 1e-3)){
  num.iter <- length(k1); ki <- rev(k1)/ktotal
  lpv <- lapply(num.iter, function(ii){c(k1[ii], rep(ki[ii], ktotal-1))})
  return(lpv)
}