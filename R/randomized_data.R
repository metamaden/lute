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
#' @param k1 Vector of first type proportions. If NULL, uses seq(1e-3, 1, 1e-3).
#' @returns lpv, a list of proportions vectors for simulation iterations.
#' @examples
#' make_lpv(k1 = c(0, 0.5, 1))
#' @export
make_lpv <- function(ktotal = 2, k1 = NULL){
  if(is(k1, "NULL")){k1 <- seq(1e-3, 1, 1e-3)}
  num.iter <- length(k1); ki <- rev(k1)/ktotal
  lpv <- lapply(seq(num.iter), function(ii){
    c(k1[ii], rep(ki[ii], ktotal-1))})
  return(lpv)
}

#' random_lgv
#'
#' Get randomized markers using Poisson distribution sampling. For a given K,
#' we assume "positive" markers have higher values than for non-K types, and 
#' thus we sample from 2 different Poisson distributions defined by different
#' lambda values (e.g. arguments lambda.pos, lambda.neg). WE also use argument 
#' gindexv to define total markers as length(gindexv) and the marker balance as
#' relative counts of each type index.
#' 
#' For example, if gindex is c(1, 1, 2), we define 3 total markers, 2 positive
#' markers for type 1 (negative for type 2) and a single positive marker for 
#' type 2 (negative for type 1).
#' 
#' @param randomize.type Type of randomization function to use. Can be "poisson"
#' for Poisson distribution sampling. 
#' @param ktotal Total types to simulate.
#' @param num.iter Total simulation iterations.
#' @param ... Additional arguments passed to randomization function. For 
#' randomize.type == "poisson", this is the mpoisson function.
#' @returns Listed lgv object containing the randomized marker values across 
#' types.
#' @examples 
#' set.seed(0)
#' random_lgv(gindexv = c(rep(1, 10), rep(2, 5)))
#' @export
random_lgv <- function(gindexv, ktotal = 2, num.iter = 1, 
                       lambda.pos = 25, lambda.neg = 2, 
                       seed.num = 0){
  set.seed(seed.num)
  lgv <- lapply(seq(ktotal), function(ki){
    gmarkerv <- gindexv
    which.pos <- which(gindexv==ki)
    which.neg <- which(!gindexv==ki)
    gmarkerv[which.pos] <- rpois(lambda = lambda.pos, n = length(which.pos))
    gmarkerv[which.neg] <- rpois(lambda = lambda.neg, n = length(which.neg))
    gmarkerv
  })
  return(lapply(seq(num.iter), function(ii){lgv}))
}
