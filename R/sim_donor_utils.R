#!/usr/bin/env R

# Author: Sean Maden
#
# Utilities to perform simulations and randomizations with one or more donors or
# cell type data sources.
#

#' get_donor_marker_flattable
#' 
#' Get a flat table of random donor marker signals by types.
#' 
#' @param ndonor Number of donors to simulate.
#' @param mean.offset.pos Poisson dist mean for randomization of offsets for
#' positive marker signals.
#' @param mean.offset.neg Poisson dist mean for randomization of offsets for
#' negative marker signals.
#' @param seed.num Token to set the random seed.
#' @param ... Additional parameters passed to `random_lgv()` to get marker 
#' signals.
#' @returns return 
#' @examples
#' 
#' get_donor_marker_flattable(ndonor = 2, gindexv = c(1,2))
#' 
#' get_donor_marker_flattable(ndonor = 10, gindexv = c(1,1,2))
#' 
#' get_donor_marker_flattable(ndonor = 10, gindexv = c(rep(1, 10), rep(2, 20)))
#' 
#' @seealso decon_results, supported_strict_methods
#' @export
new_donor_marker_table <- function(ndonor, gindexv = c(1, 2), 
                                       mean.offset.pos = 10, 
                                       mean.offset.neg = 2, 
                                       seed.num = 0, ...){
  set.seed(seed.num)
  nmarkers <- length(gindexv)
  # draw random offsets from normal dist
  offposv <- rnorm(n = ndonor, mean = mean.offset.pos)
  offnegv <- rnorm(n = ndonor, mean = mean.offset.neg)
  # get value vectors
  meanv.pos <- offposv + lambda.pos
  meanv.neg <- offnegv + lambda.pos
  # convert negative means
  meanv.pos[meanv.pos < 0] <- -1*meanv.pos
  meanv.neg[meanv.neg < 0] <- -1*meanv.neg
  # get matrix of markers (rows) by donors (cols)
  md <- do.call(cbind, lapply(seq(ndonor), function(ii){
    unlist(random_lgv(gindexv, num.iter = 1,
                      lambda.pos = meanv.pos[ii],
                      lambda.neg = meanv.neg[ii]))
  }))
  md <- as.data.frame(md)
  colnames(md) <- paste0("donor", seq(ndonor))
  md$type <- paste0("type", rep(seq(ktotal), each = nmarkers))
  md$marker <- paste0("marker", rep(seq(nmarkers), times = ktotal))
  md$marker.type <- paste0("type", gindexv)
  return(md)
}