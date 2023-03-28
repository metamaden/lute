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
#' @param k1 Vector of first type proportions. If NULL, uses 
#' seq(1e-3, 1-1e-3, 1e-3).
#' @returns lpv, a list of proportions vectors for simulation iterations.
#' @examples
#' make_lpv(k1 = c(0, 0.5, 1))
#' @export
make_lpv <- function(ktotal = 2, k1 = NULL){
  if(is(k1, "NULL")){k1 <- seq(0, 1, 1e-3)}
  num.iter <- length(k1); ki <- rev(k1)/(ktotal-1)
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
#' @param gindexv Vector of marker indices. Index values correspond to the k types,
#' and each index position represents a marker (e.g. c(1,2,2) means two markers 
#' for the second type, etc.).
#' @param num.iter Total simulation iterations.
#' @param lambda.pos Value of lambda (Poisson dist. mean) for "positive" marker 
#' status (e.g. mean of dist. for k when marker is positive for k, negative for 
#' not-k). This is passed to the argument mu when method is "nbinom".
#' @param lambda.neg Value of lambda (Poisson dist. mean) for "negative" marker 
#' status (e.g. mean of dist. for k when marker is positive for not-k, negative 
#' for k). This is passed to the argument mu when method is "nbinom".
#' @param method Type of randomization method to use. Accepts either "poisson"
#' for poisson distribution (see `?rpois` for details), or "nbinom" for the 
#' negative binomial (a.k.a. gamm poisson) distribution (see `?rnbinom` for 
#' details).
#' @param gamma.size.pos The gamma distribution magnitude for "positive" markers. 
#' This is applied when the "nbinom" method is used.
#' @param gamma.size.neg The gamma distribution magnitude for "negative" markers. 
#' This is applied when the "nbinom" method is used.
#' @param seed.num Seed value for randomization.
#' @returns Listed lgv object containing the randomized marker values across 
#' types.
#' @examples 
#' set.seed(0)
#' random_lgv(gindexv = c(rep(1, 10), rep(2, 5)))
#' @export
random_lgv <- function(gindexv, num.iter = 1, lambda.pos = 25, lambda.neg = 2, 
                       method = "nbinom", gamma.size.pos = 10, 
                       gamma.size.neg = 10, seed.num = 0){
  set.seed(seed.num); ktotal <- length(unique(gindexv))
  lgv <- lapply(seq(ktotal), function(ki){
    gmarkerv <- gindexv
    which.pos <- which(gindexv==ki); which.neg <- which(!gindexv==ki)
    if(method == "poisson"){
      gmarkerv[which.pos] <- rpois(lambda = lambda.pos, n = length(which.pos))
      gmarkerv[which.neg] <- rpois(lambda = lambda.neg, n = length(which.neg))      
    } else if(method == "nbinom"){
      gmarkerv[which.pos] <- rnbinom(size = gamma.size.pos, mu = lambda.pos,
                                     n = length(which.pos))
      gmarkerv[which.neg] <- rnbinom(size = gamma.size.neg, mu = lambda.neg,
                                     n = length(which.neg))
    } else{
      stop("Error, invalid method.")
    }

    gmarkerv
  })
  return(lapply(seq(num.iter), function(ii){lgv}))
}

#' random_sce
#'
#' Make a random SingleCellExperiment object.
#'
#' @param num.genes Number of genes to randomize.
#' @param num.cells Numnber of cells to randomize.
#' @param num.types Number of cell types to annotate.
#' @param fract.types Vector of fractions by type.
#' @param dispersion Disperison of gene expression. If NULL, uses the mean from 
#' expr.mean
#' @param expr.mean Poisson dist mean for random expression data.
#' @param seed.num Seed value for randomization of expression data.
#' @param na.include Whether to include random NA values.
#' @param na.fract Fraction of NA values to include.
#' @param zero.include Whether to include random zero-count values.
#' @param zero.fract Fraction of zero-count values to include.
#' @param verbose Whether to show verbose status messages.
#' @return New randomized SingleCellExperiment object.
#' @examples 
#' sce <- random_sce()
#' @export
random_sce <- function(num.genes = 20, num.cells = 12, num.types = 2, 
                       fract.types = NULL, dispersion = NULL, 
                       expr.mean = 10, na.include = FALSE, 
                       na.fract = 0.2, zero.include = FALSE, 
                       zero.fract = 0.2, verbose = FALSE, 
                       seed.num = 0){
  require(SingleCellExperiment)
  if(verbose){message("Getting random expression data...")}
  if(is(dispersion, "NULL")){dispersion <- expr.mean}
  mdat <- rnbinom(n = (num.cells*num.genes), 
                  size = dispersion, mu = expr.mean)
  if(na.include){ # manually add NAs
    if(verbose){message("Including NA values...")}
    num.na <- round(length(mdat)*na.fract, digits = 0)
    na.index <- sample(seq(length(mdat)), num.na)
    mdat[na.index] <- NA
  }
  if(zero.include){ # manually add zero counts
    if(verbose){message("Including NA values...")}
    num.zero <- round(length(mdat)*zero.fract, digits = 0)
    zero.index <- sample(seq(length(mdat)), num.zero)
    mdat[zero.index] <- 0
  }
  expr.ct <- matrix(mdat, ncol=num.cells, nrow=num.genes)
  if(verbose){message("Getting new colData...")}
  cellv <- paste0("cell.barcode.", seq(num.cells))
  cpertype <- round(num.cells/num.types, 0)
  
  if(is(fract.types, "NULL")){
    fract.types = rep((1/num.types), num.types)}
  typev <- paste0("type", seq(num.types))
  typev <- unlist(lapply(seq(length(typev)), function(ti){
    num <- fract.types[ti]*num.cells; rep(typev[ti], num)
  }))
  
  cd <- data.frame(cell.id = cellv, celltype = typev)
  colnames(expr.ct) <- cellv
  if(verbose){message("Getting new rowData...")}
  genev <- paste0("gene", seq(nrow(expr.ct)))
  rd <- data.frame(gene.id = genev)
  rownames(expr.ct) <- genev
  if(verbose){message("Making new sce object...")}
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts=expr.ct), colData = cd, rowData = rd)
  # manage new metadata
  description.str <- "random SingleCellExperiment made using random_sce()"
  lmd <- list(description = description.str)
  metadata(sce) <- lmd
  return(sce)
}
