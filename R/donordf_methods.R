#!/usr/bin/env R

# Author: Sean Maden
#
# Methods for table objects of type `donordf`.
#
# Sections:
#
# 1. checks: functions that check the object classes
#
# 2. conversions: functions facilitating conversion between object classes
#

#----------
# 1. checks
#----------
# functions that check the object classes.

#' check_donordf
#'
#' Checks whether data.frame is a valid simulated donor signals data.frame by
#' evaluating the column names.
#' 
#' @param df Data.frame to check.
#' @returns boolean, TRUE if df is valid, FALSE otherwise
#' @export
check_donordf <- function(df){
  cnv <- colnames(df); lcond <- list()
  lcond[["cond.donorcol"]] <- grepl("^donor\\d", cnv)
  lcond[["cond.typecol"]] <- grepl("^type$", cnv)
  lcond[["cond.markercol"]] <- grepl("^marker$", cnv)
  lcond[["cond.markertypecol"]] <- grepl("^marker\\.type$", cnv)
  # evaluate regex
  cond.allcol <- lapply(lcond, function(ii){length(which(ii)) > 0})
  cond.allcol <- unlist(cond.allcol)
  # get final eval
  cond.output <- length(cond.allcol[cond.allcol])==4
  if(cond.output){
    return(TRUE)
  } else{
    return(FALSE)
  }
  return(NULL)
}

#---------------
# 2. conversions
#---------------
# functions facilitating conversion between object classes

#' donordf_from_mexpr
#'
#' Makes a valid donordf object from an expression matrix (e.g. rows = 
#' markers/genes, columns = samples/type data)
#'
#' @param mexpr An expression matrix (rows = markers/genes, columns = 
#' samples/type data).
#' @param verbose Whether to show verbose status messages.
#' @returns mexpr, a new expression matrix
#' @examples 
#' df <- rand_donor_marker_table()
#' madj <- donoradj_combat(df, return.type = "mexpr")
#' df.adj <- donordf_from_mexpr(mexpr = madj)
#' @export
donordf_from_mexpr <- function(mexpr, verbose = FALSE){
  typev <- unique(gsub(".*;", "", colnames(mexpr)))
  df <- do.call(rbind, lapply(seq(nrow(mexpr)), function(markeri){
    mi <- mexpr[markeri,,drop=F]
    dfi <- do.call(rbind, lapply(typev, function(typei){
      mi[,grepl(paste0(".*;", typei), colnames(mi))]
    }))
    dfi <- as.data.frame(dfi); colnames(dfi) <- gsub(";.*", "", colnames(dfi))
    dfi$type <- typev; dfi$marker <- paste0("marker",markeri)
    dfi
  }))
  # get donor summaries
  which.donor.cnv <- grepl("^donor\\d", colnames(df))
  df$donor.combn.all.mean <- rowMeans(as.matrix(df[,which.donor.cnv]))
  df$donor.combn.all.median <- rowMedians(as.matrix(df[,which.donor.cnv]))
  # get marker.type from means
  markerv <- unique(df$marker)
  mapv <- unlist(lapply(markerv, function(mi){ # get marker.type mappings
    dff <- df[df$marker == mi,]
    max.val <- max(dff[,"donor.combn.all.mean"])
    max.filt <- dff$donor.combn.all.mean==max.val
    dff[max.filt,]$type[1]
  })); names(mapv) <- markerv
  df$marker.type <- "NA"
  for(mi in markerv){df[df$marker == mi,]$marker.type <- as.character(mapv[mi])}
  # final check
  if(check_donordf(df)){
    if(verbose){message("donordf conversion success. Returning...")}
    return(df)
  } else{
    stop("Error, couldn't convert mexpr to donordf. ",
         "Do the mexpr colnames have format donor;type?")
  }
  return(NULL)
}

#' mexpr_from_donordf
#'
#' Make and expression matrix from a donordf data.frame. The donordf contains 
#' one row per marker;type, while the expression matrix contains one row per
#' marker and one column per donor;type.
#'
#' @param df A donordf type data.frame.
#' @returns mexpr, a new expression matrix (rows = genes/markers, columns = 
#' samples/types)
#' @export
mexpr_from_donordf <- function(df){
  filt.donor <- grepl("donor\\d", colnames(df)) # get donor signals
  mexpr <- do.call(rbind, lapply(unique(df$marker), function(mi){
    dff <- df[df$marker==mi, ]
    ld <- lapply(unique(dff[dff$marker==mi,]$type), function(ti){
      datv <- dff[dff$type==ti, filt.donor]
      if(nrow(datv) > 1){
        stop("Error, more than one unique marker value for each type. ",
             "Does donordf contain multiple experiment groups?")}
      names(datv) <- paste0(colnames(dff[,filt.donor]), ";", ti)
      datv
    })
    unlist(ld)
  }))
  rownames(mexpr) <- unique(df$marker)
  return(mexpr)
}

#----------------------
# 3. generators
#----------------------
# functions to generate donor.data.frame objects

#' random_donordf
#' 
#' Get a flat table of random donor marker signals by types.
#' 
#' @param ndonor Number of donors to simulate.
#' @param gindexv Vector of marker indices. Index values correspond to the k types,
#' and each index position represents a marker (e.g. c(1,2,2) means two markers 
#' for the second type, etc.).
#' @param method Randomization method passed to random_lgv(). Supports either 
#' "nbinom" for negative binomial distribution (a.k.a. gamma poisson 
#' distribution) or "poisson" for poisson distribution.
#' @param lambda.pos The mean or mu value when marker status is positive.
#' @param lambda.neg The mean or mu value when marker status is negative.
#' @param gamma.pos Magnitude of gamma dispersion value when marker status is 
#' positive.
#' @param gamma.neg Magnitude of gamma dispersion value when marker status is
#' negative.
#' @param lambda.sdoff.pos Offset SD for lambda when marker status positive.
#' @param lambda.sdoff.neg Offset SD for lambda when marker status negative.
#' @param gamma.sdoff.pos Offset SD for gamma when marker status positive.
#' @param gamma.sdoff.neg Offset SD for gamma when marker status negative.
#' @param seed.num Token to set the random seed.
#' @param vebose Whether to return verbose status messages.
#' @param ... Additional parameters passed to `random_lgv()` to get marker 
#' signals.
#' @details This function returns random donor marker signal and marker design
#' details. It supports random marker distributions drawn from either a 
#' Poisson (e.g. when method == "poisson") or a Negative Binomial/Gamma Poisson
#' (e.g. when method == "nbinom", the default) distribution. Arguments 
#' `lambda.pos` and `lambda.neg` correspond to either the means (if method is 
#' "poisson") or the mu's (if method is "nbinom") of the distributions when the 
#' marker status is positive or negative.
#' 
#' The randomization scheme is to draw a series of random offset values for 
#' lambda's and gamma's, or one offset for each simulated donor. Values are 
#' drawn from normal distributions having means of 0 and having SD's as set by
#' the arguments `lambda.sdoff.pos`, `lambda.sdoff.neg`, `gamma.sdoff.pos`, and
#' `gamma.sdoff.neg`. These options are then applied to the lambda and gamma 
#' values specified by arguments `lambda.pos`, `lambda.neg`, `gamma.pos`, and 
#' `gamma.neg`.
#' 
#' The resulting table includes one row for each marker and type. Its columns 
#' are as follows: 
#' 
#' * `donor[0-99]`: Simulated donor signals with one column for eachdonor;
#' * `donor.combn.all.mean`: Mean donor signals; 
#' * `donor.combn.all.median`: Median donor signals; 
#' * `type`: The type label (e.g. a cell type name); 
#' * `marker`: The marker label (e.g. a gene marker name); 
#' * `marker.type`: The marker's type (e.g. "type1" means this is a marker of 
#'    type1, so type1 should have high signal and remaining types should have 
#'    low signal at this marker).
#' 
#' @returns Table (data.frame) of donor marker signal and marker details.
#' @examples
#' 
#' # simulate with defaults (two donors, two marker, two types)
#' rand_donor_marker_table()
#' 
#' # simulate 10 donors, 2 types, and 3 markers (2 for type1, 1 for type2)
#' rand_donor_marker_table(ndonor = 10, gindexv = c(1,1,2))
#' 
#' # simulate 10 donors, 2 types, and 30 makers (10 for type1, 20 for type2) 
#' rand_donor_marker_table(ndonor = 10, gindexv = c(rep(1, 10), rep(2, 20)))
#' 
#' @seealso random_lgv
#' @export
random_donordf <- function(ndonor = 2, gindexv = c(1, 2), 
                                    method = "nbinom",
                                    lambda.pos = 20, lambda.neg = 2,
                                    lambda.sdoff.pos = 0, lambda.sdoff.neg = 0, 
                                    gamma.pos = 10, gamma.neg = 10,
                                    seed.num = 0, verbose = FALSE, ...){
  set.seed(seed.num)
  nmarkers <- length(gindexv); ktotal <- length(unique(gindexv))
  
  # get random hyperparameter offsets from normal dist
  # get random offsets for means
  offposv <- rnorm(n = ndonor, mean = 0, sd = lambda.sdoff.pos)
  offnegv <- rnorm(n = ndonor, mean = 0, sd = lambda.sdoff.neg)
  
  # get new hyperparameter values
  # get new means
  meanv.pos <- offposv + lambda.pos
  meanv.neg <- offnegv + lambda.neg
  
  # convert negative values
  # convert means
  meanv.pos[meanv.pos < 0] <- -1*meanv.pos[meanv.pos < 0]
  meanv.neg[meanv.neg < 0] <- -1*meanv.neg[meanv.neg < 0]
  
  # get matrix of markers (rows) by donors (cols)
  md <- do.call(cbind, lapply(seq(ndonor), function(ii){
    unlist(random_lgv(gindexv, num.iter = 1, 
                      lambda.pos = meanv.pos[ii],
                      lambda.neg = meanv.neg[ii], 
                      gamma.size.pos = gamma.pos,
                      gamma.size.neg = gamma.neg, 
                      method = method, seed.num = ii,
                      ...))
  }))
  md <- as.data.frame(md); colnames(md) <- paste0("donor", seq(ndonor))
  if(ndonor > 1){
    if(verbose){message("Getting donor summary columns...")}
    which.cnv.donor <- which(grepl("donor", colnames(md)))
    md$donor.combn.all.mean <- apply(md[,which.cnv.donor], 1, mean)
    md$donor.combn.all.median <- apply(md[,which.cnv.donor], 1, median)
  }
  md$type <- paste0("type", rep(seq(ktotal), each = nmarkers))
  md$marker <- paste0("marker", rep(seq(nmarkers), times = ktotal))
  md$marker.type <- paste0("type", gindexv)
  # final check
  
  return(md)
}

