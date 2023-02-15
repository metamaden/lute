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
#' @returns Check outcome as boolean; TRUE if `df` passes, or FALSE otherwise.
#' @export
check_donordf <- function(df){
  # check class
  if(!is(df, "donor.data.frame")){
    stop("Error, df is not an object of class `donor.data.frame`.")
  }
  # check type variable
  if(!is(df$type, "factor")){
    message("Warning, type variable is not an ordered factor.")
  }
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
#' samples/type data). Column names correspond to the types used.
#' @param verbose Whether to show verbose status messages.
#' @returns New table of type `donor.data.frame`.
#' @examples 
#' sce = random_sce()
#' groupv.cd = c(rep("donor1", 7), rep("donor2", 3))
#' typev.cd = sce[["celltype"]]
#' mexpr <- assays(sce)$counts
#' df <- donordf_from_mexpr(mexpr = mexpr, groupv.cd = groupv.cd, 
#' typev.cd = typev.cd)
#' 
#' @export
donordf_from_mexpr <- function(mexpr, groupv.cd, typev.cd, typev.rd = NULL, 
                               verbose = FALSE){
  utypev.cd <- unique(typev.cd)
  ugroupv.cd <- unique(groupv.cd)
  
  if(is(typev.rd, "NULL")){ 
    # get marker.type var
    mgk <- do.call(cbind, lapply(utypev.cd, function(ti){
      matrix(rowMeans(mexpr[,typev.cd==ti]), ncol = 1)
    }))
    typev.rd <- unlist(apply(mgk, 1, function(ri){
      utypev.cd[ri==max(ri)]}))
  }
  
  df <- do.call(cbind, lapply(ugroupv.cd, function(di){
    filt.donor <- groupv.cd==di
    do.call(rbind, lapply(utypev.cd, function(ti){
      filt.type <- filt.donor & typev.cd==ti
      mdat <- rowMeans(mexpr[,filt.type], na.rm = T)
      matrix(mdat, ncol = 1)
    }))
  }))
  colnames(df) <- as.character(ugroupv.cd)
  donor.mean <- rowMeans(df, na.rm = T)
  donor.median <- rowMedians(df, na.rm = T)
  df <- as.data.frame(df)
  df$donor.mean <- donor.mean
  df$donor.median <- donor.median
  df$type <- rep(utypev.cd, each = nrow(mexpr))
  df$type <- factor(df$type, levels = unique(df$type)[order(unique(df$type))])
  df$marker <- rep(rownames(mexpr), length(utypev.cd))
  df$marker.type <- rep(typev.rd, length(utypev.cd))
  df <- as(df, "donor.data.frame")
  if(!check_donordf(df)){stop("Error, couldn't make new valid donor.data.frame")}
  return(df)
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
  if(!check_donordf(df)){stop("Error, table isn't a valid donor.data.frame.")}
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

#' donordf_from_df
#'
#' Make a new table of type `donor.data.frame` from a `data.frame`.
#' 
#' @param df Table of type `data.frame`.
#' @returns Table of type `donor.data.frame`.
#' @details Conversion function for `donor.data.frame` class. If missing columns
#' aren't present in the provided `data.frame`, add them and set their values to
#' NA.
#' @examples
#' df <- data.frame(donor1 = sample(5), marker = rep("marker1", 5))
#' donordf_from_df(df)
#' 
#' @export
donordf_from_df <- function(df){
  if(!is(df, "data.frame")){stop("Error, table must be a data.frame.")}
  cnv <- colnames(df); filt.donor <- grepl("^donor\\d", cnv)
  if(length(which(filt.donor))==0){df$donor1 <- NA}
  if(!"marker" %in% cnv){df$marker <- df$marker.type <- NA}
  if(!"type" %in% cnv){df$type <- df$marker.type <- NA}
  donordf <- as(df, "donor.data.frame")
  if(check_donordf(donordf)){ 
    return(donordf)
  } else{
    stop("Error, couldn't make new donor.data.frame from provided tables.")
  }
  return(NULL)
}


#--------------
# 3. generators
#--------------
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
#' @returns New table of type `donor.data.frame` containing donor marker signal 
#' and marker label details.
#' @examples
#' # simulate with defaults (two donors, two marker, two types)
#' random_donordf()
#' 
#' # simulate 10 donors, 2 types, and 3 markers (2 for type1, 1 for type2)
#' rand_donor_marker_table(ndonor = 10, gindexv = c(1,1,2))
#' 
#' # simulate 10 donors, 2 types, and 30 makers (10 for type1, 20 for type2) 
#' rand_donor_marker_table(ndonor = 10, gindexv = c(rep(1, 10), rep(2, 20)))
#' 
#' @seealso random_lgv, biasexpt, donor_marker_biasexpt
#' @export
random_donordf <- function(ndonor = 2, gindexv = c(1, 2), method = "nbinom",
                           lambda.pos = 20, lambda.neg = 2,lambda.sdoff.pos = 0, 
                           lambda.sdoff.neg = 0, gamma.pos = 10, gamma.neg = 10,
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
    md$donor.mean <- apply(md[,which.cnv.donor], 1, mean)
    md$donor.median <- apply(md[,which.cnv.donor], 1, median)
  }
  # get type label variable
  type.vector <- paste0("type", rep(seq(ktotal), each = nmarkers))
  type.levels <- unique(type.vector)
  type.levels <- marker.levels[order(type.levels)]
  type.factor <- factor(type.vector, levels = type.levels)
  md$type <- type.factor
  
  # get marker label variable
  marker.vector <- paste0("marker", rep(seq(nmarkers), times = ktotal))
  marker.levels <- unique(marker.vector)
  marker.levels <- marker.levels[order(marker.levels)]
  marker.factor <- factor(marker.vector, levels = marker.levels)
  md$marker <- marker.factor
  
  # get marker type
  marker.type.vector <- paste0("type", gindexv)
  marker.type.levels <- unique(marker.type.vector)
  marker.type.levels <- marker.type.levels[order(marker.type.levels)]
  marker.type.factor <- factor(marker.type.vector, levels = marker.type.levels)
  md$marker.type <- marker.type.factor
  
  # convert to donor.data.frame
  md <- as(md, "donor.data.frame")
  if(check_donordf(md)){ # final check
    return(md)
  } else{
    stop("Error, couldn't make new donor.data.frame object.")
  }
  return(NULL)
}

