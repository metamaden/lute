#!/usr/bin/env R

# Author: Sean Maden
#
# Methods for SingleCellExperiment objects.
#

#-------------------
# summary statistics
#-------------------
#' sce_groupstat
#' 
#' Get group-level summary statistics, either on rowData or collapsed colData.
#' 
#' @param sce Filtered `SingleCellExperiment` object. Should reflect the data
#' for a specific type.
#' @param group.variable Variable containing group information on scef colData.
#' @param ugroupv Vector of all unique groups to consider. If NULL, considers
#' all available group levels in group.variable.
#' @param summarytype Whether to summarize data on rowData (e.g. one entry per 
#' row/gene), or otherwise use the colData collapsed on means (e.g. take the 
#' mean expression across cells for each gene, returning one row per type).
#' @param groupstat Summary statistics to calculate. Can be either of "count", 
#' "mean", "median", "var", "sd", or "numzero" (i.e. number of zero-value 
#' entries).
#' @param verbose Whether to show verbose status messages.
#' @returns data.frame containing group-level summary statistics for all groups 
#' specified in ugroupv.
#' @details Computes summary statistics for either rowData (e.g. genes, markers,
#' etc.) or colData (e.g. samples, cells, etc.) for an object of class 
#' SingleCellExperiment.
#' 
#' @examples 
#' sce <- random_sce()
#' colData(sce)$donor <- c(rep("donor1", 7), rep("donor2", 3))
#' dfg <- sce_groupstat(sce, "donor", c("group1", "group2"), "counts", "colData")
#' @export
sce_groupstat <- function(sce, group.variable, ugroupv = NULL, 
                          assayname = "counts", summarytype = "colData", 
                          groupstat = c("count", "numzero", "var"),
                          verbose = FALSE){
  if(!is(sce, "SingleCellExperiment")){
    stop("Error, sce must be a SingleCellExperiment.")}
  cd <- colData(sce)
  # function
  groupstat.filt <- groupstat %in% c("count", "mean", "median", 
                                     "var", "sd", "numzero")
  groupstat <- groupstat[groupstat.filt]
  # parse group variable
  if(group.variable %in% colnames(cd)){
    groupv <- cd[, group.variable]
    if(is(ugroupv, "NULL")){
      ugroupv <- unique(groupv)
    } else{
      ugroupv <- ugroupv[ugroupv %in% groupv]
      if(length(ugroupv)==0){
        stop("Error, no labels in ugroupv found in sce colData")}
    }
  } else{
    stop("Error, didn't find group.variable in sce data.")
  }
  # parse type variable
  if(type.variable %in% colnames(cd)){
    typev <- cd[, type.variable]
    if(is(utypev, "NULL")){
      utypev <- unique(typev)
    } else{
      utypev <- utypev[utypev %in% typev]
      if(length(utypev)==0){
        stop("Error, no labels in ugroupv found in sce colData")}
    }
  } else{
    stop("Error, didn't find group.variable in sce data.")
  }
  # get group stats df
  filtvar <- ugroupv; filtv <- groupv
  if(!is(type.variable, "NULL")){
    ufiltvar <- paste0(rep(ugroupv, each = length(utypev))
                       , ";", utypev)
    filtv <- paste0(groupv, ";", typev)
  }
  dfg <- do.call(cbind, lapply(ugroupv, function(groupi){
    ufiltvari <- ufiltv[grepl(paste0(groupi, ";"), ufiltv)]
    dfgi <- do.call(rbind, lapply(ufiltvari, function(filti){
      scef <- sce[, filtv==filti]
      if(ncol(scef) > 0){
        exprf <- assays(scef)[[assayname]]
        if(summarytype == "colData"){
          exprf <- matrix(colMeans(exprf), nrow = 1)
        } else{exprf <- t(exprf)}
      } else{
        exprf <- NULL
      }
      get_groupstat_df(exprf, groupstat = groupstat, 
                       summarytype = summarytype)
    }))
    colnames(dfgi) <- paste0(groupi, ";", colnames(dfgi)); dfgi
  }))
  if(!is(type.variable, "NULL")){dfg$type <- utypev}
  return(dfg)
}

#' get_groupstat_df
#'
#' Get grouped statistics from an expression matrix.
#' 
#' @param exprf Filtered expression matrix.
#' @param groupstat Vector of valid group statistics.
#' @param summarytype Type of summary. Either colData or rowData.
#' @returns Table of type data.frame containing summary statistics.
#' @example
#' mexpr <- matrix(sample(10), nrow = 1, ncol = 10)
#' get_groupstat_df(mexpr, c("count", "var"), "colData")
#' @export
get_groupstat_df <- function(exprf, groupstat, summarytype){
  # make na matrix
  mna <- matrix(NA, ncol = length(groupstat), 
                nrow = ifelse(summarytype=="colData", 1, nrow(sce)))
  dfti <- as.data.frame(mna); colnames(dfti) <- groupstat
  if(!is(exprf, "NULL")){
    # get group stats
    if("count" %in% groupstat){
      dfti$count <- ifelse(summarytype == "colData", ncol(exprf), nrow(exprf))
    }
    if("mean" %in% groupstat){dfti$mean <- rowMeans(exprf)}
    if("median" %in% groupstat){dfti$median <- rowMedians(exprf)}
    if("var" %in% groupstat){dfti$var <- rowVars(exprf)}
    if("sd" %in% groupstat){dfti$sd <- rowSds(exprf)}
    if("numzero" %in% groupstat){
      dfti$numzero <- unlist(
        apply(exprf, 1, function(ri){length(which(ri==0))}))
    }
  }
  return(dfti)
}
