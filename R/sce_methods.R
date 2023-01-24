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
#' @param scef Filtered `SingleCellExperiment` object. Should reflect the data
#' for a specific type.
#' @param group.variable Variable containing group information on scef colData.
#' @param ugroupv Vector of all unique groups to consider, of which levels in
#' group.variable entries should be identical or a subset.
#' @param summarytype Whether to summarize data on rowData (e.g. one entry per 
#' row/gene), or otherwise use the colData collapsed on means (e.g. take the 
#' mean expression across cells for each gene, returning one row per type).
#' @param groupstat Summary statistics to calculate. Can be either of "mean",
#' "median", "var", "sd", or "numzero" (i.e. number of zero-value entries).
#' @param verbose Whether to show verbose status messages.
#' @returns data.frame containing group-level summary statistics for all groups 
#' specified in ugroupv.
#' @details Computes summary statistics for either rowData (e.g. genes, markers,
#' etc.) or colData (e.g. samples, cells, etc.) for an object of class 
#' SingleCellExperiment.
#' 
#' 
#' @examples 
#' sce <- random_sce()
#' colData(sce)$donor <- c(rep("donor1", 7), rep("donor2", 3))
#' dfg <- sce_groupstat(sce, "donor", c("group1", "group2"), "counts", "colData")
#' @export
sce_groupstat <- function(scef, group.variable, ugroupv, 
                          assayname = "counts", summarytype = "rowData", 
                          groupstat = c("numzero", "var"),
                          verbose = FALSE){
  groupstat.filt <- groupstat %in% c("mean", "median", "var", "sd", "numzero")
  groupstat <- groupstat[groupstat.filt]
  if(group.variable %in% colnames(colData(scef))){
    if(verbose){
      message("Appending group-level statistics for variable: ", 
              group.variable, "...")}
    dfg <- do.call(cbind, lapply(ugroupv, function(groupi){
      # make NA matrix
      dfgi <- as.data.frame(matrix(NA, 
                                   ncol = length(groupstat)+1, 
                                   nrow = ifelse(
                                     summarytype=="colData", 1, 
                                     nrow(scef))))
      colnames(dfgi) <- c("num.entries", groupstat)
      if(groupi %in% scef[[group.variable]]){
        sceff <- scef[,scef[[group.variable]]==groupi]
        exprff <- assays(sceff)[[assayname]]
        if(summarytype == "colData"){
          exprff <- matrix(colMeans(t(exprff)), nrow = 1)
        }
        dfgi$num.entries <- ncol(exprff)
        if("mean" %in% groupstat){dfgi$mean <- rowMeans(exprff)}
        if("median" %in% groupstat){dfgi$median <- rowMedians(exprff)}
        if("var" %in% groupstat){dfgi$var <- rowVars(exprff)}
        if("sd" %in% groupstat){dfgi$sd <- rowSds(exprff)}
        if("numzero" %in% groupstat){
          dfgi$numzero <- unlist(
            apply(exprff, 1, function(ri){length(which(ri==0))}))
        }
      } else{
        if(verbose){
          message("Skipping summary for group: ", 
                  groupi, "...")}
      }
      cnstr <- ";"
      if(summarytype == "colData"){
        cnstr <- paste0(cnstr, "colData_means;")}
      colnames(dfgi) <- paste0(groupi, cnstr, colnames(dfgi))
      dfgi
    }))
    return(dfg)
  } else{
    message("Warning: variable ",group.variable," not found in sce coldata.")
  }
  return(NULL)
}
