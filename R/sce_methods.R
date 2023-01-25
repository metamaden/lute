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
#' @param summarytype Whether to summarize data on rowData (e.g. one entry per 
#' row/gene), or otherwise use the colData collapsed on means (e.g. take the 
#' mean expression across cells for each gene, returning one row per type).
#' @param groupstat Summary statistics to calculate. Can be either of "count", 
#' "mean", "median", "var", "sd", or "numzero" (i.e. number of zero-value 
#' entries).
#' @param type.variable Variable containing type labels.
#' @param return.tall Whether to return a tall table. If FALSE, returns the 
#' rowData/colData-compatible format with groups in column names.
#' @param verbose Whether to show verbose status messages.
#' @returns data.frame containing group-level summary statistics for all groups 
#' specified in ugroupv.
#' @details Computes summary statistics for either rowData (e.g. genes, markers,
#' etc.) or colData (e.g. samples, cells, etc.) for an object of class 
#' SingleCellExperiment.
#' 
#' Dimensions in the returned object are determined by the argument 
#' `return.tall`. When this is TRUE, rows in the returned table correspond to
#' donors or donor-types (when `type.variable` is provided) if 
#' `summarytype`=="colData", and donor-genes otherwise. When this is FALSE, a 
#' colData/rowData-compatible object is returned as specified by `summarytype`.
#'  
#' @examples 
#' 
#' # make random data
#' sce = random_sce(zero.include = TRUE, zero.fract = 0.5)
#' colData(sce)$donor <- c(rep("donor1", 7), rep("donor2", 3))
#'
#' # get group summaries across types
#' sce_groupstat(sce, summarytype = "colData", return.tall = TRUE)
#'
#' # get group-gene summaries across types
#' sce_groupstat(sce, summarytype = "rowData", return.tall = TRUE)
#'
#' # get tall table of group-type summaries
#' sce_groupstat(sce, type.variable = "celltype", 
#'              summarytype = "colData", return.tall = TRUE)
#'
#' # get coldata-compatible group-type summaries
#' sce_groupstat(sce, type.variable = "celltype", 
#'              summarytype = "colData", return.tall = FALSE)
#'
#' # get rowdata-compatible group-type summaries
#' sce_groupstat(sce, type.variable = "celltype", 
#'              summarytype = "rowData", return.tall = FALSE)
#' 
#' @export
sce_groupstat <- function(sce, group.variable = "donor",  
                          assayname = "counts", summarytype = "colData", 
                          groupstat = c("count", "numzero", "var"),
                          return.tall = FALSE, type.variable = NULL, 
                          verbose = FALSE){
  # run checks
  # check sce
  if(!is(sce, "SingleCellExperiment")){
    stop("Error, sce must be a SingleCellExperiment.")}
  # check assay signals
  expr <- assays(sce)[[assayname]]
  cond <- is(expr, "matrix")|is(expr, "DelayedArray")
  if(!cond){
    stop("Error, assay signals should be either a matrix or DelayedArray.")}
  # filter group stats
  groupstat.filt <- groupstat %in% c("count", "mean", "median", 
                                     "var", "sd", "numzero")
  groupstat <- groupstat[groupstat.filt]
  if(verbose){message("Checking colData variables...")}; cd <- colData(sce)
  lvar <- check_coldata(cd = cd, var = c(group.variable, type.variable))
  # get group stats df
  ugroupv <- lvar[[group.variable]]$uvec
  groupv <- lvar[[group.variable]]$vec
  if(!is(type.variable, "NULL")){
    utypev <- lvar[[type.variable]]$uvec
    typev <- lvar[[type.variable]]$vec
  }
  
  ldf <- lapply(ugroupv, function(groupi){
    if(is(type.variable, "NULL")){
      exprf <- assays(sce[,groupv==groupi])[[assayname]]
      dfi <- get_groupstat_df(exprf, groupstat = groupstat, 
                              summarytype = summarytype)
      
    } else{
      dfi <- do.call(rbind, lapply(utypev, function(typei){
        sce.filt <- groupv==groupi & typev==typei
        exprff <- assays(sce[,sce.filt])[[assayname]]
        dfii <- get_groupstat_df(exprff, groupstat = groupstat, 
                                 summarytype = summarytype)
        dfii$type <- typei; return(dfii)
      }))
    }
    dfi$group <- groupi
    dfi
  })
  if(return.tall){
    dfr <- do.call(rbind, lapply(ldf, function(dfi){dfi}))
  } else{
    dfr <- do.call(cbind, lapply(ldf, function(dfi){
      groupi <- unique(dfi$group)[1]
      if(summarytype == "rowData" & !is(type.variable, "NULL")){
        dfi <- do.call(cbind, lapply(utypev, function(typei){
          dfii <- dfi[dfi$type == typei,,drop = F]
          dfii <- dfii[,!colnames(dfii) %in% c("group", "marker", "type")]
          colnames(dfii) <- paste0(groupi, ";", typei, ";", colnames(dfii))
          dfii
        }))
      } else{
        dfi <- dfi[,!colnames(dfi) %in% c("group", "type", "marker")]
        colnames(dfi) <- paste0(groupi, ";", colnames(dfi))
      }
      dfi
    }))
    cond <- summarytype == 'colData' & !is(type.variable, "NULL")
    if(cond){dfr$type <- utypev}
  }
  return(dfr)
}

#' check_coldata
#'
#' Check coldata columns.
#' 
#' @param cd Matrix of colData.
#' @param varv Vector of variable names to check
#' @returns list, contains vector of variable values and vector of unique levels
#' for variables found in cd.
#' @export
check_coldata <- function(cd, varv){
  varv <- varv[!is.null(varv)]
  lr <- lapply(varv, function(vari){
    if(!vari %in% colnames(cd)){
      stop("Error, didn't find variable ", vari, "in colData colnames")}
    lri <- list()
    lri[["vec"]] <- cd[,vari]
    lri[["uvec"]] <- unique(cd[,vari])
    return(lri)
  })
  names(lr) <- varv
  return(lr)
}

#' get_groupstat_df
#'
#' Get grouped statistics from an expression matrix.
#' 
#' @param exprf Filtered expression matrix.
#' @param groupstat Vector of valid group statistics.
#' @param summarytype Type of summary. Either colData or rowData.
#' @param round.digits Digits to round.
#' @returns Table of type data.frame containing summary statistics.
#' @examples
#' sce = random_sce(zero.include = TRUE, zero.fract = 0.5)
#' colData(sce)$donor <- c(rep("donor1", 7), rep("donor2", 3))
#' exprf <- assays(sce)$counts
#' get_groupstat_df(exprf)
#' 
#' @export
get_groupstat_df <- function(exprf, groupstat = c("count", "var"), 
                             summarytype = "colData", round.digits = 3){
  # parse summarytype
  ngene <- nrow(exprf); ncell <- ncol(exprf); exprff <- exprf
  if(summarytype == "colData"){exprff <- t(as.matrix(colMeans(exprf)))}
  # make na matrix
  groupstatf <- groupstat[!grepl("^count.*|^numzero.*", groupstat)]
  mna <- matrix(NA, ncol = length(groupstatf), nrow = nrow(exprff))
  dfti <- as.data.frame(mna); colnames(dfti) <- groupstatf
  cond <- length(which(grepl("^count.*", groupstat)))
  if(cond > 0){
    dfti$count.cells <- ncell; dfti$count.genes <- ngene
  }
  for(ri in seq(nrow(exprff))){
    if(!is(exprff, "NULL")){
      # get group stats
      ei <- exprff[ri,,drop = F]
      if("mean" %in% groupstat){
        meani <- rowMeans(ei)
        dfti[ri,]$mean <- round(meani, digits = round.digits)
      }
      if("median" %in% groupstat){
        mediani <- rowMedians(ei)
        dfti[ri,]$median <- round(mediani, digits = round.digits)
      }
      if("var" %in% groupstat){
        vari <- rowVars(ei)
        dfti[ri,]$var <- round(vari, digits = round.digits)
      }
      if("sd" %in% groupstat){
        sdi <- rowSds(ei)
        dfti[ri,]$sd <- round(sdi, digits = round.digits)
      }
    }
  }
  # parse zero counts
  cond <- length(which(grepl("^numzero.*", groupstat)))
  if(cond > 0){
    if(summarytype == "colData"){
      # get summaries across cells
      zct.cell <- median(apply(exprf, 1, function(ri){length(which(ri==0))}))
      zfr.cell <- median(zct.cell/ncell)
      dfti$median.count.cells.zero <- round(zct.cell, digits = round.digits)
      dfti$median.fract.cells.zero <- round(zfr.cell, digits = round.digits)
      # get summaries across genes
      zct.gene <- median(apply(exprf, 2, function(ci){length(which(ci==0))}))
      zfr.gene <- median(zct.gene/ngene)
      dfti$median.count.genes.zero <- round(zct.gene, digits = round.digits)
      dfti$median.fract.genes.zero <- round(zfr.gene, digits = round.digits)
    } else{
      dfti$count.cells.zero <- apply(exprf, 1, function(ri){
        length(which(ri == 0))})
    }
  }
  if(summarytype == "rowData"){dfti$marker <- rownames(exprf)}
  return(dfti)
}
