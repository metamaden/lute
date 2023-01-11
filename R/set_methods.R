#!/usr/bin/env R

# Author: Sean Maden
#
# Methods for SummarizedExperimentType objects. 
#
# Sections:
#
#   1. define class functions: Class definitions and methods to make new 
#       objects.
#
#   2. methods for converting between classes: Between-object class conversions.
#
#   3. methods for adjustments: Adjustment and weighting methods.
#
#   4. set simulations: Methods supporting simultions using 
#       SummarizedExperimentTypes objects.
#

#--------------------------
# 1. define class functions
#--------------------------
#' sce_groupstat
#' 
#' Get group-level summary statistics, either on rowData or collapsed colData.
#' 
#' @param scef Filtered `SingleCellExperiment` object. Should reflect the data
#' for a specific type.
#' @param groupvar Variable containing group information on scef colData.
#' @param ugroupv Vector of all unique groups to consider, of which levels in
#' groupvar entries should be identical or a subset.
#' @param summarytype Whether to summarize data on rowData (e.g. one entry per 
#' row/gene), or otherwise use the colData collapsed on means (e.g. take the 
#' mean expression across cells for each gene, returning one row per type).
#' @param groupstat Summary statistics to calculate. Can be either of "mean",
#' "median", "var", "sd", or "numzero" (i.e. number of zero-value entries).
#' @param verbose Whether to show verbose status messages.
#' @returns data.frame containing group-level summary statistics for all groups 
#' specified in ugroupv.
#' @examples 
#' sce <- random_sce()
#' colData(sce)$donor <- c(rep("donor1", 7), rep("donor2", 3))
#' dfg <- sce_groupstat(sce, "donor", c("group1", "group2"), "counts", "colData")
#' @export
sce_groupstat <- function(scef, groupvar, ugroupv, 
                          assayname = "counts", summarytype = "rowData", 
                          groupstat = c("numzero", "var"),
                          verbose = FALSE){
  groupstat.filt <- groupstat %in% c("mean", "median", "var", "sd", "numzero")
  groupstat <- groupstat[groupstat.filt]
  if(groupvar %in% colnames(colData(scef))){
    if(verbose){
      message("Appending group-level statistics for variable: ", 
              groupvar, "...")}
    dfg <- do.call(cbind, lapply(ugroupv, function(groupi){
      # make NA matrix
      dfgi <- as.data.frame(matrix(NA, 
                                   ncol = length(groupstat)+1, 
                                   nrow = ifelse(
                                     summarytype=="colData", 1, 
                                     nrow(scef))))
      colnames(dfgi) <- c("num.entries", groupstat)
      if(groupi %in% scef[[groupvar]]){
        sceff <- scef[,scef[[groupvar]]==groupi]
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
    message("Warning: variable ",groupvar," not found in sce coldata.")
  }
  return(NULL)
}

#------------------------------------------
# 2. methods for converting between classes
#------------------------------------------

#' set_from_sce
#'
#' Make an object of class `SummarizedExperimentTypes` from an object of class
#' `SingleCellExperiment` or `SummarizedExperiment`.
#'
#' @param sce `SingleCellExperiment` object.
#' @param method Method to summarize assays on types.
#' @param typevar Variable containing the type-level labels.
#' @param assayname Name of assays data in assays(sce).
#' @param groupvar Variable containing group-level labels used for group-wise
#' summaries. If NULL, skip group summaries.
#' @param verbose Whether to show verbose status messages.
#' @param ... Additional arguments specified for `sce_groupstat()` summaries by 
#' groups (see `?sce_groupstat` for details).
#' @returns `SummarizedExperimentTypes` object.
#' @examples
#' sce <- random_sce()
#' set <- set_from_sce(sce)
#' @export
#'
set_from_sce <- function(sce, groupvar = NULL, method = "mean", 
                         typevar = "celltype", assayname = "counts", 
                         verbose = FALSE, ...){
  if(!(is(sce, "SingleCellExperiment")|is(sce, "SummarizedExperiment"))){
    stop("sce must be of class SingleCellExperiment or SummarizedExperiment.")}
  typev <- unique(sce[[typevar]])
  expr.sce <- assays(sce)$counts
  if(!is(groupvar, "NULL")){
    ugroupv <- unique(sce[[groupvar]]) # get all possible group lvls
  }
  expr.set <- do.call(cbind, lapply(typev, function(typei){
    if(verbose){message("Summarizing type: ", typei, "...")}
    scef <- sce[,sce[[typevar]]==typei]
    exprf <- assays(scef)[[assayname]]
    gene.varv <- apply(exprf,1,var)
    gene.sdv <- apply(exprf,1,sd)
    gene.max <- apply(exprf,1,max)
    gene.min <- apply(exprf,1,min)
    if(method == "mean"){
      exprnew <- matrix(rowMeans(exprf), ncol = 1)
    } else if(method == "median"){
      exprnew <- matrix(rowMedians(exprf), ncol = 1)
    } else{
      exprnew <- exprf[,1,drop=F]
    }
    colnames(exprnew) <- paste0(typei, ";expr")
    de <- data.frame(var = gene.varv, sdv = gene.sdv, max = gene.max,
                     min = gene.min)
    # parse group-level statistics
    if(!is(groupvar, "NULL")){
      dfg <- sce_groupstat(scef = scef, groupvar = groupvar, ugroupv = ugroupv,
                           summarytype = "rowData", assayname = assayname, 
                           verbose = verbose, ...)
      condv <- is(dfg, "data.frame") & nrow(dfg) == nrow(de)
      if(condv){if(verbose){message("Binding group-level data.")}
        de <- cbind(de, dfg)}
    }
    colnames(de) <- paste0(typei, ";", colnames(de))
    return(cbind(exprnew, de))
  }))
  which.mexpr <- grepl(".*;expr$", colnames(expr.set))
  mexpr <- expr.set[,which.mexpr] # expr data
  colnames(mexpr) <- gsub(";expr.*", "", colnames(mexpr))
  rd <- expr.set[,!which.mexpr] # rowdata
  
  # get coldata
  cd.sce <- colData(sce)
  cd <- do.call(rbind, lapply(typev, function(typei){
    if(verbose){message("Working on type: ", typei, "...")}
    scef <- sce[,sce[[typevar]]==typei]
    num.cells <- ncol(scef)
    # parse 0-expr genes/cells
    count.zeroexpr <- apply(assays(scef)$counts, 1, 
                            function(ri){length(ri[ri==0])})
    num.allzeroexpr <- length(which(count.zeroexpr==num.cells))
    mean.zerocount <- mean(count.zeroexpr)
    median.zerocount <- median(count.zeroexpr)
    var.zerocount <- var(count.zeroexpr)
    sd.zerocount <- sd(count.zeroexpr)
    # get return df
    dfr <- data.frame(type = typei,
                      num.cells = num.cells,
                      num.allzeroexpr = num.allzeroexpr,
                      mean.zerocount = mean.zerocount,
                      median.zerocount = median.zerocount,
                      var.zerocount = var.zerocount,
                      sd.zerocount = sd.zerocount)
    # parse group-level statistics
    if(!is(groupvar, "NULL")){
      dfg <- sce_groupstat(scef = scef, groupvar = groupvar, ugroupv = ugroupv, 
                           assayname = assayname, summarytype = "colData", 
                           verbose = verbose, ...)
      condv <- is(dfg, "data.frame") & nrow(dfg) == nrow(dfr)
      if(condv){
        if(verbose){message("Binding group-level data.")};dfr <- cbind(dfr, dfg)
      }
    }
    return(dfr)
  }))
  # make new set object
  lassays <- list(mexpr); names(lassays) <- paste0(assayname, "_bytype")
  new.set.md <- list(assay.info = list(typevar = typevar))
  set <- SummarizedExperimentTypes(assays = lassays, rowData = rd, colData = cd)
  return(set)
}

#' set_from_set
#'
#' Expects a set object with type;group nesting. Returns a set object summarized
#' by each type across any available groups.
#'
#' @param 
#' @param 
#' 
#' @examples 
#' sce <- random_sce()
#' sce[["donor"]] <- c(rep("donor1", 2), rep("donor2", 8))
#' sce[["typevar"]] <- paste0(sce[["celltype"]], ";", sce[["donor"]])
#' set <- set_from_sce(sce, groupvar = "donor", typevar = "typevar")
#'
set_from_set <- function(set, groupvar = "donor", typevar = "celltype", 
                         method = "mean", assayname = "counts_bytype",  
                         verbose = FALSE, ...){
  # check set class
  setclass.cond <- (is(set, "SummarizedExperimentTypes")|
                      is(set, "RangedSummarizedExperiment Types"))
  if(!setclass.cond){
    stop("set must be of class SummarizedExperimentTypes or similar.")}
  
  # get params
  typev <- unique(sce[[typevar]])
  expr.sce <- assays(sce)$counts
  if(!is(groupvar, "NULL")){
    ugroupv <- unique(sce[[groupvar]]) # get all possible group lvls
  }
  
  # make new assays data
  cnv <- colnames(set)
  ldat <- lapply(typev, function(typei){
    type.filt <- grepl(paste0("^", typei, "(;.*|$)"), cnv)
    mei <- assays(set[,type.filt])[[assayname]]
    mei <- as.matrix(mei)
    unique.group.lvl <- unique(gsub(paste0("^",typei,";"), "", colnames(mei)))
    unique.group.lvl <- paste0(unique.group.lvl, collapse = ";")
    num.group <- ncol(mei)
    group.var <- rowVars(mei)
    group.mean <- rowMeans(mei)
    list(rd = data.frame(group.var = group.var,
                         group.mean = group.mean),
         cd = data.frame(type = typei, num.group = num.group, 
                         group.lvl = unique.group.lvl),
         mexpr = rowMeans(mei))
  })
  names(ldat) <- typev
  
  # assemble set components
  # assay expr table
  new.assay <- do.call(cbind, lapply(ldat, function(li){li$mexpr}))
  colnames(new.assay) <- names(ldat)
  # coldata
  new.cd <- do.call(rbind, lapply(ldat, function(li){li$cd}))
  rownames(new.cd) <- colnames(new.assay)
  # rowdata
  new.rd <- do.call(cbind, lapply(seq(length(ldat)), function(ii){
    typei <- names(ldat)[ii]; rdi <- ldat[[ii]]$rd
    colnames(rdi) <- paste0(typei, ";", colnames(rdi))
    rdi
  }))
  # metadata
  new.md <- list(assay.info = list(type = typevar, groupvar = groupvar),
                 set.original = set)
  
  # make new set object
  set <- SummarizedExperimentTypes(assays = list(assayname = new.assay), 
                                   rowData = rd, colData = cd,
                                   metadata = new.md)
  return(set)
}


#' convert_sce
#'
#' Manage conversions to a `SummarizedExperimentTypes` object.
#'
#' @param
#' @param 
#' 
#'

#---------------------------
# 3. methods for adjustments
#---------------------------

#' append_groupvar_adj
#'
#' Perform group-wise adjustments on a SummarizedExperimentTypes object. Appends 
#' group adjustments to rowData and group-adjusted signals to assays.
#'
#' @param set A SummarizedExperimentTypes object
#' @param groupvar_pattern String pattern to identify group variables from 
#' rowData, if present.
#' @param assayname Name of original assays data to adjust, as identifiable from
#' names(assays(set)).
#' @param rd.method.str Character string to append for new assays and rowData 
#' matrices.
#' @param type Adjustment type (options: "mgvdenom").
#' @param verbose Whether to show verbose status messages.
#' @param ... Additional arguments passed for specific adjustment type function.
#' @return set object with adjusted marker signal assay data appended.
#' @examples
#' sce <- random_sce()
#' sce[["donor"]] <- c(rep("donor1", 2), rep("donor2", 8))
#' set <- set_from_sce(sce, groupvar = "donor")
#' set <- append_groupvar_adj(set)
#' @export
append_groupvar_adj <- function(set, groupvar_pattern = "donor.*", 
                            assayname = "counts", type = "mgvdenom",
                            verbose = FALSE, ...){
  if(!(is(set, "SummarizedExperimentTypes"))){
    stop("set must be of class SummarizedExperimentTypes.")}
  if(type == "mgvdenom"){
    set <- groupadj_mgvdenom_fromrd(set = set, 
                                    groupvar_pattern = groupvar_pattern,
                                    rd.method.str = type, assayname = assayname,
                                    verbose = verbose, ...)
  } else{stop("Error: invalid adjustment type.")}
  return(set)
}

#' groupadj_mgvdenom_fromrd
#'
#' Get adjusted assays and rowadata, calculated as the counts over the mean 
#' group-wise variances ("mgvdenom"), or: 
#' 
#' $$counts/mgvdenom$$
#' 
#' @param set A SummarizedExperimentTypes
#' @param groupvar_pattern String pattern to identify group-level summaries from
#' rowData column labels with grepl.
#' @param assayname Name of original assays data to adjust, as identifiable from
#' names(assays(set)).
#' @param rd.method.str Character string to append for new assays and rowData 
#' matrices.
#' @param verbose Whether to return verbose status messages.
#' @returns set with updated rowData columns and new assays matrix.
#' @export
groupadj_mgvdenom_fromrd <- function(set, 
                                     groupvar_pattern = "donor.*",
                                     assayname = "counts",
                                     rd.method.str = "mgvdenom", 
                                     verbose = FALSE){
  if(verbose){message("Checking for groups in rowData...")}
  typev <- colnames(set); rd <- rowData(set); rd.cnv <- colnames(rd)
  str.patt <- paste0("^type[0-9];", groupvar_pattern, ";var$")
  rd.cnvf <- rd.cnv[grepl(str.patt, rd.cnv)]
  if(length(rd.cnvf)==0){stop("Error, no groups found in rowData.")}
  if(verbose){message("Group data found. Getting marker adjustments...")}
  rd.new <- do.call(cbind, lapply(typev, function(typei){
    type.str.patt <- paste0("^", typei, ";.*")
    which.rd.filt <- which(grepl(type.str.patt, rd.cnvf))
    rd.cnvf.filt <- rd.cnvf[which.rd.filt]
    rdf <- as.data.frame(rd[,rd.cnvf.filt])
    rowMeans(rdf, na.rm = T)
  }))
  colnames(rd.new) <- paste0(typev, ";", rd.method.str)
  if(verbose){message("Updating rowData...")}
  rowData(set) <- cbind(rd, rd.new) # update rowdata
  if(verbose){message("Getting new adjusted assay data...")}
  new.assay.name <- paste0(assayname,"_adj_",rd.method.str)
  assays(set)[[new.assay.name]] <- assays(set)[[assayname]]/rd.new
  return(set)
}

#-------------------
# 4. set simulations
#-------------------
