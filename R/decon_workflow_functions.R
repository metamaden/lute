#!/usr/bin/env R

# Author: Sean Maden
#
# Defines functions used by NextFlow processes for lute workflow.
#
#
#
#

#-----------------------
# sce summary statistics
#-----------------------

#' mexpr_na_freq
#'
#' Get NA frequency from a 2d matrix.
#' 
#' @param mexpr 2d matrix.
#' @param dim.index Integer for apply index to do NA counts. Either 
#' @returns list of vectors of NA statistics, including NA counts ("na.countv")
#' and NA frequencies ("na.freqv")
#' @export
mexpr_na_freq <- function(mexpr, dim.index = 2){
  na.countv <- apply(mexpr, 1, function(ci){length(which(is.na(ci)))})
  na.freqv <- na.countv/nrow(mexpr)
  return(list(na.countv = na.countv, na.freqv = na.freqv))
}

#-------------------------
# preprocessing -- filters
#-------------------------

#' filter_value_cells
#'
#' Filter cells on NA frequency.
#'
#' @param sce A SingleCellExperiment object.
#' @param remove.cells Whether to remove cells exceeding NA filter from sce 
#' prior to returning.
#' @param max.na.freq The upper threshold for tolerated NA frequency. Cells with 
#' NA frequency exceeding this are flagged for removal.
#' @param assayname Name of assay in assays(sce) to check for NAs.
#' @param append.metadata Whether to append metadata summary of preprocessing 
#' parameters and results.
#' @param verbose Whether to include verbose status messages.
#' @returns scef, filtered SingleCellExperiment with optional appended 
#' preprocessing metadata.
#' @examples 
#' sce <- random_sce(na.include = T, na.fract = 0.4, 
#' zero.include = T, zero.fract = 0.4)
#' scef <- filter_value_cells(sce, verbose = T)
#' @export
filter_value_cells <- function(sce, remove.cells = TRUE, max.na.freq = 0.25, 
                            assayname = "counts", append.metadata = TRUE, 
                            verbose = FALSE){
  if(!is(sce, "SingleCellExperiment")){
    stop("Error, sce must be a SingleCellExperiment.")}
  cd <- colData(sce); mexpr <- assays(sce)[[assayname]]
  # count nas
  na.countv <- apply(mexpr, 2, function(ci){length(which(is.na(ci)))})
  na.freqv <- na.countv/nrow(mexpr); filt.cellv <- na.freqv > max.na.freq
  scef <- sce[,!filt.cellv] # filter sce
  if(verbose){message("Filter on cell NA values removed ",
                      ncol(mexpr)-ncol(mexprf)," cells.")}
  if(append.metadata){
    if(verbose)(message("Appending new metadata."))
    dfmd <- data.frame(
      cell.uid = colnames(mexpr), na.count = na.countv,
      na.freq = na.freqv, above.max.na.freq = filt.cellv
    )
    lmd <- list(dfmd = dfmd, max.na.freq = max.na.freq, assayname = assayname)
    metadata(scef)$cell.na.filt <- lmd
  }
  return(scef)
}

#' filter_value_type
#'
#' Filter cell types (or other variable groupings) on frequency of genes or
#' features with cells exceeding some frequency threshold (see max.cells.na.freq).
#'
#' @param sce A SingleCellExperiment object.
#' @param filter.term Term for filter. Either "zerocount" or "NA".
#' @param type.variable Name of the type variable in sce colData.
#' @param remove.types Whether to remove cells exceeding NA filter from sce 
#' prior to returning.
#' @param max.gene.na.freq The upper threshold for frequency of NA cells 
#' (see max.na.freq) in a given type, where genes above this count towards the
#' type gene NA frequency.
#' @param max.cells.na.freq The upper threshold for tolerated NA frequency. Types with 
#' NA frequency exceeding this are flagged for removal.
#' @param assayname Name of assay in assays(sce) to check for NAs.
#' @param append.metadata Whether to append metadata summary of preprocessing 
#' parameters and results.
#' @param new.metadata.name Name of new metadata object to append (see 
#' append.metadata argument).
#' @param verbose Whether to include verbose status messages.
#' @returns scef, filtered SingleCellExperiment with optional appended 
#' preprocessing metadata.
#' @examples 
#' set.seed(0)
#' sce <- random_sce(zero.include = T, zero.fract = 0.3)
#' scef <- filter_value_type(sce, filter.term = "zerocount", verbose = T,
#' max.gene.value.freq = 0.15)
#' metadata(scef)$filter.zerocount.by.type$df.type # view new metadata
#' @export
filter_value_type <- function(sce, filter.term = c("zerocount", "NA"), 
                              type.variable = "celltype", remove.types = TRUE,
                              max.gene.value.freq = 0.25, max.cells.value.freq = 0.25,
                              assayname = "counts", append.metadata = TRUE,
                              new.metadata.name = NULL,
                              verbose = FALSE){
  if(!is(sce, "SingleCellExperiment")){
    stop("Error, sce must be a SingleCellExperiment.")}
  cd <- colData(sce)
  if(!type.variable %in% colnames(cd)){
    stop("Error, type.variable must be a column name in sce colData.")
  }
  mexpr <- assays(sce)[[assayname]]
  # count nas by type
  typev <- unique(cd[,type.variable])
  # get cell value stats
  df.cell.ct <- do.call(rbind, lapply(typev, function(ti){
    filt.type <- cd[,type.variable]==ti; mexprf <- mexpr[,filt.type]
    if(filter.term == "NA"){
      apply(mexprf, 1, function(ci){length(which(is.na(ci)))})
    } else{
      apply(mexprf, 1, function(ci){length(which(ci==0))})
    }
  }))
  df.cell.freq <- df.cell.ct/ncol(mexpr)
  df.cell.thresh <- t(apply(df.cell.freq, 1, function(ri){
    ri > max.cells.value.freq}))
  rownames(df.cell.ct) <- rownames(df.cell.freq) <- 
    rownames(df.cell.thresh) <- typev
  # get stats by type
  typev.gene.count <- unlist(apply(df.cell.thresh, 1, function(ri){
    length(which(ri))}))
  typev.gene.freq <- typev.gene.count/ncol(df.cell.freq)
  # get filters and apply them
  typev.gene.filt <- typev.gene.freq > max.gene.value.freq
  typev.remove <- typev[typev.gene.filt]
  coldata.filt <- cd[,type.variable] %in% typev.remove
  cell.uid.filt <- colnames(sce) %in% rownames(cd[coldata.filt,])
  scef <- sce[,!cell.uid.filt]
  if(verbose){message("After applying NA filter, removed ",length(typev.remove),
                      " types and ",ncol(sce)-ncol(scef)," cells...")}
  if(append.metadata){
    if(verbose){message("Appending new metadata...")}
    df.type <- data.frame(type = typev, gene.na.count = typev.gene.count,
                          gene.na.freq = typev.gene.freq,
                          eval.max.gene.value.freq = typev.gene.filt)
    lparam <- list(
      filter.term = filter.term,
      num.genes.sce.original = nrow(sce),
      num.cells.sce.original = ncol(sce),
      max.cells.value.freq = max.cells.value.freq,
      max.gene.value.freq = max.gene.value.freq
    )
    lmd <- list(parameters = lparam, df.type = df.type,
                df.cell.ct = df.cell.ct, df.cell.freq = df.cell.freq,
                df.cell.thresh = df.cell.thresh)
    if(is(new.metadata.name, "NULL")){
      new.metadata.name <- paste0("filter.", filter.term, ".by.type")
    }
    metadata(scef)[[new.metadata.name]] <- lmd
  }
  return(scef)
}

#--------------
# normalization
#--------------

#-----------------
# bias adjustments
#-----------------