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

#' filter_na_cells
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
#' sce <- random_sce(na.include = T, na.fract = 0.4)
#' scef <- filter_na_cells(sce, verbose = T)
#' @export
filter_na_cells <- function(sce, remove.cells = TRUE, max.na.freq = 0.25, 
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

#' filter_na_type
#'
#' Filter cell types (or other variable groupings) on frequency of genes or
#' features with cells exceeding some frequency threshold (see max.cells.na.freq).
#'
#' @param sce A SingleCellExperiment object.
#' @param filterterm Term for filter. Either "zerocount" or "NA".
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
#' # use to filter cell types
#' sce <- random_sce(na.include = T, na.fract = 0.4)
#' scef <- filter_na_type(sce, verbose = T)
#' metadata(scef)$filter.na.type$df.type # view summary by type
#' # can also be used to filter on grouping, such as by donor
#' colData(sce)$donor <- c(
#'    rep("donor1", 3), 
#'    rep("donor2", 5), 
#'    rep("donor3", 2)
#' )
#' scef <- filter_na_cells(sce, type.variable = "donor", 
#'    new.metadata.name = "filter.na.donor", verbose = T)
#' @export
filter_value_type <- function(sce, type.variable = "celltype", remove.types = TRUE, 
                           max.gene.na.freq = 0.25, max.cells.na.freq = 0.25, 
                           assayname = "counts", append.metadata = TRUE, 
                           new.metadata.name = "filter.na.type",
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
  # get cell na stats
  dfna.ct <- do.call(rbind, lapply(typev, function(ti){
    filt.type <- cd[,type.variable]==ti; mexprf <- mexpr[,filt.type]
    apply(mexprf, 1, function(ci){length(which(is.na(ci)))})
  }))
  dfna.freq <- dfna.ct/ncol(mexpr)
  dfna.thresh <- t(apply(dfna.freq, 1, function(ri){ri > max.cells.na.freq}))
  rownames(dfna.ct) <- rownames(dfna.freq) <- rownames(dfna.thresh) <- typev
  # get stats by type
  typev.gene.count <- unlist(apply(dfna.thresh, 1, function(ri){
    length(which(ri))}))
  typev.gene.freq <- typev.gene.count/ncol(dfna.freq)
  # get filters and apply them
  typev.gene.filt <- typev.gene.freq > max.gene.na.freq
  typev.remove <- typev[typev.gene.filt]
  coldata.filt <- cd[,type.variable] %in% typev.remove
  cell.uid.filt <- colnames(sce) %in% rownames(cd[coldata.filt,])
  scef <- sce[,!cell.uid.filt]
  if(verbose){message("After applying NA filter, removed ",length(typev.remove),
                      " types and ",ncol(sce)-ncol(scef)," cells...")}
  if(append.metadata){
    if(verbose){message("Appending new metadata...")}
    dftype <- data.frame(type = typev,
                         gene.na.count = typev.gene.count,
                         gene.na.freq = typev.gene.freq,
                         eval.max.gene.na.freq = typev.gene.filt)
    lparam <- list(
      num.genes.sce.original = nrow(sce),
      num.cells.sce.original = ncol(sce),
      max.cells.na.freq = max.cells.na.freq,
      max.gene.na.freq = max.gene.na.freq
    )
    lmd <- list(parameters = lparam, df.type = dftype,
                dfna.cell.ct = dfna.ct, dfna.cell.freq = dfna.freq,
                dfna.cell.thresh = dfna.thresh)
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