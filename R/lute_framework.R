#!/usr/bin/env r

# Author: Sean Maden
#
# Framework functions and utilities for the lute package.

#' lute framework
#'
#' Perform a deconvolution experiment.
#' 
#' @param sce Object of type SingleCellExperiment. Optional (see z).
#' @param z Signature matrix of cell type-specific signals. Optional (see sce).
#' @param y Bulk mixed signals matrix of samples, which can be matched to 
#' single-cell samples. Optional (see y.se).
#' @param y.se SummarizedExperiment or similar data type containing the bulk signals
#' matrix in its assays (e.g. accessible with assays(y.se)[[assay.name]] using
#' the provided assay.name argument). Optional (see y).
#' @param s Cell size factor transformations of length equal to the K cell types 
#' to deconvolve. Optional, if not provided, uses equal weights for types.
#' @param assay.name Name of expression matrix in sce, and optionally y.se, 
#' assays. Optional (see sce, y.se).
#' @param celltype.variable Name of cell type variable in sce coldata.
#' @param markers.per.type Number of top markers to get per cell type.
#' @param return.info Whether to return metadata and original method outputs 
#' with predicted proportions.
#' @param typemarker.algorithm Which type-specific marker selection algorthm to 
#' use. If NULL, skips type marker analyses.
#' @param deconvolution.algorithm Where deconvolution algorithm to use. If NULL, 
#' skips deconvolution.
#' @param verbose Whether to show verbose status messages.
#' 
#' @returns A list containing results returned from type marker selection and
#' deconvolution runs, with additional information returned if
#'  \code{return.info == TRUE}.
#' 
#' @details Main function to use the lute deconvolution framework. Manages data
#' conversions and mappings to deconvolution experiment steps, including setup,
#' gene marker identification, and main deconvolution runs.
#' 
#' Support is provided for \linkS4class{SummarizedExperiment}-type or 
#' matrix-type inputs for the Z signature matrix and Y bulk signals matrix. Note,
#' both Z and Y need to be provided or derivable in order to run deconvolution.
#' 
#' @examples 
#' # get example bulk data
#' y <- lute:::.get_decon_example_data()$y
#' 
#' # get example sce
#' sce <- random_sce()[seq(10),]
#' 
#' # get framework results
#' experiment.results <- lute(sce = sce, y = y)
#' 
#' @export
lute <- function(sce = NULL, z = NULL, y = NULL, y.se = NULL, s = NULL, 
                 return.info = FALSE, markers.per.type = 20,
                 assay.name = "counts", celltype.variable = "celltype",
                 typemarker.algorithm = "meanratio", 
                 deconvolution.algorithm = "nnls",
                 verbose = TRUE){
  results.list <- list()
  if(!is(sce, "NULL")){marker.argument.valid <- TRUE}
  if(marker.argument.valid & !is(typemarker.algorithm, "NULL")){
    if(is(sce, "NULL")){
      stop("Error, provide sce to perform typemarkers analyses.")}
    if(verbose){message("Parsing marker gene arguments...")}
    typemarker.results <- marker.vector <- map_typemarker_algorithm(
      algorithm = typemarker.algorithm,
      sce = sce,
      assay.name = assay.name,
      celltype.variable = celltype.variable,
      markers.per.type = markers.per.type,
      return.info = return.info)
    if(is(typemarker.results, "list")){
      marker.vector <- typemarker.results[["markers"]]}
    if(verbose){message("Filtering sce...")}
    sce.filter <- rownames(sce) %in% marker.vector
    sce <- sce[sce.filter,]
    z <- .get_z_from_sce(sce, assay.name, celltype.variable)
    results.list[["typemarker.results"]] <- typemarker.results
  }
  y.cond <- is(y, "NULL") & !is(y.se, "NULL")
  if(y.cond){y <- assays(y.se)[[assay.name]]}
  if(!is(deconvolution.algorithm, "NULL")){
    if(is(z, "NULL")){
      if(is(sce, "NULL")){
        stop("Error, provide either sce or z to perform deconvolution.")
      } else{
        z <- .get_z_from_sce(sce, assay.name, celltype.variable)
      }
    }
    if(is(y, "NULL")){
      stop("Error, provide y to perform deconvolution.")
    }
    if(is(s, "NULL")){s <- rep(1, ncol(z))}
    if(verbose){message("Parsing deconvolution arguments...")}
    deconvolution.results <- map_deconvolution_algorithm(
      algorithm = deconvolution.algorithm,
      z = z, y = y, s = s, return.info = return.info)
    results.list[["deconvolution.results"]] <- deconvolution.results
  }
  return(results.list)
}

#'
#'
map_typemarker_algorithm <- function(algorithm, sce, assay.name, 
                                     celltype.variable, markers.per.type, 
                                     return.info){
  if(algorithm %in% c("mr", "MR", "Mr", "deconvobuddies", "DeconvoBuddies", 
                         "meanratio", "meanRatio", "Meanratio", "MeanRatio",
                         "meanratios", "meanRatios", "Meanratio", "MeanRatio")){
    message("Using meanratiosParam...")
    typemarker.string <- "meanratiosParam"
  } else{
    message("Warning, unidentified marker selection algorithm provided. ",
            "Skipping marker selection")
    return.string <- "FALSE"
  }
  typemarker.string <- paste0(typemarker.string, 
                              "(sce = sce, assay.name = '",assay.name,"', ",
                              "markers.per.type = ", markers.per.type, ", ",
                              "celltype.variable = '",celltype.variable,"',",
                              "return.info = ", return.info, ")")
  new.param <- eval(parse(text = typemarker.string))
  return(typemarkers(new.param))
}

#'
#'
map_deconvolution_algorithm <- function(algorithm, y = y, z = z, s = s, 
                                        return.info = return.info){
  if(algorithm %in% c("nnls", "Nnls", "NNLS", "nnlsParam", "NNLSParam", "NnlsParam")){
    message("Using NNLS...")
    deconvolution.string <- "nnlsParam"
  } else if(algorithm %in% c("music", "Music", "MuSiC", "musicParam", "MusicParam", "MuSiCParam")){
    message("Using MuSiC...")
    deconvolution.string <- "musicParam"
  } else if(algorithm %in% c("music2", "Music2", "MuSiC2", "musicParam2", "MusicParam2", "MuSiCParam2")){
    message("Using MuSiC2...")
    deconvolution.string <- "music2Param"
  } else if(algorithm %in% c("epic", "Epic", "EPIC", "epicParam", "EpicParam", "EPICParam")){
    message("Using EPIC...")
    deconvolution.string <- "epicParam"
  } else if(algorithm %in% c("bisque", "Bisque", "BISQUE", "bisqueParam", "BisqueParam", "BISQUEParam")){
    message("Using Bisque...")
    deconvolution.string <- "bisqueParam"
  } else if(algorithm %in% c("deconrnaseq", "Deconrnaseq", "DECONRNASEQ", 
                                "deconrnaseqParam", "DeconrnaseqParam", "DECONRNASEQParam")){
    message("Using DeconRNASeq...")
    deconvolution.string <- "deconrnaseqParam"
  } else if(algorithm %in% c("scdc", "Scdc", "SCDC", 
                                "scdcParam", "ScdcParam", "SCDCParam")){
    message("Using SCDC...")
    deconvolution.string <- "scdcParam"
  } else{
    stop("Error, unidentified deconvolution algorithm provided. ")
  }
  deconvolution.string <- paste0(deconvolution.string,
                                 "(y = y, z = z, s = s, ",
                                 "return.info = ", return.info, ")")
  new.param <- eval(parse(text = deconvolution.string))
  return(deconvolution(new.param))
}
