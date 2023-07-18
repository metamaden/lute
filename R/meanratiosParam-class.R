#!/usr/bin/env R

# Author: Sean Maden

#' nnlsParam-class
#'
#' Uses DeconvoBuddies::get_mean_ratio2()
#' 
#' @include lute_generics.R
#' @include typemarkersParam-class.R
#' 
#' @details Main constructor for class \linkS4class{meanratiosParam}.
#' @rdname meanratiosParam-class
#' @seealso \linkS4class{typemarkersParam}
#' 
#' @param assay.name Name of expression matrix in sce assays.
#' @param sce Object of type SingleCellExperiment.
#' @param celltype.variable Name of cell type variable in sce coldata.
#' 
#' @examples 
#' lexample <- lute:::.get_decon_example_data()
#' sce.example <- random_sce()
#' #new.param <- meanratiosParam(sce = sce.example, celltype.variable = "celltype", markers.per.type = 20)
#' #markers <- typemarkers(new.param)
#' 
#' @aliases 
#' MeanratiosParam-class, MeanRatiosParam-class
#' 
setClass("meanratiosParam", contains="typemarkersParam", 
         slots = c(assay.name = "character", 
                   sce = "SingleCellExperiment",
                   celltype.variable = "character"))

#' Make new object of class meanratiosParam
#'
#' Main constructor for class \linkS4class{meanratiosParam}.
#'
#' @param assay.name Name of expression matrix in sce assays.
#' @param sce Object of type SingleCellExperiment.
#' @param celltype.variable Name of cell type variable in sce coldata.
#' @param markers.per.type Number of top markers to get per cell type.
#' @param return.info Whether to return metadata and original method outputs 
#' with predicted proportions.
#' 
#' @returns Object of class \linkS4class{meanratiosParam}
#' 
#' @seealso \linkS4class{typemarkersParam}
#'
#' @details Main class for mapping arguments to the mean ratios method 
#' implemented as \code{DeconvoBuddies::get_mean_ratio2()}.
#' 
#' @export
meanratiosParam <- function(sce, assay.name = "counts", 
                            celltype.variable = "cellType", 
                            markers.per.type = 20, 
                            return.info = FALSE) {
  new("meanratiosParam", 
      sce = sce, assay.name = assay.name, celltype.variable = celltype.variable,
      markers.per.type = markers.per.type, return.info = return.info)
}

#' Cell type markers method for meanratiosParam
#'
#' Defines the typemarkers method for \linkS4class{meanratiosParam}.
#'
#' @param object An object of class \linkS4class{meanratiosParam}.
#'
#' @details Takes an object of class \linkS4class{meanratiosParam} as input, 
#' returning either a vector of cell type gene markers, or (if 
#' \code{return.info == TRUE}) a list containing such a vector along with 
#' original function outputs.
#'
#' @returns Either a vector of gene markers, or a list containing such a vector 
#' with the original method outputs.
#'
#' @export
setMethod("typemarkers", signature(object = "meanratiosParam"), function(object){
  require(DeconvoBuddies)
  sce <- object[["sce"]]
  celltype.variable <- object[["celltype.variable"]]
  assay.name = object[["assay.name"]]
  markers.per.type <- object[["markers.per.type"]]
  # get marker results
  marker.table <- DeconvoBuddies::get_mean_ratio2(sce = sce,
                                                  cellType_col = celltype.variable,
                                                  assay_name = assay.name,
                                                  add_symbol = FALSE) %>%
    as.data.frame()
  # filter top markers
  unique.cell.types <- sce[[celltype.variable]] %>% as.character() %>% unique()
  
  top.markers.list <- lapply(unique.cell.types, function(unique.type.id){
    marker.table %>% 
      dplyr::filter(cellType.target == unique.type.id) %>% 
      dplyr::arrange(rank_ratio) %>% 
      dplyr::top_n(n = markers.per.type)
  })
  top.markers.table <- do.call(rbind, top.markers.list)
  top.markers.vector <- top.markers.table$gene
  # parse return.info
  return.list <- top.markers.vector
  if(object[["return.info"]]){
    return.list <- list(markers = top.markers.vector,
                        result.info = marker.table,
                        metadata = object)}
  return(return.list)
})
