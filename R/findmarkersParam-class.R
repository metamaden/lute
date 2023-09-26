#!/usr/bin/env R

# Author: Sean Maden

#' findmarkersParam-class
#'
#' class definition for findmarkersParam, which uses scran::findMarkers()
#' 
#' @include lute_generics.R
#' @include typemarkersParam-class.R
#' 
#' @details Main constructor for class \linkS4class{findmarkersParam}.
#' @rdname findmarkersParam-class
#' @seealso \linkS4class{typemarkersParam}
#' 
#' @param assay.name Name of expression matrix in sce assays.
#' @param sce Object of type SingleCellExperiment.
#' @param celltype.variable Name of cell type variable in sce coldata.
#' @param test.type Test type (see ?findMarkers for options).
#' 
#' @examples 
#' lexample <- get_decon_example_data()
#' sce.example <- random_sce()
#' new.param <- findmarkersParam(sce = sce.example, 
#' celltype.variable = "celltype", markers.per.type = 5)
#' markers <- typemarkers(new.param)
#' 
#' @aliases 
#' FindmarkersParam-class, findMarkersParam-class
#' 
#' @returns New object.
#'
setClass("findmarkersParam", contains="typemarkersParam", 
         slots = c(assay.name = "character", 
                   sce = "SingleCellExperiment",
                   celltype.variable = "character",
                   test.type = "character"))

#' Make new object of class findmarkersParam
#'
#' Main constructor for class \linkS4class{findmarkersParam}.
#'
#' @param assay.name Name of expression matrix in sce assays.
#' @param sce Object of type SingleCellExperiment.
#' @param celltype.variable Name of cell type variable in sce coldata.
#' @param markers.per.type Number of top markers to get per cell type.
#' @param test.type Test type (see ?findMarkers for options).
#' @param return.info Whether to return metadata and original method outputs 
#' with predicted proportions.
#' 
#' @returns Object of class \linkS4class{findmarkersParam}
#' 
#' @seealso \linkS4class{typemarkersParam}
#'
#' @details Main class for mapping arguments to the findMarkers method 
#' implemented as \code{scran::findMarkers()}.
#' 
#' @export
findmarkersParam <- function(sce, assay.name = "counts", 
                             celltype.variable = "cellType",
                             test.type = "wilcox",
                             markers.per.type = 20, 
                             return.info = FALSE) {
  new("findmarkersParam", 
      sce = sce, assay.name = assay.name, celltype.variable = celltype.variable,
      test.type = test.type, markers.per.type = markers.per.type, return.info = return.info)
}

#' Cell type markers method for findmarkersParam
#'
#' Defines the typemarkers method for \linkS4class{findmarkersParam}.
#'
#' @param object An object of class \linkS4class{findmarkersParam}.
#'
#' @details Takes an object of class \linkS4class{findmarkersParam} as input, 
#' returning either a vector of cell type gene markers, or (if 
#' \code{return.info == TRUE}) a list containing such a vector along with 
#' original function outputs.
#' 
#' @importFrom scran findMarkers
#' @importFrom dplyr %>%
#'
#' @returns Returns the top available markers, with type-specific marker filters,
#' as either a vector of marker IDs or a results list.
#'
#' @export
setMethod("typemarkers", signature(object = "findmarkersParam"), function(object){
  sce <- object[["sce"]]
  celltype.variable <- object[["celltype.variable"]]
  assay.name <- object[["assay.name"]]
  markers.per.type <- object[["markers.per.type"]]
  test.type <- object[["test.type"]]
  # get marker results
  list.result <- list.marker.table <- list(); markers.filter <- c()
  unique.cell.types <- sce[[celltype.variable]] %>% as.character() %>% unique()
  for(type in unique.cell.types){
    scef <- sce[!rownames(sce) %in% markers.filter,]
    message("selecting among ",nrow(scef)," genes for markers of type: ", type, "...")
    list.result[[type]] <- findMarkers(x = scef, 
                                        group = sce[[celltype.variable]],
                                        assay.type = assay.name,
                                        test.type = test.type)[[type]]
    dfi <- list.result[[type]][,seq(4)]
    summary.colname <- colnames(dfi)[grepl("summary\\..*", colnames(dfi))]
    new.colname <- paste0("abs.", summary.colname)
    dfi$abs.summary <- dfi[,summary.colname] %>% as.numeric() %>% abs()
    dfi$cellType.target <- type
    dfi$gene <- rownames(dfi)
    dfi <- as.data.frame(dfi)
    # get filtered markers
    list.marker.table[[type]] <- dfi %>% dplyr::arrange(abs.summary) %>% 
      dplyr::top_n(n = markers.per.type) %>% as.data.frame()
    markers.filter <- list.marker.table[[type]]$gene
  }
  top.marker.table <- do.call(rbind, list.marker.table) %>% as.data.frame()
  top.markers.vector <- top.marker.table$gene
  # parse return.info
  return.list <- top.markers.vector %>% unique()
  if(object[["return.info"]]){
    return.list <- list(markers = top.markers.vector,
                        result.info = marker.table,
                        metadata = object)}
  return(return.list)
})
