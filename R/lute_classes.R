#!/usr/bin/env R

### Author: Sean Maden
###
### Classes supporting generic and framework functions.
###

#' cellProportionsPredictions-class
#'
#' Class for cell type predictions.
#' 
#' @details Main constructor for class \linkS4class{cellProportionsPredictions}.
#' @rdname cellProportionsPredictions-class
#' @param predictionsTable Table containing cell type predictions.
#' @param cellTypeVector Character vector of cell type labels.
#' @param sampleIdVector Character vector of sample id labels.
#' @returns New cellProportionsPredictions object.
#' @examples 
#' new("cellProportionsPredictions")
#' predictionsTable <- matrix(sample(100,50),nrow=10)
#' colnames(predictionsTable) <- paste0("cell_type",seq(ncol(predictionsTable)))
#' rownames(predictionsTable) <- paste0("sample", seq(nrow(predictionsTable)))
#' cellProportionsPredictions(predictionsTable)
setClass("cellProportionsPredictions", slots=c(predictionsTable="data.frame",
                                                 cellTypeVector="character",
                                                 sampleIdVector="character"))

#' Make new cellProportionsPredictions object.
#' 
#' @param predictionsTable Table of cell type predictions.
#' @param cellTypeVector Character vector of cell type labels.
#' @param sampleIdVector Character vector of sample id labels.
#' @returns New cellProportionsPredictions object.
#' @importFrom methods new
#' @returns New cellProportionsPredictions object.
#' 
#' @examples
#' exampleData <- get_decon_example_data()
#' 
#' @export
cellProportionsPredictions <- function(predictionsTable, 
                                       cellTypeVector=NULL, 
                                       sampleIdVector=NULL) {
  if(is(cellTypeVector, "NULL")){
    cellTypeVector <- colnames(predictionsTable) 
  }
  if(is(sampleIdVector, "NULL")){
    sampleIdVector <- rownames(predictionsTable)
  }
  dfPredictions <- as.data.frame(predictionsTable)
  new("cellProportionsPredictions", 
      predictionsTable=dfPredictions,
      cellTypeVector=as.character(cellTypeVector),
      sampleIdVector=as.character(sampleIdVector))
}

#' Inspect cellProportionsPredictions object.
#' 
#' @param object Object of type cellProportionsPredictions (see 
#' \code{?cellProportionsPredictions}).
#' @importFrom methods show
#' @importFrom utils head
#' @details Method behavior for show.
#' @returns Shows object summaries.
#' 
#' @examples
#' exampleData <- get_decon_example_data()
#' 
#' @export
setMethod("show", "cellProportionsPredictions", function(object) {
  predictionsTable <- object@predictionsTable
  ## metadata summaries
  message("Number of bulk samples (J): ", 
          paste0(length(object@sampleIdVector), collapse="; "))
  message("Number of cell types (K): ", 
          paste0(length(object@cellTypeVector), collapse="; "))
  message("Cell type labels:\n", 
          paste0("\t", object@cellTypeVector, collapse="; "))
  ## table summary
  print(
    message("predictionsTable summary:\n"))
  print(head(predictionsTable))
})
