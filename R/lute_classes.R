#!/usr/bin/env R

#
# Classes supporting generic and framework functions.
#

#' cellProportionsPredictions-class
#'
#' Class for cell type predictions.
#' 
#' @details Main constructor for class \linkS4class{cellProportionsPredictions}.
#' @rdname cellProportionsPredictions-class
#' @param predictions.table Table containing cell type predictions.
#' @param cell.type.vector Character vector of cell type labels.
#' @param sample.id.vector Character vector of sample id labels.
#' @returns New cellProportionsPredictions object.
#' @examples 
#' new("cellProportionsPredictions")
#' ptable <- matrix(sample(100,50),nrow=10)
#' colnames(ptable) <- paste0("cell_type",seq(ncol(ptable)))
#' rownames(ptable) <- paste0("sample", seq(nrow(ptable)))
#' pred <- cellProportionsPredictions(ptable)
#' pred
setClass("cellProportionsPredictions", slots=c(predictions.table="data.frame",
                                                 cell.type.vector="character",
                                                 sample.id.vector="character"))

#' Make new cellProportionsPredictions object.
#' 
#' @param predictions.table Table of cell type predictions.
#' @param cell.type.vector Character vector of cell type labels.
#' @param sample.id.vector Character vector of sample id labels.
#' @returns cellProportionsPredictions object.
#' @importFrom methods new
#' @returns New cellProportionsPredictions object.
#' 
#' @examples
#' example.data <- get_decon_example_data()
#' 
#' @export
cellProportionsPredictions <- function(predictions.table, 
                                       cell.type.vector=NULL, 
                                       sample.id.vector=NULL) {
  if(is(cell.type.vector, "NULL")){
    cell.type.vector <- colnames(predictions.table) 
  }
  if(is(sample.id.vector, "NULL")){
    sample.id.vector <- rownames(predictions.table)
  }
  predictions.df <- as.data.frame(predictions.table)
  new("cellProportionsPredictions", 
      predictions.table=predictions.df,
      cell.type.vector=as.character(cell.type.vector),
      sample.id.vector=as.character(sample.id.vector))
}

#' Inspect cellProportionsPredictions object.
#' 
#' @param object cellProportionsPredictions object.
#' @importFrom methods show
#' @importFrom utils head
#' @details Method behavior for show.
#' @returns Shows object summaries.
#' 
#' @examples
#' example.data <- get_decon_example_data()
#' 
#' @export
setMethod("show", "cellProportionsPredictions", function(object) {
  ptable <- object@predictions.table
  # metadata summaries
  unique.cell.types.vector <- colnames(ptable)
  unique.sample.id.vector <- rownames(ptable)
  message("Number of bulk samples (J): ", 
          paste0(length(object@sample.id.vector), collapse="; "))
  message("Number of cell types (K): ", 
          paste0(length(object@cell.type.vector), collapse="; "))
  message("Cell type labels:\n", 
          paste0("\t", object@cell.type.vector, collapse="; "))
  # table summary
  print(
    message("predictions.table summary:\n"))
  print(head(ptable))
})
