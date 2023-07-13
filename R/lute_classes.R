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
#' @examples 
#' new("cellProportionsPredictions")
#' # 
#' ptable <- matrix(sample(100,50),nrow=10)
#' colnames(ptable) <- paste0("cell_type",seq(ncol(ptable)))
#' rownames(ptable) <- paste0("sample", seq(nrow(ptable)))
#' ptable <- as_tibble(ptable)
#' pred <- cellProportionsPredictions(ptable)
#' pred
#' @aliases 
#' MeanratiosParam-class, MeanRatiosParam-class
#' 
setClass("cellProportionsPredictions", slots = c(predictions.table = "tbl"))

#' Make new cellProportionsPredictions object.
#' 
#' @param predictions.table Table of cell type predictions.
#' @returns cellProportionsPredictions object.
#' @export
cellProportionsPredictions <- function(predictions.table) {
  predictions.table <- predictions.table %>% as_tibble()
  new("cellProportionsPredictions", predictions.table=predictions.table)
}

#' Inspect cellProportionsPredictions object.
#' @param object cellProportionsPredictions object.
#' @details Method behavior for show.
#' @export
setMethod("show", "cellProportionsPredictions", function(object) {
  ptable <- object@predictions.table
  # metadata summaries
  unique.cell.types.vector <- colnames(ptable)
  unique.sample.id.vector <- rownames(ptable)
  message("Number of bulk samples (J): ", 
          paste0(length(unique.sample.id.vector), collapse = "; "))
  message("Number of cell types (K): ", 
          paste0(length(unique.cell.types.vector), collapse = "; "))
  message("Cell type labels:\n", 
          paste0("\t", unique.cell.types.vector, collapse = "; "))
  # table summary
  print(
    message("predictions.table summary:\n"))
  print(head(ptable))
})
