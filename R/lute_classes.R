#!/usr/bin/env R

#
# Classes supporting generic and framework functions.
#

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
#' new("cellProportionsPredictions")
#' 
#' @aliases 
#' MeanratiosParam-class, MeanRatiosParam-class
#' 
setClass("cellProportionsPredictions", contains = "tbl", 
         slots = c(predictions.table = "tbl"))

cellProportionsPredictions <- function(predictions.table, return.info = FALSE) {
  predictions.table <- predictions.table %>% as_tibble()
  new("cellProportionsPredictions", predictions.table=predictions.table)
}

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

ptable <- matrix(sample(100,50),nrow=10)
colnames(ptable) <- paste0("cell_type",seq(ncol(ptable)))
rownames(ptable) <- paste0("sample", seq(nrow(ptable)))
ptable <- as_tibble(ptable)
pred <- cellProportionsPredictions(ptable)
pred
