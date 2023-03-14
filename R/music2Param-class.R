#' music2Param-class
#' 
#' Main constructor for class to manage mappings to the deconvolution method function \code{MuSiC::music2_prop()}.
#' 
#' @include lute_generics.R
#' @include deconParam-class.R
#' @include referencebasedParam-class.R
#' @include independentbulkParam-class.R
#' 
#' @details Main constructor for class \linkS4class{music2Param}.
#' @rdname music2Param-class
#' @seealso \linkS4class{deconParam}, \linkS4class{referencebasedParam}, \linkS4class{independentbulkParam}
#'
#' @examples 
#' lexample <- .get_decon_example_data()
#' 
#' @aliases 
#' music2Param-class, MuSiC2Param-class, Music2Param-class
#'
setClass("music2Param", contains="independentbulkParam", slots=c(
	sc.eset = "ExpressionSet", assay.name = "character", batch.variable = "character", 
	celltype.variable = "character", iter.max = "numeric", nu = "numeric", epsilon = "numeric", 
	truep = "numeric", ct.cell.size = "numeric", celltype.subset = "character"))