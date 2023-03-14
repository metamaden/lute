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

#' Deconvolution method for scdcParam
#'
#' Main method to access the SCDC deconvolution method from the main lute deconvolution genetic.
#'
#' @details Takes an object of class scdcParam as input, returning a list or vector of predicted 
#' cell type proportions.
#'
#' @returns Either a vector of predicted proportions, or a list containing predictions, metadata, 
#' and original outputs.
#'
#' @export
setMethod("deconvolution", signature(object = "scdcParam"), function(object){
  require(SCDC); require(Biobase)
  # load data
  lparam <- callNextMethod()
  object <- lparam[["object"]]

  # instantiate function objects
  y.eset <- lparam[["y.eset"]]
  sc.eset <- object@sc.eset
  celltype.subset <- object@celltype.subset
  batch.variable <- object@batch.variable
  iter.max <- object@iter.max
  nu <- object@nu
  epsilon <- object@epsilon
  truep <- object@truep
  s <- object@s
  weight.basis <- object@weight.basis
  transform.bisque <- object@transform.bisque

  # get result according to method type
  source.libarary <- object@method.type
  if(source.library == "MuSiC"){
  	result <- MuSiC::music2_prop(bulk.control.mtx = y, bulk.case.mtx = yi, sc.sce = sc.sce, 
  		clusters = celltype.variable, samples = batch.variable, cell_size = cell_size, select.ct = NULL)
  } else{
	result <- MuSiC2::music2_prop(bulk.eset = y.eset, sc.eset = sc.sce, condition = condition.variable,
	                      control = control.label, case = case.label, clusters = celltype.variable,
	                      samples = batch.variable, cell_size = cell_size, select.ct = unique.types)
  }

  # get predictions
  result <- SCDC::SCDC_prop(bulk.eset = y.eset, sc.eset = sc.eset,
  	ct.varname = celltype.variable, sample = batch.variable,
  	iter.max = iter.max, nu = nu, epsilon = epsilon, truep = truep,
  	ct.cell.size = s, ct.sub = celltype.subset, weight.basis = weight.basis,
  	Transform_bisque = transform_bisque)

  # return results
  lr <- predictions <- result$bulk.props
  if(object[["return.info"]]){
    lr <- list(predictions = predictions, result.info = result, 
               metadata = list(lmd = lparam[["metadata"]], 
                y.eset = y.eset, sc.eset = sc.eset))}
  return(lr)
})