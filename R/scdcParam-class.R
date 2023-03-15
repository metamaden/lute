#' scdcParam-class
#' 
#' Main constructor for class to manage mappings to the deconvolution method function \code{SCDC::SCDC_prop()}.
#' 
#' @include lute_generics.R
#' @include independentbulkParam-class.R
#' 
#' @details Main constructor for class \linkS4class{scdcParam}.
#' @rdname scdcParam-class
#' @seealso \linkS4class{deconParam}, \linkS4class{referencebasedParam}, \linkS4class{independentbulkParam}
#'
#' @examples 
#' lexample <- .get_decon_example_data()
#' 
#' @references 
#' 
#' Meichen Dong. SCDC: Bulk RNA-seq deconvolution by scRNA-seq with multi-reference 
#' datasets. (2023) GitHub, R package version 0.0.0.9000. URL: https://github.com/meichendong/SCDC
#' 
#' Meichen Dong, Aatish Thennavan, Eugene Urrutia, Yun Li, Charles M Perou, Fei 
#' Zou, Yuchao Jiang, SCDC: bulk gene expression deconvolution by multiple 
#' single-cell RNA sequencing references, Briefings in Bioinformatics, Volume 
#' 22, Issue 1, January 2021, Pages 416–427, 
#' https://doi-org.proxy1.library.jhu.edu/10.1093/bib/bbz166
#' 
#' @aliases 
#' SCDCParam-class, ScdcParam-class, ScdcParam-class
#'
setClass("scdcParam", contains="independentbulkParam", slots=c(y.eset = "ExpressionSet", 
	sc.eset = "ExpressionSet", assay.name = "character", batch.variable = "character", 
	celltype.variable = "character", iter.max = "numeric", nu = "numeric", epsilon = "numeric", 
	truep = "numeric", ct.cell.size = "numeric", celltype.subset = "character"))

#' Make new object of class scdcParam
#'
#' Main constructor for class \linkS4class{scdcParam}.
#'
#' @param y Bulk mixed signals matrix of samples, which can be matched to single-cell samples.
#' @param yi Bulk mixed signals matrix of independent samples, which should not overlap samples in y.
#' @param z Signature matrix of cell type-specific signals. If not provided, can be computed from a
#' provided ExpressionSet containing single-cell data.
#' @param s Cell size factor transformations of length equal to the K cell types to deconvolve.
#' @param y.eset ExpressionSet of bulk mixed expression signals.
#' @param sc.eset ExpressionSet of single-cell transcriptomics data.
#' @param celltype.subset Vector of cell types to use for basis matrix.
#' @param assay.name Expression data type (e.g. counts, logcounts, tpm, etc.).
#' @param batch.variable Name of variable identifying the batches in sc.eset pData/coldata.
#' @param celltype.variable Name of cell type labels variable in sc.eset pData/coldata.
#' @param weight.basis Argument weight.basis for SCDC_prop().
#' @param iter.max Argument iter.max for SCDC_prop().
#' @param epsilon Argument epsilon for SCDC_prop().
#' @param nu Argument nu for SCDC_prop().
#' @param truep Argument truep for SCDC_prop().
#' @param return.info Whether to return metadata and original method outputs with predicted proportions.
#'
#' @details Takes standard inputs for the Bisque method. If user provides matrices, will convert these
#' into ExpressionSet objects compatible with the main bisque method.
#' 
#' @returns Object of class \linkS4class{scdcParam}.
#' 
#' @export
scdcParam <- function(y = NULL, yi = NULL, z = NULL, s = NULL, y.eset = NULL, sc.eset = NULL,
					  celltype.subset = NULL, assay.name = "counts", batch.variable = "batch.id", 
					  celltype.variable = "celltype", iter.max = 1000, nu = 1e-4, epsilon = 0.01,
					  weight.basis = TRUE, transform.bisque = FALSE, truep = NULL, return.info = FALSE) {
  require(Biobase)
  # check y.eset/y
  if(is(y, "NULL")){
    if(is(y.eset, "NULL")){
      stop("Error, need to provide either y or bulk.eset.")
    } else{
      message("Getting y from provided bulk.eset...")
      y <- as.matrix(exprs(y.eset))
    }
  } else{
      if(is(y.eset, "NULL")){
      message("Making ExpressionSet from provided y...")
      y.eset <- .make_eset_from_matrix(mat = y, batch.id = "SubjectName")
      # need at least 2 columns/samples to pass to bisque
      if(ncol(y.eset) == 1){
        sample.name <- colnames(y.eset)
        y.eset <- cbind(y.eset, y.eset)
        colnames(y.eset) <- c(sample.name, paste0(sample.name, "_rep1"))
      }
    }
  }

  # check sc.eset
  if(is(sc.eset, "NULL")){
    stop("Error, no single-cell ExpressionSet provided.")  
    # add condition to call splatter simulations by default?
  } else{
    if(!batch.variable %in% colnames(pData(sc.eset))){
    stop("Error, didn't find batch id variable ",batch.variable,
         " in sc.eset pData/coldata.")
    } else{
      id.sc <- unique(sc.eset[[batch.variable]])
    }
    if(!celltype.variable %in% colnames(pData(sc.eset))){
      stop("Error, didn't find celltype id variable ", celltype.variable, 
           " in sc.eset pData/coldata.")
    }
    if(is(z, "NULL")){
      message("Getting z from sc.eset...")
      sce <- .get_sce_from_eset(sc.eset)
      z <- .get_z_from_sce(sce = sce, celltype.variable = celltype.variable)
    }
  }

  # parse s
  unique.types <- colnames(z)
  unique.types <- unique.types[order(unique.types)]
  if(is(s, "NULL")){
  	message("Setting equal cell size factors...")
  	s <- rep(1, ncol(z))
  	names(s) <- unique.types
 }
  # parse ct.sub
  if(is(celltype.subset, "NULL")){
  	message("Using cell type labels for basis matrix.")
  	celltype.subset <- unique.types
  }
  
  # parse batch ids in bulk and sc
  message("Checking batch ids in bulk and sc eset...")
  if(cond <- !batch.variable %in% colnames(pData(y.eset))){
    stop("Error, didn't find batch variable in y.eset pData: ", batch.variable)
  } else{
    id.bulk <- unique(y.eset[[batch.variable]])
  }
  id.overlap <- intersect(id.sc, id.bulk)
  id.unique <- unique(c(id.sc, id.bulk))
  id.onlybulk <- id.bulk[!id.bulk %in% id.overlap]
  id.onlysc <- id.sc[!id.sc %in% id.overlap]
  message("Found ", length(id.unique), " unique batch ids...")
  message("Found ", length(id.overlap), " overlapping batch ids...")
  message("Found ", length(id.onlybulk), " bulk-only batch ids...")
  message("Found ", length(id.onlysc), " sc-only batch ids...")
  if(length(id.overlap) == 0){stop("Error, no overlapping markers in y.eset and sc.eset.")}
  
  # parse independent bulk samples
  if(length(id.onlybulk)==0){
    if(is(yi, "NULL")){
      stop("Error, no independent bulk samples found. ",
        "Provide either yi, or additional y samples.")
    } else{
      message("Using provided yi for independent bulk samples...")
    }
  } else{
    if(is(yi, "NULL")){
      message("Making yi from provided y bulk...")
      yi <- exprs(y.eset)[,colnames(y.eset) %in% id.onlybulk]
    } else{
      message("Using provided yi for independent bulk samples...")
    }
  }

  # check rowSums for basis matrix
  eset.basis <- scdc_basis_eset(sc.eset = sc.eset, ct.sub = ct.sub, 
  	ct.varname = ct.varname, min.sum = 0)
  if(nrow(eset.basis)==0){
  	stop("Error, no genes pass a minimum sum expression of 0 for the basis cell types")}

  # check number of cell types
  if(length(unique.types) < 4){
  	stop("Error, need at least 2 cell types for predictions")}
  if(length(intersect(unique.types, celltype.subset)) < 4){
  	stop("Error, need at least 2 cell types in celltype.subset, for basis matrix")}
  	
  new("scdcParam", y = y, yi = yi, z = z, s = s, y.eset = y.eset, sc.eset = sc.eset, 
  	  celltype.subset = celltype.subset, assay.name = assay.name, batch.variable = batch.variable, 
      celltype.variable = celltype.variable, return.info = return.info,
      iter.max = iter.max, nu = nu, epsilon = epsilon, truep = truep,
      weight.basis = weight.basis, transform.bisque = transform.bisque)
}

#' @export
scdc_basis_eset <- function(sc.eset, ct.sub, ct.varname, min.sum = 0){
	ct.sub <- ct.sub[!is.na(ct.sub)] # filter missing/NA types
	filter <- pData(sc.eset)[,ct.varname] %in% ct.sub
	eset.sub <- sc.eset[,filter]
	sums.vector <- rowSums(exprs(eset.sub))
	sums.filter <- sums.vector > min.sum
	return(eset.sub[sums.filter,])
}

#' Deconvolution method for scdcParam
#'
#' Main method to access the SCDC deconvolution method from the main lute deconvolution genetic.
#'
#' @details Takes an object of class scdcParam as input, returning a list or vector of predicted 
#' cell type proportions.
#'
#' @returns Either a vector of predicted proportions, or a list containing 
#' the predictions, metadata, and original outputs.
#' 
#' @references 
#' 
#' Meichen Dong. SCDC: Bulk RNA-seq deconvolution by scRNA-seq with multi-reference 
#' datasets. (2023) GitHub, R package version 0.0.0.9000. URL: https://github.com/meichendong/SCDC
#' 
#' Meichen Dong, Aatish Thennavan, Eugene Urrutia, Yun Li, Charles M Perou, Fei 
#' Zou, Yuchao Jiang, SCDC: bulk gene expression deconvolution by multiple 
#' single-cell RNA sequencing references, Briefings in Bioinformatics, Volume 
#' 22, Issue 1, January 2021, Pages 416–427, 
#' https://doi-org.proxy1.library.jhu.edu/10.1093/bib/bbz166
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
