test_that("randomSingleCellExperiment produces expected SingleCellExperiment", {

	expect_true(is(randomSingleCellExperiment(), "SingleCellExperiment"))
	
	expect_equal(names(assays(randomSingleCellExperiment())), "counts")

})

test_that("nnlsParam produces expected results", {

	exampleList <- getDeconvolutionExampleData()
	param <- nnlsParam(
	  cellScaleFactors=exampleList[["cellScaleFactors"]], 
	  bulkExpression=exampleList[["bulkExpression"]], 
	  referenceExpression=exampleList[["referenceExpression"]]
	)

	expect_true(is(param, "nnlsParam"))

	expect_true(is(param@referenceExpression, "matrix"))

	expect_true(is(param@bulkExpression, "matrix"))

	expect_true(is(param@cellScaleFactors, "numeric"))

	expect_true(is(param, "referencebasedParam"))

	expect_true(is(param, "deconvolutionParam"))
	
	expect_equal(names(assays(randomSingleCellExperiment())), "counts")

})

test_that("cellScaleFactors apply to the correct referenceExpression columns", {
  
  exampleList <- getDeconvolutionExampleData()
  param <- nnlsParam(
    cellScaleFactors=exampleList[["cellScaleFactors"]], 
    bulkExpression=exampleList[["bulkExpression"]], 
    referenceExpression=exampleList[["referenceExpression"]]
  )
  deconvolutionResults <- deconvolution(param)
  
})