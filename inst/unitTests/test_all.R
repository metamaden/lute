test_that("randomSingleCellExperiment produces expected SingleCellExperiment", {

	expect_equal(class(randomSingleCellExperiment())[1], "SingleCellExperiment")
	
	expect_equal(names(assays(randomSingleCellExperiment())), "counts")

})

test_that("nnlsParam produces expected results", {

	exampleList <- getDeconvolutionExampleData()
	param <- nnlsParam(
	  cellScaleFactors=exampleList[["cellScaleFactors"]], 
	  bulkExpression=exampleList[["bulkExpression"]], 
	  referenceExpression=exampleList[["referenceExpression"]]
	)
	
	expect_equal(class(param)[1], "nnlsParam")

	expect_equal(class(param@referenceExpression)[1], "matrix")

	expect_equal(class(param@bulkExpression)[1], "matrix")

	expect_equal(class(param@cellScaleFactors)[1], "numeric")

	expect_equal(names(assays(randomSingleCellExperiment())), "counts")
	
	expect_true(inherits(param, c("referencebasedParam", "deconvolutionParam")))
	
})

test_that("deconvolution() results have expected structure.", {
  
  exampleList <- getDeconvolutionExampleData()
  param <- nnlsParam(
    cellScaleFactors=exampleList[["cellScaleFactors"]], 
    bulkExpression=exampleList[["bulkExpression"]], 
    referenceExpression=exampleList[["referenceExpression"]]
  )
  deconvolutionResults <- deconvolution(param)
  
  expect_equal(class(deconvolutionResults)[1], "cellProportionsPredictions")
  expect_equal(names(attributes(deconvolutionResults)[1]), "predictionsTable")
  expect_equal(names(attributes(deconvolutionResults)[2]), "cellTypeVector")
  expect_equal(names(attributes(deconvolutionResults)[3]), "sampleIdVector")
  expect_equal(names(attributes(deconvolutionResults)[4]), "class")
  
})

test_that("cellScaleFactors apply to the correct referenceExpression columns", {
  
  exampleList <- getDeconvolutionExampleData()
  param <- nnlsParam(
    cellScaleFactors=exampleList[["cellScaleFactors"]], 
    bulkExpression=exampleList[["bulkExpression"]], 
    referenceExpression=exampleList[["referenceExpression"]]
  )
  
  expect.vector.marker.1
  
  expect_true(identical())
  
  param@referenceExpression[,2]*param@cellScaleFactors[2]
  
})



