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

test_that("Example function returns expected results for K2,K3,K30", {
  
})

test_that("cellScaleFactors apply to the correct referenceExpression columns", {
  
  
  #---------
  #
  # TEST K2
  #
  #---------
  #
  # get transform
  exampleList <- getDeconvolutionExampleData()
  referenceExpression <- exampleList[["referenceExpression"]]
  cellScaleFactors <- exampleList[["cellScaleFactors"]]
  transformResult <- lute:::.zstransform(referenceExpression, cellScaleFactors)
  #
  # get expected values
  factorType1 <- exampleList$cellScaleFactors["type1"]
  referenceExpressionType1 <- exampleList$referenceExpression[,"type1"]
  expectedProductType1 <- as.numeric(factorType1*referenceExpressionType1)
  factorType2 <- exampleList$cellScaleFactors["type2"]
  referenceExpressionType2 <- exampleList$referenceExpression[,"type2"]
  expectedProductType2 <- as.numeric(factorType2*referenceExpressionType2)
  #
  # run tests
  expect_equal(transformResult[1,1], expectedProductType1[1])
  expect_equal(transformResult[2,1], expectedProductType1[2])
  expect_equal(transformResult[1,2], expectedProductType2[1])
  expect_equal(transformResult[2,2], expectedProductType2[2])
  
  #---------
  #
  # TEST K3
  #
  #---------
  #
  # get transform
  exampleList <- getDeconvolutionExampleData(numberTypes=3)
  referenceExpression <- exampleList[["referenceExpression"]]
  cellScaleFactors <- exampleList[["cellScaleFactors"]]
  transformResult <- lute:::.zstransform(referenceExpression, cellScaleFactors)
  #
  # get expected values
  factorType1 <- exampleList$cellScaleFactors["type1"]
  referenceExpressionType1 <- exampleList$referenceExpression[,"type1"]
  expectedProductType1 <- as.numeric(factorType1*referenceExpressionType1)
  factorType2 <- exampleList$cellScaleFactors["type2"]
  referenceExpressionType2 <- exampleList$referenceExpression[,"type2"]
  expectedProductType2 <- as.numeric(factorType2*referenceExpressionType2)
  #
  # run tests
  expect_equal(transformResult[1,1], expectedProductType1[1])
  expect_equal(transformResult[2,1], expectedProductType1[2])
  expect_equal(transformResult[1,2], expectedProductType2[1])
  expect_equal(transformResult[2,2], expectedProductType2[2])
  
})



