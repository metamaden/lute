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

test_that("bisqueParam produces expected results", {
  
  exampleList <- getDeconvolutionExampleDataBisque()
  scData <- exampleList[["singleCellExpressionSet"]]
  scData$celltype <- scData$cellType
  scData$batch.id <- scData$SubjectName
  bulkExpressionSet <- exampleList$bulkExpressionSet
  phenoData(bulkExpressionSet)$batch.id <- colnames(bulkExpressionSet)
  
  param <- bisqueParam(
    cellScaleFactors=c("type1" = 1, "type2" = 10), 
    bulkExpressionSet=bulkExpressionSet, 
    scData=scData
  )
  
  expect_equal(class(param)[1], "bisqueParam")
  
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

test_that("Example function returns expected results for K=2, 3, and 30", {
  
  # test k2
  numTypes <- 2
  exampleK2 <- getDeconvolutionExampleData(seq(numTypes), numberTypes=2)
  expect_equal(ncol(exampleK2$referenceExpression), numTypes)
  
  # test k3
  numTypes <- 3
  exampleK3 <- getDeconvolutionExampleData(seq(numTypes), numberTypes=3)
  expect_equal(ncol(exampleK3$referenceExpression), numTypes)
  
  # test k30
  numTypes <- 30
  exampleK30 <- getDeconvolutionExampleData(seq(numTypes), numberTypes=30)
  expect_equal(ncol(exampleK30$referenceExpression), numTypes)
  
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
  exampleList <- getDeconvolutionExampleData(seq(3), numberTypes=3)
  transformResult <- lute:::.zstransform(
    exampleList[["referenceExpression"]], exampleList[["cellScaleFactors"]])
  #
  # get expected values
  expectedProductType1 <- as.numeric(
    exampleList$cellScaleFactors["type1"]*exampleList$referenceExpression[,"type1"])
  expectedProductType2 <- as.numeric(
    exampleList$cellScaleFactors["type2"]*exampleList$referenceExpression[,"type2"])
  #
  # run tests
  expect_equal(transformResult[1,1], expectedProductType1[1])
  expect_equal(transformResult[2,1], expectedProductType1[2])
  expect_equal(transformResult[1,2], expectedProductType2[1])
  expect_equal(transformResult[2,2], expectedProductType2[2])
  
})
