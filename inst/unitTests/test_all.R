test_that("randomSingleCellExperiment produces expected SingleCellExperiment", {

	expect_equal(class(randomSingleCellExperiment())[1], "SingleCellExperiment")
	
	expect_equal(names(assays(randomSingleCellExperiment())), "counts")

})

test_that("nnlsParam produces expected results", {
  
  set.seed(0)
  
	exampleList <- getDeconvolutionExampleData()
	
	# evaluate param
	param <- nnlsParam(
	  cellScaleFactors=exampleList[["cellScaleFactors"]], 
	  bulkExpression=exampleList[["bulkExpression"]], 
	  referenceExpression=exampleList[["referenceExpression"]]
	)
	#
	# run tests
	expect_equal(class(param)[1], "nnlsParam")
	expect_equal(class(param@referenceExpression)[1], "matrix")
	expect_equal(class(param@bulkExpression)[1], "matrix")
	expect_equal(class(param@cellScaleFactors)[1], "numeric")
	expect_equal(names(assays(randomSingleCellExperiment())), "counts")
	expect_true(inherits(param, c("referencebasedParam", "deconvolutionParam")))
	
	# evalute deconvolution results
	deconvolutionResults <- deconvolution(param)
	predictionsTable <- deconvolutionResults@predictionsTable
	digitsRound <- 5
	# expect values
	expectType1Sample1 <- as.numeric(round(0.9985854, digitsRound))
	expectType2Sample1 <- round(0.001414605, digitsRound)
	expectType1Sample2 <- round(0.2231840, digitsRound)
	expectType2Sample2 <- round(0.776815983, digitsRound)
	# observed values
	observeType1Sample1 <- round(predictionsTable[1,1], digitsRound)
	observeType2Sample1 <- round(predictionsTable[1,2], digitsRound)
	observeType1Sample2 <- round(predictionsTable[2,1], digitsRound)
	observeType2Sample2 <- round(predictionsTable[2,2], digitsRound)
	#
	# run tests
	expect_equal(nrow(predictionsTable), 2)
	expect_equal(ncol(predictionsTable), 2)
	expect_equal(expectType1Sample1, observeType1Sample1)
	expect_equal(expectType2Sample1, observeType2Sample1)
	expect_equal(expectType1Sample2, observeType1Sample2)
	expect_equal(expectType2Sample2, observeType2Sample2)
	
})

test_that("bisqueParam produces expected results", {
  
  set.seed(0)
  
  exampleList <- getDeconvolutionExampleDataBisque()
  scData <- exampleList[["singleCellExpressionSet"]]
  scData$celltype <- scData$cellType
  scData$batch.id <- scData$SubjectName
  bulkExpressionSet <- exampleList$bulkExpressionSet
  phenoData(bulkExpressionSet)$batch.id <- colnames(bulkExpressionSet)
  
  # test param properties
  param <- bisqueParam(
    cellScaleFactors=c("type1" = 1, "type2" = 10), 
    bulkExpressionSet=bulkExpressionSet, 
    scData=scData
  )
  #
  # run tests
  expect_equal(class(param)[1], "bisqueParam")
  expect_equal(class(param@referenceExpression)[1], "matrix")
  expect_equal(class(param@bulkExpression)[1], "matrix")
  expect_equal(class(param@cellScaleFactors)[1], "numeric")
  expect_equal(names(assays(randomSingleCellExperiment())), "counts")
  expect_true(inherits(param, c("referencebasedParam", "deconvolutionParam")))
  
  # evalute deconvolution results
  deconvolutionResults <- deconvolution(param)
  predictionsTable <- deconvolutionResults@predictionsTable
  digitsRound <- 5
  # expect values
  expectType1Sample1 <- as.numeric(round(0.4417268, digitsRound))
  expectType2Sample1 <- round(0.5582732, digitsRound)
  expectType1Sample2 <- round(0.5099963, digitsRound)
  expectType2Sample2 <- round(0.4900037, digitsRound)
  # observed values
  observeType1Sample1 <- round(predictionsTable[1,1], digitsRound)
  observeType2Sample1 <- round(predictionsTable[1,2], digitsRound)
  observeType1Sample2 <- round(predictionsTable[2,1], digitsRound)
  observeType2Sample2 <- round(predictionsTable[2,2], digitsRound)
  #
  # run tests
  expect_equal(nrow(predictionsTable), 100)
  expect_equal(ncol(predictionsTable), 2)
  expect_equal(expectType1Sample1, observeType1Sample1)
  expect_equal(expectType2Sample1, observeType2Sample1)
  expect_equal(expectType1Sample2, observeType1Sample2)
  expect_equal(expectType2Sample2, observeType2Sample2)
  
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

test_that(
  paste0(
    "cellScaleFactors (k2, k3, k30) apply to the correct referenceExpression ",
    "columns"), 
  {
  
  
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
  
  #---------
  #
  # TEST K30
  #
  #---------
  K = 30
  # get transform
  exampleList <- getDeconvolutionExampleData(seq(K), numberTypes=K)
  transformResult <- lute:::.zstransform(
    exampleList[["referenceExpression"]], exampleList[["cellScaleFactors"]])
  #
  # get expected values
  expectedProductType1 <- as.numeric(
    exampleList$cellScaleFactors["type1"]*exampleList$referenceExpression[,"type1"])
  expectedProductType2 <- as.numeric(
    exampleList$cellScaleFactors["type2"]*exampleList$referenceExpression[,"type2"])
  expectedProductType30 <- as.numeric(
    exampleList$cellScaleFactors["type30"]*exampleList$referenceExpression[,"type30"])
  #
  # run tests
  expect_equal(transformResult[1,1], expectedProductType1[1])
  expect_equal(transformResult[2,1], expectedProductType1[2])
  expect_equal(transformResult[1,2], expectedProductType2[1])
  expect_equal(transformResult[2,2], expectedProductType2[2])
  expect_equal(transformResult[1,30], expectedProductType30[1])
  expect_equal(transformResult[2,30], expectedProductType30[2])
  
  #-------------------
  #
  # shuffle labels, K3
  #
  #-------------------
  K = 3
  # get transform
  exampleList <- getDeconvolutionExampleData(seq(K), numberTypes=K)
  s.vector.shuffle <- exampleList[["cellScaleFactors"]][rev(seq(3))]
  transformResult <- lute:::.zstransform(
    exampleList[["referenceExpression"]], s.vector.shuffle)
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

test_that(
  paste0("ypb_from_sce() performs properly."), {
  exampleData <- getDeconvolutionExampleData()
  exampleSingleCellExperiment <- randomSingleCellExperiment()
  examplePseudobulk <- ypb_from_sce(
    singleCellExperiment = exampleSingleCellExperiment)
  
  # run tests
  expect_equal(colnames(examplePseudobulk), "singleCellExperiment.pseudobulk")
  expect_equal(nrow(examplePseudobulk), nrow(exampleSingleCellExperiment))
  expect_equal(
    rownames(examplePseudobulk), rownames(exampleSingleCellExperiment))
  
})

