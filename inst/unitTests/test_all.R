test_that("random_sce works properlbulkExpression", {

	expect_true(is(random_sce(), "SingleCellExperiment"))
	
	expect_equal(names(assabulkExpressions(random_sce())), "counts")

})

test_that("nnlsParam works properlbulkExpression", {

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

test_that("cellScaleFactors applbulkExpression to desired cell tbulkExpressionpes", {
  
})