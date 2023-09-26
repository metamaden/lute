test_that("random_sce works properly", {

	checkTrue(is(random_sce(), "SingleCellExperiment"))
	
	checkEquals(names(assays(random_sce())), "counts")

})

test_that("nnlsParam works properly", {

	lexample <- get_decon_example_data()
	param <- nnlsParam(s = lexample[["s"]], y = lexample[["y"]], z = lexample[["z"]])

	checkTrue(is(param, "nnlsParam"))

	checkTrue(is(param@z, "matrix"))

	checkTrue(is(param@y, "matrix"))

	checkTrue(is(param@s, "numeric"))

	checkTrue(is(param, "referencebasedParam"))

	checkTrue(is(param, "deconvolutionParam"))

	checkFalse(is(param, "independentbulkParam"))
	
	checkEquals(names(assays(random_sce())), "counts")

})