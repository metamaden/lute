test_that("random_sce works properly", {

	expect_true(is(random_sce(), "SingleCellExperiment"))
	
	expect_equal(names(assays(random_sce())), "counts")

})

test_that("nnlsParam works properly", {

	lexample <- get_decon_example_data()
	param <- nnlsParam(s = lexample[["s"]], y = lexample[["y"]], z = lexample[["z"]])

	expect_true(is(param, "nnlsParam"))

	expect_true(is(param@z, "matrix"))

	expect_true(is(param@y, "matrix"))

	expect_true(is(param@s, "numeric"))

	expect_true(is(param, "referencebasedParam"))

	expect_true(is(param, "deconvolutionParam"))
	
	expect_equal(names(assays(random_sce())), "counts")

})