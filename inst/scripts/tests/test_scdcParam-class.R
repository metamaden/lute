require(lute)
# get data
lexample <- lute:::.get_decon_example_data_scdc()
sc.eset <- lexample[["sc.eset"]]
y.eset <- lexample[["y.eset"]]
# get yi
id.sc <- unique(sc.eset[["SubjectName"]])
id.bulk <- colnames(y.eset)
filter <- id.bulk[!id.bulk %in% id.sc]
yi <- exprs(y.eset)[,filter]
# get param object
param <- scdcParam(yi = yi, y.eset = y.eset, sc.eset = sc.eset,
                   batch.variable = "SubjectName", 
                   celltype.variable = "cellType",
                   celltype.subset = NULL)
# get just predictions
result <- deconvolution(param)
