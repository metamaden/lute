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
                   celltype.subset = c("type1", "type2"))
# get just predictions
result <- deconvolution(param)
# get full results
param@return.info <- TRUE
result <- deconvolution(param)

#-------------------
# test scdc function
#-------------------
batch.variable <- "SubjectName"
celltype.variable <- "cellType"

lexample <- lute:::.get_decon_example_data_scdc()
sc.eset <- lexample[["sc.eset"]]
y.eset <- lexample[["y.eset"]]
nu <- 1e-4
epsilon <- 0.01
iter.max <- 1000
truep <- NULL
weight.basis <- T
unique.types <- unique(sc.eset[[celltype.variable]])
celltype.subset <- unique.types[1:2]
s <- rep(1, length(unique.types))
names(s) <- unique.types

result <- SCDC::SCDC_prop(bulk.eset = y.eset, 
                          sc.eset = sc.eset,
                          ct.varname = celltype.variable, 
                          sample = batch.variable,
                          iter.max = iter.max, 
                          nu = nu, epsilon = epsilon, 
                          truep = truep,
                          ct.cell.size = s, 
                          ct.sub = celltype.subset,
                          Transform_bisque = F)


