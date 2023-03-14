library(lute)
# source("./R/lute_utilities.R")

#------------------
# load example data
#------------------
# example params
batch.variable <- "SubjectName"
celltype.variable <- "cellType"

# get data
lexample <- lute:::.get_decon_example_data_bisque()
sc.eset <- lexample[["sc.eset"]]
y.eset <- lexample[["y.eset"]]

#-----------------
# test example data
#-----------------



# get z from sc.eset

sce <- SingleCellExperiment(assays = list(counts = exprs(sc.eset)))
colData(sce) <- DataFrame(pData(sc.eset))
z <- lute:::.get_z_from_sce(sce = sce, assay.name = "counts", 
                            celltype.variable = celltype.variable)

#-----------------------------------
# test bisqueParam with example data
#-----------------------------------
param <- bisqueParam(y.eset = y.eset, sc.eset = sc.eset,
                     batch.variable = "SubjectName",
                     celltype.variable = "cellType")
res <- deconvolution(param)

#-------------------------------------
# test example data on original method
#-------------------------------------
results <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset = lexample[["y.eset"]], 
                                                  sc.eset = lexample[["sc.eset"]])
# proportions <- results$bulk.props
# proportions.sample1 <- proportions$bulk.props[,1]
