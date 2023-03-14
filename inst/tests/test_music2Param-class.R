library(lute)
# source("./R/lute_utilities.R")

#--------------------
# final class example
#--------------------
# example params
batch.variable <- "SubjectName"
celltype.variable <- "cellType"

# get data
lexample <- lute:::.get_decon_example_data_music2()
sc.eset <- lexample[["sc.eset"]]
sc.sce <- lute:::.get_sce_from_eset(sc.eset)
y.eset <- lexample[["y.eset"]]
y <- exprs(y.eset)
unique.types <- unique(sc.sce[[celltype.variable]])
unique.types <- unique.types[order(unique.types)]
s <- rep(1, length(unique.types))
names(s) <- unique.types
cell_size <- data.frame(cell_type = unique.types, cell_size = s)

# get yi
id.sc <- unique(sc.eset[[batch.variable]])
id.bulk <- colnames(y.eset)
filter <- id.bulk[!id.bulk %in% id.sc]
yi <- exprs(y.eset)[,filter]

result <- MuSiC::music2_prop(bulk.control.mtx = y, bulk.case.mtx = yi,
                             sc.sce = sc.sce, clusters = celltype.variable,
                             samples = batch.variable, cell_size = cell_size,
                             select.ct = NULL)


# get param object
param <- music2Param(bulk.control.mtx = y, bulk.case.mtx = yi,
                     sc.sce = sc.sce, clusters = celltype.variable,
                     samples = batch.variable, cell_size = s)

# get just predictions
res <- deconvolution(param)

# get full results
param@return.info <- TRUE
res <- deconvolution(param)

#---------------------
# test music2 function
#---------------------

result <- MuSiC::music2_prop()
