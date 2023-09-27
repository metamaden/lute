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

# set condition.variable
condition.variable <- "condition"
control.label <- "control"
case.label <- "case"
# assign controls
pData(y.eset)[,condition.variable] <- control.label
colData(sc.sce)[,condition.variable] <- control.label
# assign cases
pData(y.eset)[filter, condition.variable] <- case.label

#---------------------------
# test MuSiC2 implementation
#---------------------------
result <- MuSiC2::music2_prop(bulk.eset = y.eset, 
                              sc.eset = sc.sce, 
                              condition = condition.variable,
                              control = control.label,
                              case = case.label,
                              clusters = celltype.variable,
                              samples = batch.variable, 
                              cell_size = cell_size,
                              select.ct = unique.types)

#--------------------------
# test MuSiC implementation
#--------------------------
result <- MuSiC::music2_prop(bulk.control.mtx = y, bulk.case.mtx = yi,
                             sc.sce = sc.sce, clusters = celltype.variable,
                             samples = batch.variable, cell_size = cell_size,
                             select.ct = NULL)

#----------------------
# test vignette example
#----------------------
# source: https://jiaxin-fan.github.io/MuSiC2/articles/introduction.html

benchmark.eset = readRDS("./ignore/bulk-eset.rds")
seger.eset = readRDS("./ignore/single-eset.rds")

est.prop = music2_prop(bulk.eset = benchmark.eset, 
                       sc.eset = seger.eset, condition='group', 
                       control='healthy',case='t2d', 
                       clusters = 'cellType', 
                       samples = 'sampleID', 
                       select.ct = c('acinar','alpha','beta','delta','ductal','gamma'), 
                       n_resample=20, 
                       sample_prop=0.5,cutoff_c=0.05,cutoff_r=0.01)$Est.prop

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
