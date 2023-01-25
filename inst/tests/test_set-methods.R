#!/usr/bin/env R

# Author: Sean Maden
#
# Testing methods for SummarizedExperimentTypes objects
#

require(lute)

#-------------
# set_from_sce
#-------------

sce = random_sce()
colData(sce)$donor <- c(rep("donor1", 7), rep("donor2", 3))

group.variable = "donor"
method = "mean"
type.variable = "celltype"
assayname = "counts"
make.set.plots = TRUE
verbose = FALSE


# run checks
if(!(is(sce, "SingleCellExperiment")|is(sce, "SummarizedExperiment"))){
  stop("sce must be of class SingleCellExperiment or SummarizedExperiment.")}

# get new assays matrix
typev <- unique(sce[[type.variable]])
ma <- do.call(cbind, lapply(typev, function(typei){
  if(verbose){message("Summarizing type: ", typei, "...")}
  type.filt <- sce[[type.variable]]==typei; scef <- sce[,type.filt]
  exprf <- assays(scef)[[assayname]]
  if(ncol(exprf) > 0){
    mai <- make_new_assaydata(exprf, method = "mean", na.rm = TRUE, 
                              verbose = verbose)
    mai <- matrix(mai, ncol = 1); colnames(mai) <- typei
    return(mai)
  }
}))
rownames(ma) <- rownames(sce)
lma <- list(ma); names(lma) <- paste0("summarized_", assayname)

# parse metadata
# get rowdata
rd <- sce_groupstat(sce, group.variable = type.variable, assayname = assayname,
                    summarytype = "rowData", return.tall = FALSE)
rownames(rd) <- rownames(ma)

# get coldata
cd <- sce_groupstat(sce, group.variable = type.variable, assayname = assayname,
                    summarytype = "colData", return.tall = TRUE)
rownames(cd) <- colnames(ma)

# metadata
lmd <- list(assay.info = list(
  stat.method = method, 
  sce.assayname = assayname, 
  type.variable = type.variable
))

# make new set object
set <- SummarizedExperimentTypes(assays = lma, rowData = rd, 
                                 colData = cd, metadata = lmd)



