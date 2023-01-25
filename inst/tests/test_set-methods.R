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


typev <- unique(sce[[type.variable]])

# get new assays matrix
ma <- do.call(cbind, lapply(typev, function(typei){
  if(verbose){message("Summarizing type: ", typei, "...")}
  type.filt <- sce[[type.variable]]==typei
  scef <- sce[,type.filt]
  exprf <- assays(scef)[[assayname]]
  mai <- matrix(rowMeans(exprf), ncol = 1)
  colnames(mai) <- typei
  return(mai)
}))
rownames(ma) <- rownames(sce)
lma <- list(ma); names(lma) <- paste0("summarized_", assayname)

# parse metadata
# get rowdata
rd <- sce_groupstat(sce, group.variable = type.variable, 
                    summarytype = "rowData", return.tall = FALSE)
marker.filt <- grepl(".*;marker$", colnames(rd))
rownames(rd) <- rd[,which(marker.filt)[1]]
rd$marker <- rd[,which(marker.filt)[1]]
rd <- rd[,!marker.filt]
rd <- rd[order(match(rd$marker, rownames(ma))),]
cond.rd <- identical(rd$type, rownames(ma))
if(!cond.rd){
  message("Warning, couldn't match coldata types to ma types.")
  rd <- matrix(nrow = nrow(ma), ncol = 0)
}
# get coldata
cd <- sce_groupstat(sce, group.variable = type.variable, 
                    summarytype = "colData", return.tall = TRUE)
rownames(cd) <- cd$group; cd$type <- cd$group
cd <- cd[,!colnames(cd)=="group"]
cd <- cd[order(match(cd$type, colnames(ma))),]
cond.cd <- identical(cd$type, colnames(ma))
if(!cond.cd){
  message("Warning, couldn't match coldata types to ma types.")
  cd <- matrix(nrow = ncol(ma), ncol = 0)
}

# metadata
lmd <- list(assay.info = list(
  stat.method = method, sce.assayname = assayname, 
  type.variable = type.variable, group.variable = group.variable
))

# make new set object
set <- SummarizedExperimentTypes(assays = lma, rowData = rd, 
                                 colData = cd, metadata = lmd)

return(set)



