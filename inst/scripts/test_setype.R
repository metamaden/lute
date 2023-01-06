#!/usr/bin/env R

# Author: Sean Maden
#
# Test methods for SummarizedExperimentTypes class
#
#

libv <- c("SingleCellExperiment", "lute")
sapply(libv, library, character.only = TRUE)

#----------------
# make sce object
#----------------
sce <- random_sce()

#------------------
# make set from sce
#------------------
set <- set_from_sce(sce, typevar = "celltype", method = "mean")

# script
method = "mean"
typevar <- "celltype"
typev <- unique(sce[[typevar]])
expr.sce <- assays(sce)$counts
expr.set <- do.call(cbind, lapply(typev, function(typei){
  exprf <- expr.sce[,sce[[typevar]]==typei]
  gene.varv <- apply(exprf,1,var)
  gene.sdv <- apply(exprf,1,sd)
  gene.max <- apply(exprf,1,max)
  gene.min <- apply(exprf,1,min)
  if(method == "mean"){
    exprnew <- matrix(rowMeans(exprf), ncol = 1)
  } else if(method == "median"){
    exprnew <- matrix(rowMedians(exprf), ncol = 1)
  } else{
    exprnew <- exprf[,1,drop=F]
  }
  colnames(exprnew) <- paste0(typei, ";expr")
  de <- data.frame(var = gene.varv,
                   sdv = gene.sdv,
                   max = gene.max,
                   min = gene.min)
  colnames(de) <- paste0(typei, ";", colnames(de))
  return(cbind(exprnew, de))
}))
which.mexpr <- grepl(".*;expr$", colnames(expr.set))
mexpr <- expr.set[,which.mexpr] # expr data
colnames(mexpr) <- gsub(";.*", "", colnames(mexpr))
rd <- expr.set[,!which.mexpr] # rowdata

# get coldata
cd.sce <- colData(sce)
cd <- do.call(rbind, lapply(typev, function(typei){
  scef <- sce[,sce[[typevar]]==typei]
  num.cells <- ncol(scef)
  # parse 0-expr genes/cells
  count.zeroexpr <- apply(assays(scef)$counts, 1, function(ri){length(ri[ri==0])})
  num.allzeroexpr <- length(which(count.zeroexpr==num.cells))
  mean.zerocount <- mean(count.zeroexpr)
  median.zerocount <- median(count.zeroexpr)
  var.zerocount <- var(count.zeroexpr)
  sd.zerocount <- sd(count.zeroexpr)
  # insert group-wise summary code here
  group.wise.mean <- group.wise.var <- num.groups <- NA
  data.frame(type = typei,
             num.cells = num.cells,
             num.allzeroexpr = num.allzeroexpr,
             mean.zerocount = mean.zerocount,
             median.zerocount = median.zerocount,
             var.zerocount = var.zerocount,
             sd.zerocount = sd.zerocount,
             num.groups = num.groups,
             group.wise.mean = group.wise.mean,
             group.wise.var = group.wise.var)
}))

# get summary stats
do.call(rbind, lapply())

set.assays <- Assays(SimpleList(as.matrix(mexpr)))

set <- SummarizedExperimentTypes(assays = list(mexpr),
                                 rowData = rd, colData = cd)



