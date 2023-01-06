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
# set random seed
set.seed(0)

# example se object
# get params
num.genes <- 20
num.cells <- 10
num.types <- 2
expr.mean <- 10
# get expr
expr.ct <- matrix(rpois(num.cells*num.genes, lambda = expr.mean), 
                 ncol=num.cells, nrow=num.genes)
# get coldata
cellv <- paste0("cell.barcode.", seq(num.cells))
cpertype <- round(num.cells/num.types, 0)
typev <- c(rep("type1", 5), rep("type2", cpertype))
cd <- data.frame(cell.id = cellv,
                 celltype = typev)
colnames(expr.ct) <- cellv
# get rowdata
genev <- paste0("gene", seq(nrow(expr.ct)))
rd <- data.frame(gene.id = genev)
rownames(expr.ct) <- genev

# make sce object
sce <- SingleCellExperiment(assays = list(counts=expr.ct), 
                            colData = cd,
                            rowData = rd)

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



