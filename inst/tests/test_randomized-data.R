#!/usr/bin/env R

# Author: Sean Maden
#
#

require(lute)

#------------------
# test random_sce()
#------------------
num.genes = 20
num.cells = 100
num.types = 2
fract.types = NULL
expr.mean = 10 
na.include = FALSE
na.fract = 0.2
zero.include = FALSE
zero.fract = 0.2
verbose = FALSE
seed.num = 0



if(verbose){message("Getting random expression data...")}
mdat <- rnbinom(n = num.cells*num.genes, size = expr.mean, mu = expr.mean)
if(na.include){ # manually add NAs
  if(verbose){message("Including NA values...")}
  num.na <- round(length(mdat)*na.fract, digits = 0)
  na.index <- sample(seq(length(mdat)), num.na)
  mdat[na.index] <- NA
}
if(zero.include){ # manually add zero counts
  if(verbose){message("Including NA values...")}
  num.zero <- round(length(mdat)*zero.fract, digits = 0)
  zero.index <- sample(seq(length(mdat)), num.zero)
  mdat[zero.index] <- 0
}
expr.ct <- matrix(mdat, ncol=num.cells, nrow=num.genes)
if(verbose){message("Getting new colData...")}
cellv <- paste0("cell.barcode.", seq(num.cells))
cpertype <- round(num.cells/num.types, 0)

if(is(fract.types, "NULL")){
  fract.types = rep((1/num.types), num.types)}
typev <- paste0("type", seq(num.types))
typev <- unlist(lapply(seq(length(typev)), function(ti){
  num <- fract.types[ti]*num.cells; message(num)
  rep(typev[ti], num)
}))

cd <- data.frame(cell.id = cellv, celltype = typev)
colnames(expr.ct) <- cellv
if(verbose){message("Getting new rowData...")}
genev <- paste0("gene", seq(nrow(expr.ct)))
rd <- data.frame(gene.id = genev)
rownames(expr.ct) <- genev
if(verbose){message("Making new sce object...")}
sce <- SingleCellExperiment(assays = list(counts=expr.ct), colData = cd, rowData = rd)


