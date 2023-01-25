#!/usr/bin/env R

# Author: Sean Maden
#
# Testing methods for SummarizedExperimentTypes objects
#


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



if(!(is(sce, "SingleCellExperiment")|is(sce, "SummarizedExperiment"))){
  stop("sce must be of class SingleCellExperiment or SummarizedExperiment.")}
typev <- unique(sce[[type.variable]])
expr.sce <- assays(sce)$counts
if(!is(group.variable, "NULL")){
  ugroupv <- unique(sce[[group.variable]]) # get all possible group lvls
}


expr.set <- do.call(cbind, lapply(typev, function(typei){
  if(verbose){message("Summarizing type: ", typei, "...")}
  scef <- sce[,sce[[type.variable]]==typei]
  exprf <- assays(scef)[[assayname]]
  gene.varv <- apply(exprf,1,var)
  gene.sdv <- apply(exprf,1,sd)
  gene.max <- apply(exprf,1,max)
  gene.min <- apply(exprf,1,min)
  # get new assay data
  exprnew <- make_new_assaydata(exprf, method = method, 
                                na.rm = T, verbose = verbose)
  exprnew <- matrix(exprnew, ncol = 1)
  colnames(exprnew) <- paste0(typei, ";expr")
  de <- data.frame(var = gene.varv, sdv = gene.sdv, max = gene.max,
                   min = gene.min)
  # parse group-level statistics
  if(!is(group.variable, "NULL")){
    dfg <- sce_groupstat(sce = sce, group.variable = group.variable,
                         type.variable = type.variable, 
                         summarytype = "rowData", return.tall = FALSE, 
                         assayname = assayname, verbose = verbose)
    condv <- is(dfg, "data.frame") & nrow(dfg) == nrow(de)
    if(condv){if(verbose){message("Binding group-level data.")}
      de <- cbind(de, dfg)}
  }
  colnames(de) <- paste0(typei, ";", colnames(de))
  return(cbind(exprnew, de))
}))

which.mexpr <- grepl(".*;expr$", colnames(expr.set))
mexpr <- expr.set[,which.mexpr] # expr data
colnames(mexpr) <- gsub(";expr.*", "", colnames(mexpr))
rd <- expr.set[,!which.mexpr] # rowdata

# get coldata
cd.sce <- colData(sce)
cd <- do.call(rbind, lapply(typev, function(typei){
  if(verbose){message("Working on type: ", typei, "...")}
  scef <- sce[,sce[[type.variable]]==typei]
  num.cells <- ncol(scef)
  # parse 0-expr genes/cells
  count.zeroexpr <- apply(assays(scef)$counts, 1, 
                          function(ri){length(ri[ri==0])})
  num.allzeroexpr <- length(which(count.zeroexpr==num.cells))
  mean.zerocount <- mean(count.zeroexpr)
  median.zerocount <- median(count.zeroexpr)
  var.zerocount <- var(count.zeroexpr)
  sd.zerocount <- sd(count.zeroexpr)
  # get return df
  dfr <- data.frame(type = typei,
                    num.cells = num.cells,
                    num.allzeroexpr = num.allzeroexpr,
                    mean.zerocount = mean.zerocount,
                    median.zerocount = median.zerocount,
                    var.zerocount = var.zerocount,
                    sd.zerocount = sd.zerocount)
  return(dfr)
}))
# parse group-level statistics
if(!is(group.variable, "NULL")){
  dfg <- sce_groupstat(sce = sce, group.variable = group.variable,
                       type.variable = type.variable, return.tall = FALSE,
                       assayname = assayname, summarytype = "colData", 
                       verbose = verbose)
  cd <- cbind(dfr, dfg)
}

# metadata
lmd <- list(assay.info = list(
  stat.method = method, sce.assayname = assayname, 
  type.variable = type.variable, group.variable = group.variable
))
# make new set object
lassays <- list(mexpr); names(lassays) <- paste0("summarized_", assayname)
new.set <- SummarizedExperimentTypes(assays = lassays, rowData = rd, 
                                     colData = cd, metadata = lmd)
# parse standard plot options
if(make.set.plots){
  lp <- get_set_plots(set = new.set, group.variable = group.variable,
                      type.variable = type.variable, verbose = verbose)
  metadata(new.set)[["set_plots"]] <- lp
}