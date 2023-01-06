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

# assign new group var "donor"
colData(sce)$donor <- rep(c("donor1", "donor2"), ncol(sce)/2)

# show donor by type
table(sce$celltype, sce$donor)
#       donor1 donor2
# type1      3      2
# type2      2      3

#------------------
# make set from sce
#------------------
# set <- set_from_sce(sce, typevar = "celltype", method = "mean")

sce_groupstat <- function(scef, groupvar, assayname = "counts", 
                          summarytype = "rowData",
                          groupstat = c("mean", "median", 
                                        "var", "sd", "numzero"),
                          verbose = FALSE){
  if(groupvar %in% colnames(colData(scef))){
    if(verbose){
      message("Appending group-level statistics for variable: ", groupvar)}
    num.groups <- length(unique(scef[[groupvar]]))
    if(num.groups > 1){
      ugroupv <- unique(scef[[groupvar]])
      dfg <- do.call(cbind, lapply(ugroupv, function(groupi){
        if(verbose){message("summarizing group: ", groupi)}
        sceff <- scef[,scef[[groupvar]]==groupi]
        exprff <- as.matrix(assays(sceff)[[assayname]])
        if(summarytype == "colData"){
          exprff <- matrix(apply(t(exprff),2,mean), nrow = 1)
        }
        dfgi <- data.frame(num.entries = rep(ncol(exprff), nrow(exprff)))
        if("mean" %in% groupstat){dfgi$mean <- rowMeans(exprff)}
        if("median" %in% groupstat){dfgi$median <- rowMedians(exprff)}
        if("var" %in% groupstat){dfgi$var <- rowVars(exprff)}
        if("sd" %in% groupstat){dfgi$sd <- rowSds(exprff)}
        if("numzero" %in% groupstat){
          dfgi$numzero <- unlist(
            apply(exprff, 1, function(ri){length(which(ri==0))}))
        }
        if(summarytype == "colData"){
          colnames(dfgi) <- paste0(groupi, ";colData_means;", colnames(dfgi))
        } else{
          colnames(dfgi) <- paste0(groupi, ";", colnames(dfgi))
        }
        dfgi
      }))
      return(dfg)
    } else{
      if(verbose){
        message("One group detected; skipping group statistics...")}
    }
  } else{
    message("Warning: variable ",groupvar," not found in sce coldata.")
  }
  return(NULL)
}

# script
groupvar <- 'donor' 
method <- "mean"
typevar <- "celltype"
assayname <- "counts"
verbose <- TRUE
typev <- unique(sce[[typevar]])
expr.sce <- assays(sce)$counts
expr.set <- do.call(cbind, lapply(typev, function(typei){
  if(verbose){message("Summarizing type: ", typei, "...")}
  scef <- sce[,sce[[typevar]]==typei]
  exprf <- assays(scef)[[assayname]]
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
  de <- data.frame(var = gene.varv, sdv = gene.sdv, max = gene.max,
                   min = gene.min)
  # parse group-level statistics
  if(!is(groupvar, "NULL")){
    dfg <- sce_groupstat(scef = scef, groupvar = groupvar, 
                         summarytype = "rowData",
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
colnames(mexpr) <- gsub(";.*", "", colnames(mexpr))
rd <- expr.set[,!which.mexpr] # rowdata

# get coldata
cd.sce <- colData(sce)
cd <- do.call(rbind, lapply(typev, function(typei){
  if(verbose){message("Working on type: ", typei, "...")}
  scef <- sce[,sce[[typevar]]==typei]
  num.cells <- ncol(scef)
  # parse 0-expr genes/cells
  count.zeroexpr <- apply(assays(scef)$counts, 1, function(ri){length(ri[ri==0])})
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
  # parse group-level statistics
  if(!is(groupvar, "NULL")){
    dfg <- sce_groupstat(scef = scef, groupvar = groupvar, assayname = assayname, 
                         summarytype = "colData", verbose = verbose)
    condv <- is(dfg, "data.frame") & nrow(dfg) == nrow(dfr)
    if(condv){
      if(verbose){message("Binding group-level data.")};dfr <- cbind(dfr, dfg)
    }
  }
  return(dfr)
}))


# get summary stats
do.call(rbind, lapply())

set.assays <- Assays(SimpleList(as.matrix(mexpr)))

set <- SummarizedExperimentTypes(assays = list(mexpr),
                                 rowData = rd, colData = cd)



