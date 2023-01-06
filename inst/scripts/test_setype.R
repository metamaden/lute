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
# two donors type1, one donor type2:
colData(sce)$donor <- c(rep("donor1", 7), rep("donor2", 3))
# make set from sce
set <- set_from_sce(sce, typevar = "celltype", method = "mean")

sce_groupstat <- function(scef, groupvar, ugroupv, assayname = "counts", 
                          summarytype = "rowData",
                          groupstat = c("mean", "median", 
                                        "var", "sd", "numzero"),
                          verbose = FALSE){
  groupstat.filt <- groupstat %in% c("mean", "median", "var", "sd", "numzero")
  groupstat <- groupstat[groupstat.filt]
  if(groupvar %in% colnames(colData(scef))){
    if(verbose){
      message("Appending group-level statistics for variable: ", 
              groupvar, "...")}
    dfg <- do.call(cbind, lapply(ugroupv, function(groupi){
      # make NA matrix
      dfgi <- as.data.frame(matrix(NA, 
                                   ncol = length(groupstat)+1, 
                                   nrow = ifelse(
                                     summarytype=="colData", 1, 
                                     nrow(scef))))
      colnames(dfgi) <- c("num.entries", groupstat)
      if(groupi %in% scef[[groupvar]]){
        sceff <- scef[,scef[[groupvar]]==groupi]
        exprff <- assays(sceff)[[assayname]]
        if(summarytype == "colData"){
          exprff <- matrix(colMeans(t(exprff)), nrow = 1)
        }
        dfgi$num.entries <- ncol(exprff)
        if("mean" %in% groupstat){dfgi$mean <- rowMeans(exprff)}
        if("median" %in% groupstat){dfgi$median <- rowMedians(exprff)}
        if("var" %in% groupstat){dfgi$var <- rowVars(exprff)}
        if("sd" %in% groupstat){dfgi$sd <- rowSds(exprff)}
        if("numzero" %in% groupstat){
          dfgi$numzero <- unlist(
            apply(exprff, 1, function(ri){length(which(ri==0))}))
        }
      } else{
        if(verbose){
          message("Skipping summary for group: ", 
                  groupi, "...")}
      }
      cnstr <- ";"
      if(summarytype == "colData"){
        cnstr <- paste0(cnstr, "colData_means;")}
      colnames(dfgi) <- paste0(groupi, cnstr, colnames(dfgi))
      dfgi
    }))
    return(dfg)
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
ugroupv <- unique(sce[[groupvar]]) # get all possible group lvls
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
    dfg <- sce_groupstat(scef = scef, groupvar = groupvar, ugroupv = ugroupv,
                         summarytype = "rowData", assayname = assayname, 
                         verbose = verbose)
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
    dfg <- sce_groupstat(scef = scef, groupvar = groupvar, ugroupv = ugroupv, 
                         assayname = assayname, summarytype = "colData", 
                         verbose = verbose)
    condv <- is(dfg, "data.frame") & nrow(dfg) == nrow(dfr)
    if(condv){
      if(verbose){message("Binding group-level data.")};dfr <- cbind(dfr, dfg)
    }
  }
  return(dfr)
}))


# get summary stats
set <- SummarizedExperimentTypes(assays = list(mexpr),
                                 rowData = rd, colData = cd)



