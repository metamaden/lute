libv <- c("SingleCellExperiment", "SummarizedExperiment")
sapply(libv, library, character.only = T)

sce.fpath <- file.path("data", "SCE_DLPFC-n3_tran-etal.rda")
sce <- get(load(sce.fpath))

expr.scale <- "logcounts"
marker.varname <- "marker"
type.varname <- "cellType"
marker.summary <- "mean"
return.type <- "lgv"

# check sce type
cond.sce <- is(sce, "SingleCellExperiment")
if(!cond.sce){stop("Error, sce not a SingleCellExperiment.")}
# check marker variable
cond.marker <- marker.varname %in% colnames(rowData(sce))
scef <- sce
if(cond.marker){
  which.marker <- which(rowData(scef)[,marker.varname])
  scef <- scef[which.marker,]
} else{
  if(verbose){message("marker.varname not in rowData; using all genes.")}
}
# check scale
cond.scale <- expr.scale %in% names(assays(scef))
if(!cond.scale){stop("Error, cond.scale value not found in sce assays.")}
me <- eval(parse(text = paste0("assays(scef)$",expr.scale)))
# check type variable
cd <- colData(scef); cond.type <- type.varname %in% colnames(cd)
if(!cond.type){stop("Error, type.varname not in colData.")}
typevar <- cd[,type.varname]; typev <- unique(typevar)
# parse marker summary options
if(marker.summary == "mean"){
  if(verbose){message("Getting marker mean expression by type...")}
  z <- do.call(cbind, lapply(typev, function(typei){
    which.typei <- typevar==typei
    DelayedArray::rowMeans(me[,which.typei])
  }))
} else if(marker.summary == "median"){
  if(verbose){message("Getting marker median expression by type...")}
  z <- do.call(cbind, lapply(typev, function(typei){
    which.typei <- typevar==typei
    matrixStats::rowMedians(me[,which.typei])
  }))
} else{
  stop("Error, invalid marker.summary option.")
}
colnames(z) <- typev
if(return.type == "lgv"){
  if(verbose){message("Returning lgv...")}
  lgv <- list(lapply(seq(ncol(z)), function(ci){z[,ci]}))
  names(lgv) <- colnames(z)
  return(lgv)
}
if(verbose){message("Returning Z...")}
return(z)


kexpr_sce <- function(sce, return.lgv = TRUE, expr.scale = "logcounts", 
                      marker.varname = "marker", type.varname = "cellType", 
                      marker.summary = "mean", verbose = FALSE){
  # check sce type
  cond.sce <- is(sce, "SingleCellExperiment")
  if(!cond.sce){stop("Error, sce not a SingleCellExperiment.")}
  # check marker variable
  cond.marker <- marker.varname %in% colnames(rowData(sce)); scef <- sce
  if(cond.marker){
    which.marker <- which(rowData(scef)[,marker.varname])
    scef <- scef[which.marker,]
  } else{
    if(verbose){message("marker.varname not in rowData; using all genes.")}
  }
  # check scale
  cond.scale <- expr.scale %in% names(assays(scef))
  if(!cond.scale){stop("Error, cond.scale not in sce assays.")}
  me <- eval(parse(text = paste0("assays(scef)$",expr.scale)))
  # check type variable
  cd <- colData(scef); cond.type <- type.varname %in% colnames(cd)
  if(!cond.type){stop("Error, type.varname not in colData.")}
  typevar <- cd[,type.varname]; typev <- unique(typevar)
  # parse marker summary options
  if(marker.summary == "mean"){
    if(verbose){message("Getting marker mean expression by type...")}
    z <- do.call(cbind, lapply(typev, function(typei){
      DelayedArray::rowMeans(me[,typevar==typei])
    }))
  } else if(marker.summary == "median"){
    if(verbose){message("Getting marker median expression by type...")}
    z <- do.call(cbind, lapply(typev, function(typei){
      matrixStats::rowMedians(me[,typevar==typei])
    }))
  } else{
    stop("Error, invalid marker.summary option.")
  }
  colnames(z) <- typev
  if(return.lgv){
    if(verbose){message("Returning lgv...")}
    lgv <- lapply(seq(ncol(z)), function(ci){z[,ci]})
    names(lgv) <- colnames(z)
    return(lgv)
  }
  if(verbose){message("Returning Z reference...")}
  return(z)
}


lgv <- kexpr_sce(sce, verbose = T)




# test implementation
library(lute)


sce.fpath <- file.path("data", "SCE_DLPFC-n3_tran-etal.rda")
sce <- get(load(sce.fpath))

lpv <- list(c(0.2, 0.8))
lsv <- list(c(1,10))

lda <- decon_analysis(lpv = lpv, lsv = lsv, sce = sce, type.varname = "k2", 
                      verbose = T)

lda <- decon_analysis(lpv = lpv, sce = sce, type.varname = "k2", 
                      verbose = T)

lda <- decon_analysis(lpv = lpv, lsv = lsv, sce = sce, type.varname = "cellType", 
                      verbose = T)


lda$lgg$ggpt1
lda$lgg$ggvp
lda$lgg$ggpt2








