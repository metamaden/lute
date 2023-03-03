.get_celltypes_from_sce <- function(sce, celltype.variable = "celltype"){
  require(SingleCellExperiment); require(SummarizedExperiment)
  celltype.vector <- sce[[celltype.variable]]
  celltype.char <- as.character(celltype.vector)
  unique.types <- unique(celltype.char)
  unique.types <- unique.types[order(unique.types)]
  celltype.fact <- factor(celltype.vector, levels = unique.types)
  list(variable = celltype.variable,
       unique.types = unique.types,
       character = celltype.char,
       factor = celltype.fact)
}

.get_z_from_sce <- function(sce, assay.name = "counts", 
                            celltype.variable = "celltype"){
  require(SingleCellExperiment); require(SummarizedExperiment)
  ltype <- .get_celltypes_from_sce()
  mexpr <- as.matrix(assays(sce)[[assay.name]])
  Znew <- do.call(cbind, lapply(ltype[["unique.types"]], function(typei){
    datav <- mexpr[,ltype[["character"]]==typei]; .z_operator(datav)
  }))
  colnames(Znew) <- ltype[["unique.types"]]
  rownames(Znew) <- rownames(sce)
  return(Znew)
}

.ypb_from_sce <- function(sce, assay.name, celltype.variable, S = NULL){
  ltype <- .get_celltypes_from_sce()
  if(is(S, "NULL")){
    S <- rep(1, length(ltype[["unique.types"]]))
    names(S) <- ltype[["unique.types"]]
  }
  Znew <- .get_z_from_sce(sce, assay.name, celltype.variable)
  P <- prop.table(table(ltype[["character"]]))
  P <- P[order(match(names(P), ltype[["unique.types"]]))]
  ZSnew <- sweep(Znew, 2, S, "*")
  ypb <- t(t(P) %*% t(ZSnew))
  return(ypb)
}

.z_operator <- function(datav){
  rowMeans(datav)
}


