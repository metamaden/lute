#!/usr/bin/env R

# Author: Sean Maden
#
# Utilities and miscellaneous functions supporting the lute package for deconvolution experiments.
#

#' @importFrom Biobase ExpressionSet AnnotatedDataFrame pData exprs
#' @importFrom SingleCellExperiment counts SingleCellExperiment
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#'
.get_celltypes_from_sce <- function(sce, celltype.variable = "celltype"){
  celltype.vector <- as.data.frame(SummarizedExperiment::colData(sce))[,celltype.variable]
  celltype.char <- as.character(celltype.vector)
  unique.types <- unique(celltype.char)
  unique.types <- unique.types[order(unique.types)]
  celltype.fact <- factor(celltype.vector, levels = unique.types)
  list(variable = celltype.variable, unique.types = unique.types,
       character = celltype.char, factor = celltype.fact)
}

#' ypb_from_sce
#'
#' Get pseudobulk from a SingleCellExperiment object.
#'
#' @param sce An object of type \linkS4class{SingleCellExperiment}.
#' @param assay.name Name of expression matrix in \code{sce} assays.
#' @param celltype.variable Variable name for cell type labels in \code{sce} 
#' coldata.
#' @param S Vector of cell type size scale factors. Optional.
#' @returns Matrix of simulated bulk convoluted signals.
#' @export
ypb_from_sce <- function(sce, assay.name, celltype.variable, S = NULL){
  ltype <- .get_celltypes_from_sce(sce, celltype.variable)
  if(is(S, "NULL")){
    S <- rep(1, length(ltype[["unique.types"]]))
    names(S) <- ltype[["unique.types"]]
  }
  Znew <- .get_z_from_sce(sce, assay.name, celltype.variable)
  P <- prop.table(table(ltype[["character"]]))
  P <- P[order(match(names(P), ltype[["unique.types"]]))]
  ZSnew <- .zstransform(Znew, S)
  ypb <- t(t(P) %*% t(ZSnew))
  return(ypb)
}

#' signature_matrix_from_sce
#' 
#' Calculate a Z signature matrix from object of type 
#' \linkS4class{SingleCellExperiment}.
#' 
#' @param sce An object of type \linkS4class{SingleCellExperiment}.
#' @param assay.name Name of expression matrix in \code{sce} assays.
#' @param celltype.variable Variable name for cell type labels in \code{sce} 
#' coldata. 
#' @param summary.method Summary statistic function to use.
#' @details Calculate a Z signature matrix from object of type 
#' \linkS4class{SingleCellExperiment}.
#' @returns New Z signature matrix.
#' @examples
#' sce.example <- random_sce()
#' signature_matrix_from_sce(sce.example)
#' @export
signature_matrix_from_sce <- function(sce, 
                                      celltype.variable = "celltype", 
                                      summary.method = "mean", 
                                      assay.name = "counts"){
  require(dplyr)
  # gets the z signature matrix from an sce object
  expression.matrix <- assays(sce)[[assay.name]] %>% as.matrix()
  cd <- colData(sce)
  unique.cell.types <- cd[,celltype.variable] %>% unique()
  unique.cell.types <- unique.cell.types[order(unique.cell.types)]
  z <- do.call(cbind, lapply(unique.cell.types, function(cell.type.index){
    filter.index <- cd[,celltype.variable]==cell.type.index
    if(summary.method == "mean"){
      DelayedArray::rowMeans(expression.matrix[,filter.index])
    } else{
      stop('Error, unrecognized summary.method.')}
  }))
  colnames(z) <- unique.cell.types
  return(z)
}

.get_z_from_sce <- function(sce, assay.name = "counts", 
                            celltype.variable = "celltype"){
  ltype <- .get_celltypes_from_sce(sce = sce, celltype.variable = celltype.variable)
  mexpr <- as.matrix(assays(sce)[[assay.name]])
  Znew <- do.call(cbind, lapply(ltype[["unique.types"]], function(typei){
    datav <- mexpr[,ltype[["character"]]==typei]; .z_operator(datav)
  }))
  colnames(Znew) <- ltype[["unique.types"]]
  rownames(Znew) <- rownames(sce)
  return(Znew)
}

.get_sce_from_eset <- function(eset){
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = exprs(eset)))
  SummarizedExperiment::colData(sce) <- DataFrame(Biobase::pData(eset))
  return(sce)
}

.get_eset_from_matrix <- function(mat, batch.variable = "SampleName"){
  pdata <- data.frame(new.variable = colnames(mat))
  colnames(pdata) <- batch.variable
  rownames(pdata) <- colnames(y)
  eset <- Biobase::ExpressionSet(assayData = mat, phenoData = Biobase::AnnotatedDataFrame(pdata))
  return(est)
}

.zstransform <- function(z, s){
  sweep(z, 2, s, FUN = "*")
}

.z_operator <- function(datav){
  rowMeans(datav)
}

.get_decon_example_data <- function(seed.num=0){
  set.seed(seed.num)
  set.seed(0)
  y <- matrix(rnbinom(n=10, size=10, mu=10), ncol = 1)
  z <- matrix(rnbinom(n=20, size=10, mu=10), ncol = 2)
  rownames(y) <- rownames(z) <- paste0("marker", seq(nrow(y)))
  colnames(z) <- paste0("type", seq(ncol(z)))
  colnames(y) <- paste0("sample", seq(ncol(y)))
  s <- c(1, 10)
  return(list(z = z, y = y, s = s))
}

.get_zvar <- function(z){
  z.var <- matrix(0, nrow = nrow(z), ncol = ncol(z))
  rownames(z.var) <- rownames(z)
  colnames(z.var) <- colnames(z)
  z.var
}

.get_decon_example_data_bisque <- function(seed.num = 0){
  set.seed(seed.num)
  # get y.eset
  y <- lute:::.get_decon_example_data()[["y"]]
  y <- cbind(y, y, y, y, y, y)
  colnames(y) <- c(paste0("sample", seq(2)), paste0("bulk",seq(4)))
  df.y.pheno <- data.frame(SubjectName = colnames(y))
  rownames(df.y.pheno) <- colnames(y)
  y.eset <- ExpressionSet(assayData = y, phenoData = AnnotatedDataFrame(df.y.pheno))

  # get z.eset
  sce <- random_sce(num.genes = 10, num.cells = 100, num.types = 2)
  df.z.pheno <- data.frame(cellType = sce[["celltype"]], SubjectName = paste0("sample", seq(ncol(sce))))
  rownames(df.z.pheno) <- colnames(sce)
  z.eset <- ExpressionSet(assayData = counts(sce), phenoData = AnnotatedDataFrame(df.z.pheno))
  rownames(z.eset) <- rownames(y.eset)
  
  # return
  lr <- list(y.eset = y.eset, sc.eset = z.eset)
  return(lr)
}

.get_decon_example_data_scdc <- function(seed.num = 0){
  set.seed(seed.num)
  # get y.eset
  y <- lute:::.get_decon_example_data()[["y"]]
  y <- cbind(y, y, y, y, y, y)
  colnames(y) <- c(paste0("sample", seq(2)), paste0("bulk",seq(4)))
  df.y.pheno <- data.frame(SubjectName = colnames(y))
  rownames(df.y.pheno) <- colnames(y)
  y.eset <- ExpressionSet(assayData = y, phenoData = AnnotatedDataFrame(df.y.pheno))

  # get z.eset
  sce <- random_sce(num.genes = 10, num.cells = 300, num.types = 4)
  df.z.pheno <- data.frame(cellType = sce[["celltype"]], SubjectName = paste0("sample", seq(ncol(sce))))
  rownames(df.z.pheno) <- colnames(sce)
  z.eset <- ExpressionSet(assayData = counts(sce), phenoData = AnnotatedDataFrame(df.z.pheno))
  rownames(z.eset) <- rownames(y.eset)
  
  # return
  lr <- list(y.eset = y.eset, sc.eset = z.eset)
  return(lr)
}

.get_decon_example_data_music2 <- function(seed.num = 0){
  set.seed(seed.num)
  # get y.eset
  y <- lute:::.get_decon_example_data()[["y"]]
  y <- cbind(y, y, y, y, y, y)
  colnames(y) <- c(paste0("sample", seq(2)), paste0("bulk",seq(4)))
  df.y.pheno <- data.frame(SubjectName = colnames(y))
  rownames(df.y.pheno) <- colnames(y)
  y.eset <- ExpressionSet(assayData = y, phenoData = AnnotatedDataFrame(df.y.pheno))

  # get z.eset
  sce <- random_sce(num.genes = 10, num.cells = 300, num.types = 2)
  df.z.pheno <- data.frame(cellType = sce[["celltype"]], SubjectName = paste0("sample", seq(ncol(sce))))
  rownames(df.z.pheno) <- colnames(sce)
  z.eset <- ExpressionSet(assayData = counts(sce), phenoData = AnnotatedDataFrame(df.z.pheno))
  rownames(z.eset) <- rownames(y.eset)
  
  # return
  lr <- list(y.eset = y.eset, sc.eset = z.eset)
  return(lr)
}