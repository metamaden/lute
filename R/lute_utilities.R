#!/usr/bin/env R

### Author: Sean Maden
###
### Utilities and miscellaneous functions supporting the lute package for deconvolution experiments.
###

#' get_celltypes_from_sce
#' 
#' Extract cell type values from SingleCellExperiment.
#' 
#' @param sce A SingleCellExperiment object.
#' @param celltype.variable Variable containing cell type labels (e.g. "type1", 
#' "type2", etc.).
#' @returns List of cell type variable metadata and values.
#'
#' @importFrom Biobase ExpressionSet AnnotatedDataFrame pData exprs
#' @importFrom SingleCellExperiment counts SingleCellExperiment
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#' 
#' @examples
#' lexample <- get_decon_example_data()
#'
#' @export
get_celltypes_from_sce <- function(sce, celltype.variable="celltype"){
  celltype.vector <- as.data.frame(SummarizedExperiment::colData(sce))[,celltype.variable]
  celltype.char <- as.character(celltype.vector)
  unique.types <- unique(celltype.char)
  unique.types <- unique.types[order(unique.types)]
  celltype.fact <- factor(celltype.vector, levels=unique.types)
  lr <- list(variable=celltype.variable, 
             unique.types=unique.types, 
             character=celltype.char, 
             factor=celltype.fact)
  return(lr)
}

#' ypb_from_sce
#'
#' Get pseudobulk from a SingleCellExperiment object.
#'
#' @param sce An object of type \linkS4class{SingleCellExperiment}.
#' @param assay.name Name of expression matrix in \code{sce} assays.
#' @param celltype.variable Variable name for cell type labels in \code{sce} 
#' coldata.
#' @param sample.id.variable Variable name for sample/group ID labels in 
#' \code{sce} coldata.
#' @param S Vector of cell type size scale factors. Optional.
#' @returns Matrix of simulated bulk convoluted signals.
#' 
#' @examples
#' sce.example <- random_sce()
#' ypb_from_sce(sce.example)
#' 
#' @export
ypb_from_sce <- function(sce, assay.name="counts", 
                         celltype.variable="celltype", 
                         sample.id.variable=NULL, S=NULL){
  num.groups <- 1; unique.group.id.vector <- ""
  if(!is(sample.id.variable, "NULL")){
    group.id.vector <- sce[[sample.id.variable]]
    unique.group.id.vector <- group.id.vector 
    unique.group.id.vector <- unique(unique.group.id.vector)
    unique.group.id.vector <- as.character(unique.group.id.vector)
    num.groups <- length(unique.group.id.vector)
  }
  list.cell.types <- get_celltypes_from_sce(
    sce=sce, celltype.variable=celltype.variable)
  num.types <- length(list.cell.types[["unique.types"]])
  ypb.list <- lapply(unique.group.id.vector, function(group.id){
    sce.filter <- sce
    if(num.groups > 1){
      filter.group <- sce[[sample.id.variable]]==group.id
      sce.filter <- sce[,filter.group]
    }
    
    input_s <- S
    if(is(input_s, "NULL")){
      input_s <- rep(1, num.types)
      names(input_s) <- list.cell.types[["unique.types"]]
    }
    
    input_Znew <- get_z_from_sce(sce.filter, assay.name, celltype.variable)
    input_P <- table(list.cell.types[["character"]])
    input_P <- prop.table(input_P)
    order.p <- match(names(input_P), list.cell.types[["unique.types"]])
    order.p <- order(order.p)
    input_P <- input_P[order.p]
    input_ZSnew <- .zstransform(input_Znew, input_s)
    input_ypb <- t(t(input_P) %*% t(input_ZSnew))
    return(input_ypb)
  })
  ypb.table <- do.call(cbind, ypb.list)
  ypb.table <- as.data.frame(ypb.table)
  if(num.groups > 1){
    colnames(ypb.table) <- unique.group.id.vector
  } else{
    colnames(ypb.table) <- "sce.pseudobulk"
  }
  return(ypb.table)
}

#' signature_matrix_from_sce
#' 
#' Calculate a Z signature matrix from object of type 
#' \linkS4class{SingleCellExperiment}.
#' 
#' @param sce An object of type \linkS4class{SingleCellExperiment}.
#' @param assay.name Name of expression matrix in \code{sce} assays (e.g. 
#' "counts").
#' @param celltype.variable Variable name for cell type labels in \code{sce} 
#' coldata (e.g. "type1", "type2", etc.). 
#' @param summary.method Summary statistic function to use.
#' @details Calculate a Z signature matrix from object of type 
#' \linkS4class{SingleCellExperiment}.
#' @returns New Z signature matrix.
#' 
#' @examples
#' sce.example <- random_sce()
#' signature_matrix_from_sce(sce.example)
#' 
#' @export
signature_matrix_from_sce <- function(sce, 
                                      celltype.variable="celltype", 
                                      summary.method="mean", 
                                      assay.name="counts"){
  ## gets the z signature matrix from an sce object
  expression.matrix <- assays(sce)[[assay.name]]
  expression.matrix <- as.matrix(expression.matrix)
  cd <- colData(sce)
  unique.cell.types <- unique(cd[,celltype.variable])
  unique.cell.types <- unique.cell.types[order(unique.cell.types)]
  input_z <- do.call(cbind, lapply(unique.cell.types, function(cell.type.index){
    filter.index <- cd[,celltype.variable]==cell.type.index
    if(summary.method == "mean"){
      DelayedArray::rowMeans(expression.matrix[,filter.index])
    } else{
      Biobase::rowMedians(expression.matrix[,filter.index])
    }
  }))
  colnames(input_z) <- unique.cell.types
  return(input_z)
}

#' get_z_from_sce
#' 
#' Makes the Z cell atlas reference from a SingleCellExperiment.
#' 
#' @param sce A SingleCellExperiment object.
#' @param assay.name Name of expression assay type (e.g. "counts").
#' @param celltype.variable Name of variable containing cell type labels (e.g. 
#' "type1", "type2", etc.).
#' @returns Matrix of cell summary values (Z reference atlas).
#'
#' @importFrom SummarizedExperiment assays
#' 
#' @examples
#' lexample <- get_decon_example_data()
#' @export
get_z_from_sce <- function(sce, assay.name="counts", celltype.variable="celltype"){
  ltype <- get_celltypes_from_sce(sce=sce, celltype.variable=celltype.variable)
  mexpr <- as.matrix(assays(sce)[[assay.name]])
  input_Znew <- do.call(cbind, lapply(ltype[["unique.types"]], function(typei){
    datav <- mexpr[,ltype[["character"]]==typei]; .z_operator(datav)
  }))
  colnames(input_Znew) <- ltype[["unique.types"]]
  rownames(input_Znew) <- rownames(sce)
  return(input_Znew)
}

.zstransform <- function(z, s){
  sweep(z, 2, s, FUN="*")
}

.z_operator <- function(datav){
  rowMeans(datav)
}

#' get_decon_example_data
#' 
#' Make example data for deconvolution.
#' 
#' @param num.bulk.samples Number of bulk samples.
#' @param num.markers Number of cell type markers.
#' @param num.types Number of cell types.
#' @returns Example data as list.
#' 
#' @importFrom stats rpois
#' 
#' @examples
#' example.data <- get_decon_example_data()
#' 
#' @export
get_decon_example_data <- function(num.bulk.samples=2, num.markers=10,
                                    num.types=2){
  input_y <- matrix(
    rpois(n=num.markers*num.bulk.samples, lambda=seq(0, 50, 5)), 
    ncol=num.bulk.samples)
  input_z <- matrix(
    rpois(n=num.types*num.markers, lambda=seq(0, 50, 5)), 
    ncol=num.types)
  rownames(input_y) <- rownames(input_z) <- paste0("marker", seq(num.markers))
  colnames(input_z) <- paste0("type", seq(num.types))
  colnames(input_y) <- paste0("sample", seq(num.bulk.samples))
  input_s <- c(1, 10)
  names(input_s) <- colnames(input_z)
  return(list(z=input_z, y=input_y, s=input_s))
}

.get_zvar <- function(z){
  input_z_variable <- matrix(0, nrow=nrow(z), ncol=ncol(z))
  rownames(input_z_variable) <- rownames(z)
  colnames(input_z_variable) <- colnames(z)
  input_z_variable
}

#' get_decon_example_data_bisque
#'
#' Get example data for Bisque algorithm.
#'
#' @param num.bulk.samples Number of bulk samples.
#' @param num.markers Number of cell type markers.
#' @param num.cells Number of cells.
#' @param num.types Number of cell types.
#' @returns Example data as list.
#'
#' @importFrom Biobase ExpressionSet
#' @importFrom Biobase AnnotatedDataFrame
#' @importFrom BiocGenerics counts
#' 
#' @examples
#' example.data <- get_decon_example_data()
#'
#' @export
get_decon_example_data_bisque <- function(num.bulk.samples=100,
                                           num.markers=1000, 
                                           num.cells=1000, 
                                           num.types=2){
  lexample <- get_decon_example_data(num.bulk.samples=num.bulk.samples,
                               num.markers=num.markers,
                               num.types=num.types)
  input_y <- lexample[["y"]]
  colnames(input_y) <- c(paste0("sample", seq(num.bulk.samples/2)), 
                   paste0("bulk", seq(num.bulk.samples/2)))
  df.y.pheno <- data.frame(SubjectName=colnames(input_y))
  rownames(df.y.pheno) <- colnames(input_y)
  y.eset <- ExpressionSet(assayData=input_y, 
                          phenoData=AnnotatedDataFrame(df.y.pheno))
  sce <- random_sce(num.genes=num.markers, 
                    num.cells=num.cells, 
                    num.types=num.types)
  df.z.pheno <- data.frame(cellType=sce[["celltype"]], 
                           SubjectName=
                             paste0("sample", seq(num.cells)))
  rownames(df.z.pheno) <- colnames(sce)
  z.eset <- ExpressionSet(assayData=counts(sce), 
                          phenoData=AnnotatedDataFrame(df.z.pheno))
  rownames(z.eset) <- rownames(y.eset)
  return_list <- list(y.eset=y.eset, sc.eset=z.eset)
  return(return_list)
}

#' get_decon_example_data_scdc
#' 
#' Get example data for SCDC
#'
#' @returns Example data as list.
#'
#' @importFrom Biobase ExpressionSet
#' @importFrom Biobase AnnotatedDataFrame
#' @importFrom BiocGenerics counts
#' 
#' @examples
#' example.data <- get_decon_example_data()
#' 
#' @export
get_decon_example_data_scdc <- function(){
  ## get y.eset
  input_y <- get_decon_example_data()[["y"]]
  input_y <- cbind(input_y, input_y, input_y, input_y, input_y, input_y)
  colnames(input_y) <- c(paste0("sample", seq(2)), paste0("bulk",seq(4)))
  df.y.pheno <- data.frame(SubjectName=colnames(input_y))
  rownames(df.y.pheno) <- colnames(input_y)
  y.eset <- ExpressionSet(assayData=input_y, 
                          phenoData=AnnotatedDataFrame(df.y.pheno))

  ## get z.eset
  sce <- random_sce(num.genes=10, num.cells=300, num.types=4)
  df.z.pheno <- data.frame(cellType=sce[["celltype"]], 
                           SubjectName=paste0("sample", seq(ncol(sce))))
  rownames(df.z.pheno) <- colnames(sce)
  z.eset <- ExpressionSet(assayData=counts(sce), 
                          phenoData=AnnotatedDataFrame(df.z.pheno))
  rownames(z.eset) <- rownames(y.eset)
  
  ## return
  return_list <- list(y.eset=y.eset, sc.eset=z.eset)
  return(return_list)
}

#'
#'
#' @importFrom Biobase ExpressionSet
#' @importFrom Biobase AnnotatedDataFrame
#' @importFrom BiocGenerics counts
#'
get_decon_example_data_music2 <- function(){
  ## get y.eset
  input_y <- get_decon_example_data()[["y"]]
  input_y <- cbind(input_y, input_y, input_y, input_y, input_y, input_y)
  colnames(input_y) <- c(paste0("sample", seq(2)), 
                   paste0("bulk",seq(ncol(input_y)-2)))
  
  df.y.pheno <- data.frame(SubjectName=colnames(input_y))
  rownames(df.y.pheno) <- colnames(input_y)
  y.eset <- ExpressionSet(assayData=input_y, 
                          phenoData=AnnotatedDataFrame(df.y.pheno))

  ## get z.eset
  sce <- random_sce(num.genes=10, num.cells=300, num.types=2)
  df.z.pheno <- data.frame(cellType=sce[["celltype"]], 
                           SubjectName=paste0("sample", seq(ncol(sce))))
  rownames(df.z.pheno) <- colnames(sce)
  z.eset <- ExpressionSet(assayData=counts(sce), 
                          phenoData=AnnotatedDataFrame(df.z.pheno))
  rownames(z.eset) <- rownames(y.eset)
  
  ## return
  return_list <- list(y.eset=y.eset, sc.eset=z.eset)
  return(return_list)
}


#'
#'
#' @param list.pred List of cell type proportions predictions.
#' @param column.labels Vector of cell type labels 
#' (e.g. "type1", "type2", etc.).
#' @param row.labels Vector of sample id labels 
#' (e.g. "sample1", "sample2", etc.).
#'
#'
.parse_deconvolution_predictions_results <- function(list.pred, 
                                                     column.labels, 
                                                     row.labels){
  if(is(column.labels, "NULL")){column.labels <- seq(length(list.pred[[1]]))}
  if(is(row.labels, "NULL")){row.labels <- seq(length(list.pred))}
  table.pred <- do.call(rbind, list.pred)
  table.pred <- apply(table.pred, 1, function(ri){ri/sum(ri)}) 
  table.pred <- t(table.pred)
  colnames(table.pred) <- column.labels
  rownames(table.pred) <- row.labels
  ## convert
  table.pred <- cellProportionsPredictions(table.pred, column.labels, 
                                           row.labels)
  return(table.pred)
}
