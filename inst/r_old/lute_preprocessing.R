#!/usr/bin/env R

# Author: Sean Maden
#
# Preprocessing functions.

#' sce_preprocess_groups
#'
#' Perform group adjustment on a SingleCellExperiment.
#'
#' @param group.variable Variable containing group labels ("sample.id").
#' @param celltype.variable Variable containing cell type labels ("celltype").
#' @param assay.name Name of assay to adjust ("counts").
#' @param new.assay.stem New assay character string stem ("adjusted").
#' @param downsample Whether to perform downsampling on group-type cell amounts (TRUE).
#' @param zscale.threshold Minimum threshold to filter whole groups on their 
#' z-scale values for group-type cell amounts (NULL to skip).
#' @param negative.to.zero.expression Whether to set negative expression values to 
#' zero after adjustment (TRUE).
#' @param return.type Either "list", including summary statistics metadata, or "sce" for
#'  SingleCellExperiment without metadata.
#' @param verbose Whether to show verbose status messages.
#' 
#' @returns SingleCellExperiment object with group-adjusted expression. 
#' 
#' @details Perform group adjustment on a SingleCellExperiment object. Includes options to
#' perform downsampling and filtering of whole groups on the z-scale of group-type cell 
#' amounts.
#' 
#' @examples 
#' #sce.example <- random_sce(num.genes = 1000, num.cells = 1000)
#' #sce.example[["sample.id"]] <- c(rep("sample1", 800), rep("sample2", 200))
#' #list.preprocess <- sce_preprocess_groups(sce.example)
#' #names(list.preprocess)
#' #names(list.preprocess$metadata)
#' #names(assays(sce))
#' 
sce_preprocess_groups <- function(sce, group.variable = "sample.id", 
                             celltype.variable = "celltype",
                             assay.name = "counts", new.assay.stem = "adjusted",
                             downsample = TRUE, zscale.threshold = NULL, 
                             negative.to.zero.expression = TRUE,
                             return.type = "list", verbose = FALSE){
  require(dplyr)
  coldata <- colData(sce)
  if(!group.variable %in% colnames(coldata)){
    stop("group.variable not in coldata.")}
  if(!celltype.variable %in% colnames(coldata)){
    stop("celltype.variable not in coldata.")}
  message("sce object contains ", ncol(sce), " cells.")
  if(return.type == "list"){
    list.metadata <- list()
    list.metadata[["zero;input"]] <- 
      .sce_summarize_groups(sce, assay.name, group.variable)
  }
  if(downsample){
    if(verbose){message("Performing downsampling...")}
    sce <- .downsample_sce(sce, zscale.threshold, assay.name, 
                           celltype.variable, group.variable)
    if(return.type == "list"){
      list.metadata[["first;downsample"]] <- 
        .sce_summarize_groups(sce, assay.name, group.variable)}
  }
  if(!is(zscale.threshold, "NULL")){
    if(verbose){message("Applying scale filter...")}
    sce <- .zscale_filter(sce, zscale.threshold, group.variable, celltype.variable)
    if(return.type == "list"){
      list.metadata[["second;zscale.filter"]] <- 
        .sce_summarize_groups(sce, assay.name, group.variable)}
    message("After z scale threshold filter, sce object contains ", 
            ncol(sce), " cells.")
  }
  if(verbose){message("Adjusting on groups...")}
  sce <- .sce_group_adjustment(sce, assay.name, new.assay.stem, group.variable, 
                                 celltype.variable, negative.to.zero.expression)
  if(return.type == "list"){
    list.metadata[["third;group.adjustment"]] <- 
      .sce_summarize_groups(sce, assay.name, group.variable)
    list.metadata[["data.dispersion"]] <- .data_dispersion(list.metadata)
    list.metadata[["plots.dispersion"]] <- 
      .plot_dispersion(list.metadata$data.dispersion)
    return.object <- list(sce = sce, metadata = list.metadata)
  } else{
    return.object <- sce
  }
  return(return.object)
}

#'
#'
#'
.downsample_sce <- function(sce, zscale.threshold, assay.name, 
                            celltype.variable, group.variable){
  require(dplyr); require(scuttle)
  group.vector <- sce[[group.variable]]
  unique.groups <- group.vector %>% unique()
  cell.types.vector <- sce[[celltype.variable]]
  unique.cell.types <- cell.types.vector %>% unique()
  sce <- do.call(cbind, lapply(unique.cell.types, function(type){
    message("Downsampling for type ", type, "...")
    filter.sce <- cell.types.vector==type
    sce.type <- sce[,filter.sce]
    group.vector <- sce.type[[group.variable]]
    expression.matrix <- assays(sce.type)[[assay.name]]
    expression.matrix.downsampled <- 
      downsampleBatches(expression.matrix, batch = group.vector)
    assays(sce.type)[[assay.name]] <- expression.matrix.downsampled 
    sce.type
  }))
  return(sce)
}

#'
#'
#'
.zscale_filter <- function(sce, zscale.threshold, group.variable, celltype.variable){
  require(dplyr)
  coldata <- colData(sce)
  cell.counts.table <- table(coldata[,group.variable], 
                             coldata[,celltype.variable]) %>% 
    as.data.frame()
  scale.vector <- scale(cell.counts.table$Freq)[,1]
  filter.counts.table <- which(scale.vector <= zscale.threshold)
  filter.groups.on.scale <- cell.counts.table[filter.counts.table, 1] %>% 
    as.character()
  filter.groups <- coldata[,group.variable] %in% filter.groups.on.scale
  sce <- sce[,!filter.groups]
  return(sce)
}

#'
#'
#'
.sce_group_adjustment <- function(sce, assay.name, new.assay.stem, 
                                    group.variable, celltype.variable,
                                    negative.to.zero.expression){
  require(sva)
  new.assay.name <- paste0(assay.name, "_", new.assay.stem)
  expression.matrix <- assays(sce)[[assay.name]]
  cell.id.vector <- colnames(sce)
  pheno <- data.frame(group = sce[[group.variable]], 
                      celltype = sce[[celltype.variable]])
  mod <- model.matrix(~celltype, data = pheno)
  expression.matrix.adjusted <- ComBat(dat = expression.matrix,
                                       batch = pheno$group, mod = mod)
  if(negative.to.zero.expression){
    filter.negative.expression <- expression.matrix.adjusted < 0
    message("Converting negative values to zero...")
    expression.matrix.adjusted[filter.negative.expression] <- 0
  }
  assays(sce)[[new.assay.name]] <- expression.matrix.adjusted
  return(sce)
}

#'
#'
#'
.sce_summarize_groups <- function(sce, assay.name, group.variable,
                                 statistics.vector = c("mean", "sum", 
                                                       "num.detected", 
                                                       "prop.detected", 
                                                       "median", "var")
                                 ){
  require(dplyr); require(scuttle)
  expression.matrix <- assays(sce)[[assay.name]]
  unique.groups <- sce[[group.variable]] %>% unique() %>% as.character()
  unique.groups <- c(unique.groups, "all")
  list.se.summaries <- lapply(statistics.vector, function(method){
    if(method == "var"){
      list.vars <- lapply(unique.groups, function(group){
        if(group == "all"){
          rowVars(expression.matrix) %>% as.matrix()
        } else{
          group.filter <- sce[[group.variable]] == group
          rowVars(expression.matrix[,group.filter]) %>% as.matrix()
        }
      })
      matrix.vars <- do.call(cbind, list.vars)
      colnames(matrix.vars) <- unique.groups
      rownames(matrix.vars) <- rownames(sce)
      matrix.vars
    } else{ 
      se.group <- summarizeAssayByGroup(expression.matrix, 
                                  sce[[group.variable]],
                                  statistics = method)
      se.all <- summarizeAssayByGroup(expression.matrix, 
                                          rep("all", ncol(sce)),
                                          statistics = method)
      matrix.group <- assays(se.group)[[method]] %>% as.matrix()
      matrix.all <- assays(se.all)[[method]] %>% as.matrix()
      cbind(matrix.group, matrix.all)
    }
    
  })
  names(list.se.summaries) <- statistics.vector
  return(list.se.summaries)
}

#'
#'
#'
.data_dispersion <- function(metadata){
  require(dplyr)
  steps.vector <- names(metadata)
  list.dispersion <- lapply(steps.vector, function(step){
    step.data <- metadata[[step]]
    mean.data <- step.data$mean
    colnames(mean.data) <- paste0(colnames(mean.data), "_mean")
    var.data <- step.data$var
    colnames(var.data) <- paste0(colnames(var.data), "_var")
    data.dispersion <- cbind(mean.data, var.data) %>% as.data.frame()
    data.dispersion$step <- step
    return(data.dispersion)
  })
  data.dispersion <- do.call(rbind, list.dispersion)
  return(data.dispersion)
}

#'
#'
.plot_dispersion <- function(data.dispersion){
  require(dplyr); require(ggplot2)
  lgg <- list()
  data.dispersion$step <- factor(data.dispersion$step, levels = 
                                   c("zero;input", "first;downsample", 
                                     "second;zscale.filter", 
                                     "third;group.adjustment"))
  groups.vector <- gsub("_.*", "", colnames(data.dispersion))
  unique.groups <- groups.vector %>% unique()
  unique.groups <- unique.groups[!unique.groups == "step"]
  lgg <- lapply(unique.groups, function(group){
    filter.data <- grepl(paste0(group, "_"), colnames(data.dispersion))
    plot.data <- data.dispersion[,filter.data]
    plot.data <- plot.data + 1
    colnames(plot.data) <- gsub(".*_", "", colnames(plot.data))
    plot.data <- cbind(plot.data, data.dispersion[,"step",drop=F]) %>% as.data.frame()
    ggplot(plot.data, aes(x = mean, y = var)) + 
      geom_point(alpha = 0.4) + geom_abline(slope = 1, intercept = 0) + 
      geom_smooth(method = "glm") + theme_bw() +
      scale_y_log10() + scale_x_log10() + facet_wrap(~step) +
      xlab("Mean (log10 + 1)") + ylab("Variance (log10 + 1)")
  })
  names(lgg) <- unique.groups
  # plot composite using tall data
  data.tall <- lapply(unique.groups[!unique.groups=="all"], function(group){
    filter.group <- grepl(paste0(group,"_"), colnames(data.dispersion))
    data.filter <- data.dispersion[,filter.group] %>% as.data.frame()
    data.filter <- data.filter + 1
    colnames(data.filter) <- gsub(".*_", "", colnames(data.filter))
    data.filter$group.id <- group
    data.filter$step <- data.dispersion$step
    data.filter
  })
  data.tall <- do.call(rbind, data.tall)
  lgg[["all.composite"]] <- ggplot(data.tall, aes(x = mean, y = var, 
                                 color = group.id, shape = group.id)) + 
    geom_point(alpha = 0.4) + 
    geom_abline(slope = 1, intercept = 0) +
    geom_smooth(method = "glm") + theme_bw() +
    scale_x_log10() + scale_y_log10() + 
    facet_wrap(~step)
  lgg[["data.tall"]] <- data.tall
  return(lgg)
}
