#!/usr/bin/env R

# Author: Sean Maden
#
# Methods to select markers using typemarkers() and supported methods.
#

# this is the main task for parallel selection of group markers
.get.group.markers.iter <- function(task.iter, sce, assay.name, 
                                   celltype.variable, markers.per.type,
                                   downsample.within.group, verbose = F){
  
  if(verbose){message("Getting markers for group id: ", task.iter, "...")}
  filter <- sce[[group.variable]] == task.iter; sce.filter <- sce[,filter]
  expression.matrix <- assays(sce.filter)[[assay.name]]
  
  # downsampling
  if(downsample.within.group){
    batch.variable <- sce.filter[[celltype.variable]]
    expression.matrix <- downsampleBatches(expression.matrix,
                                           batch = batch.variable) %>% 
      as.matrix()
  }
  assays(sce.filter)[[assay.name]] <- expression.matrix
  
  
  # results
  group.marker.results <- meanratiosParam(sce = sce.filter, return.info = TRUE, 
                                          assay.name = assay.name,
                                          celltype.variable = celltype.variable, 
                                          markers.per.type = markers.per.type) %>% 
    typemarkers()
  group.marker.table <- group.marker.results$result.info
  if(nrow(group.marker.table) > 0){
    group.marker.table <- group.marker.table[
      group.marker.table$rank_ratio >= markers.per.type,]
    topmarkers.vector <- group.marker.table$markers
    filter.topmarkers <- group.marker.table$gene %in% topmarkers.vector
    group.marker.table <- group.marker.table[filter.topmarkers,]
    group.marker.table$group.id <- task.iter
  }
  return(group.marker.table)
}

# this gets the markers by group, returning a list of markers by group

#' markers_by_group
#'
#' Get gene markers organized by groups.
#' 
#' @param sce Object of type \linkS4class{SingleCellExperiment}.
#' 
#' @param group.variable Variable name in \code{sce} coldata containing the 
#' group ID labels.
#' 
#' @param celltype.variable Variable name in \code{sce} coldata containing the 
#' cell type labels.
#' 
#' @param assay.name Name of the expression matrix in \code{sce} assays to use 
#' for marker selection.
#' 
#' @param markers.per.type Target maximum number of genes per cell type to select.
#' 
#' @param typemarker.algorithm Algorithm for typemarker selection, supported in 
#' \code{lute} (see \code{?typemarkers} for details).
#' 
#' @param return.type Format of return data, either "list" or "table".
#' 
#' @param downsample.within.group Whether to downsample on cell types prior to 
#' marker selection (performed *within* groups)
#' 
#' @param parallelize Whether to use parallelization (mclapply) when getting 
#' markers by group. Otherwise uses lapply.
#' 
#' @param verbose Whether to show verbose status messages.
#' 
#' @returns List or table of markers (see \code{return.type} argument), 
#' organized by group IDs.
#' 
#' @details Gets the within-group cell type markers.
#' 
#' @examples
#' sce.example <- random_sce(num.cells = 100, num.genes = 100)
#' sce.example[["sample_id"]] <- c(rep("sample1", 10), rep("sample2", 80), rep("sample1", 10))
#' group.markers <- markers_by_group(sce.example)
#' 
#' 
#' 
#' @seealso \code{filter_group_markers}
#' @export
markers_by_group <- function(sce, 
                             group.variable = "sample_id", 
                             celltype.variable = "celltype", 
                             assay.name = "counts", 
                             markers.per.type = 100, 
                             typemarker.algorithm = "meanratios",
                             return.type = "list",
                             downsample.within.group = T,
                             parallelize = T,
                             verbose = FALSE){
  require(dplyr)
  require(parallel)
  require(scuttle)
  
  # the goal is ONLY to annotate concordance and overlap
  sce = sce.example
  group.variable = "sample_id"
  celltype.variable = "celltype"
  assay.name = "counts"
  markers.per.type = 20
  typemarker.algorithm = "meanratios"
  return.type = "list"
  verbose = FALSE
  downsample.within.group = T
  parallelize = T
  
  if(verbose){
    message("Getting gene markers for each specified group")}
  unique.group.vector <- sce[[group.variable]] %>% unique() %>% as.character()
  
  
  # parse task with parallel option
  if(parallelize){
    group.markers.list <- mclapply(unique.group.vector, 
                                   lute:::.get.group.markers.iter, 
                                   sce.example,
                                   assay.name, 
                                   celltype.variable, 
                                   markers.per.type)
  } else{
    group.markers.list <- lapply(unique.group.vector, 
                                 get.group.markers.iter, 
                                 sce.example, 
                                 assay.name, 
                                 celltype.variable, 
                                 markers.per.type)
  }
  names(group.markers.list) <- unique.group.vector
  if(return.type == "list"){
    return.data <- group.markers.list
    names(return.data)
  } else{
    return.data <- do.call(rbind, group.markers.list) %>% as.data.frame()
  }
  if(verbose){message("finished marker lists. Returning...")}
}
