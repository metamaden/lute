#!/usr/bin/env R

# Author: Sean Maden
#
# Methods to select markers using typemarkers() and supported methods.
#

# this is the main task for parallel selection of group markers
get.group.markers.iter <- function(task.iter, sce, assay.name, 
                                   celltype.variable, markers.per.type,
                                   downsample.within.group){
  
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
#' @param group.variable Variable name in \code{sce} coldata containing the 
#' group ID labels.
#' @param celltype.variable Variable name in \code{sce} coldata containing the 
#' cell type labels.
#' @param assay.name Name of the expression matrix in \code{sce} assays to use 
#' for marker selection.
#' @param markers.per.type Target maximum number of genes per cell type to select.
#' @param typemarker.algorithm Algorithm for typemarker selection, supported in 
#' \code{lute} (see \code{?typemarkers} for details).
#' @param return.type Format of return data, either "list" or "table".
#' @param downsample.within.group Whether to downsample on cell types prior to 
#' marker selection (performed *within* groups)
#' @param parallelize Whether to use parallelization (mclapply) when getting 
#' markers by group. Otherwise uses lapply.
#' @param verbose Whether to show verbose status messages.
#' @returns List or table of markers (see \code{return.type} argument), 
#' organized by group IDs.
#' 
#' @details Gets the within-group cell type markers.
#' 
#' @examples
#' sce.example <- random_sce(num.cells = 100, num.genes = 100)
#' sce.example[["sample_id"]] <- c(rep("sample1", 10), rep("sample2", 80), rep("sample1", 10))
#' group.markers <- markers_by_group(sce.example)
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
                                   get.group.markers.iter, 
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















#'
#'
#'
sce_marker_filter <- function(group.markers,
                        minimum.group.overlap.rate = 0.5,
                        num.markers.final = 20,
                        group.adjust = TRUE,
                        verbose = TRUE){
  message("Getting markers to exclude..")
  list.markers.filter.info <- get_filter_group_markers(markers.by.group,
                                                minimum.group.overlap.rate = 0.5,
                                                verbose = TRUE)
  
  
  dim(sce)
  filter.sce <- rownames(sce) %in% list.markers.filter.info$markers.filter.vector
  sce.filter <- sce[!filter.sce,]
  dim(sce.filter)
  metadata(sce.filter)$list.markers.filter.info <- list.markers.filter.info
  
  sce.final <- sce[!filter.sce,]
  if(group.adjust){
    message("Performing group adjustments...")
    sce.final <- sce_preprocess_groups(sce.filter, group.variable = "sample.id", 
                                       celltype.variable = "celltype",
                                       ...)
  }
  
  message("Getting final markers")
  lute(sce.final)
  
  return(sce.filter)
}

#' get_filter_group_markers
#'
#' Get concordant, overlapping markers from a list of markers by group.
#'
#' @param markers.list List of gene markers, organized by groups studied (e.g. 
#' returned from calling \code{marker_overlaps()} (see \code{?marker_overlaps} 
#' for details).
#' @param minimum.group.overlap.rate Minimum amount of overlap by group. Should be a 
#' fraction of the total groups in the data (e.g. 0.2 means the marker is 
#' overlapping in 20% of groups).
#' @param verbose Whether to show verbose status messages.
#' @returns List containing markers for removal/filtering, with metadata about the
#' different filter parameters and type marker selection results.
#' @details Performs the following marker filter steps:
#' * Uniqueness: Removes duplicate markers across type lists.
#' * Concordance: Removes markers which aren't specific to the same type across 
#' groups.
#' * Overlap: Removes markers below an overlap rate, defined using the 
#' \code{min.group.overlap.rate} argument.
#' @seealso \code{markers_by_group}, \code{marker_overlaps}
#' @examples
#' set.seed(0)
#' sce.example <- random_sce(num.cells = 100, num.genes = 1000)
#' sce.example[["sample_id"]] <- c(rep("sample1", 50), rep("sample2", 50))
#' group.markers <- markers_by_group(sce.example, markers.per.type = 5)
#' markers.filter.vector <- get_filter_group_markers(group.markers)
#' length(markers.filter.vector)
#' 
#' @export
get_filter_group_markers <- function(group.markers,
                                     minimum.group.overlap.rate = 0.5,
                                     verbose = TRUE){
  require(dplyr)
  if(is(group.markers, "list")){
    markers.table <- do.call(rbind, group.markers) %>% as.data.frame()
  }
  num.markers.input <- markers.table$gene %>% unique() %>% length()
  if(verbose){message("Found ", num.markers.input, " total markers.")}
  
  if(verbose){message("Getting marker overlap info...")}
  overlap.info <- marker_overlaps(group.markers)
  
  markers.start.vector <- markers.table$gene
  which.non.concordant <- grepl(";", overlap.info$types)
  markers.concordant <- overlap.info[!which.non.concordant,]
  num.concordant <- nrow(markers.concordant)
  if(verbose){message("Found ", num.concordant, " concordant markers by type.")}
  filter.min.overlap <- markers.concordant$overlap.group.rate >= 
    minimum.group.overlap.rate
  markers.overlap <- markers.concordant[filter.min.overlap,]
  num.overlap <- nrow(markers.overlap)
  if(verbose){
    message("Found ", num.overlap, 
            " concordant markers with overlap rate at least ",
            minimum.group.overlap.rate,".")}
  
  # get markers for removal
  markers.non.concordant <- overlap.info[which.non.concordant,]
  markers.low.overlap <- overlap.info[!filter.min.overlap,]
  markers.filter.vector <- unique(c(markers.non.concordant$marker, 
                                  markers.low.overlap$marker))
  # return.list
  metadata.list <- list(markers.low.overlap = markers.low.overlap,
                        markers.non.concordant = markers.non.concordant,
                        num.markers.input = num.markers.input,
                        num.concordant = num.concordant,
                        num.overlap = num.overlap,
                        minimum.group.overlap.rate = minimum.group.overlap.rate)
  return.list <- list(markers.filter.vector = markers.filter.vector,
                      marker.overlap.info = markers.overlap,
                      filter.metadata = metadata.list)
  return(return.list)
}

#' marker_overlaps
#'
#' Gets marker overlap info among groups.
#'
#' @param group.markers Group-specific markers list returned from 
#' \code{markers_by_group} (see \code{?markers_by_group} for details).
#' @param marker.filter.vector Vector of markers to compare (exclude all 
#' additional markers).
#' @returns Table containing marker overlap info.
#' @examples 
#' sce.example <- random_sce(num.cells = 100, num.genes = 100)
#' sce.example[["sample_id"]] <- c(rep("sample1", 20), rep("sample2", 80))
#' group.markers <- markers_by_group(sce.example)
#' group.overlaps <- marker_overlaps(group.markers)
#' 
#' @export
marker_overlaps <- function(group.markers, marker.filter.vector = NULL){
  require(dplyr)
  if(is(group.markers, "list")){
    group.table <- do.call(rbind, group.markers) %>% as.data.frame()
  }
  if(is(group.markers, "matrix")){
    group.table <- group.table %>% as.data.frame()
  }
  if(!is(marker.filter.vector, "NULL")){
    filter.markers <- group.table$gene %in% marker.filter.vector
    group.table <- group.table[filter.markers,]
  }
  # get labels overlapping
  total.groups <- group.table$group.id %>% unique() %>% length()
  unique.markers <- group.table$gene %>% unique()
  overlap.table <- data.frame(marker = unique.markers)
  overlap.table$groups.overlapping <- 
    overlap.table$types.overlapping <- 
    overlap.table$number.overlapping.groups <- 
    overlap.table$number.overlapping.types <- 
    overlap.table$overlap.group.rate <- ""
  for(marker in unique.markers){
    filter.marker <- group.table$gene == marker
    group.table.filter <- group.table[filter.marker,]
    groups.present <- unique(group.table.filter$group.id)
    types.present <- unique(group.table.filter$cellType.target)
    
    overlap.table.filter <- overlap.table$marker==marker
    overlap.table[overlap.table.filter,]$groups.overlapping <- 
      paste0(groups.present, collapse = ";")
    overlap.table[overlap.table.filter,]$types.overlapping <- 
      paste0(types.present, collapse = ";")
    overlap.table[overlap.table.filter,]$number.overlapping.groups <- 
      length(groups.present)
    overlap.table[overlap.table.filter,]$number.overlapping.types <- 
      length(types.present)
    overlap.table[overlap.table.filter,]$overlap.group.rate <- 
      length(groups.present)/total.groups
  }
  return(overlap.table)
}



