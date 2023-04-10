#!/usr/bin/env R

# Author: Sean Maden
#
# Methods to select markers using typemarkers() and supported methods.
#

#' filter_group_markers
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
#' @returns List containing marker filter results, including the filtered 
#' markers list, the overlap info table, and filter metadata.
#' @details Performs the following marker filter steps:
#' * Uniqueness: Removes duplicate markers across type lists.
#' * Concordance: Removes markers which aren't specific to the same type across 
#' groups.
#' * Overlap: Removes markers below an overlap rate, defined using the 
#' \code{min.group.overlap.rate} argument.
#' @seealso \code{markers_by_group}, \code{marker_overlaps}
#' @examples
#' sce.example <- random_sce(num.cells = 100)
#' sce.example[["sample_id"]] <- c(rep("sample1", 10), rep("sample2", 80), rep("sample1", 10))
#' group.markers <- markers_by_group(sce.example)
#' filtered.group.markers <- filter_group_markers(group.markers)
#' @export
filter_group_markers <- function(group.markers,
                                 minimum.group.overlap.rate = 0.5,
                                 verbose = TRUE){
  require(dplyr)
  if(is(group.markers, "list")){
    markers.table <- do.call(rbind, group.markers) %>% as.data.frame()
  }
  num.markers.input <- markers.table$gene %>% unique()
  if(verbose){message("Found ", num.markers.input, " total markers.")}
  
  if(is(marker.overlaps, "NULL")){
    if(verbose){message("Getting marker overlap info...")}
    overlap.info <- marker_overlaps(group.markers)
  }
  
  markers.start.vector <- markers.table$gene
  filter.non.concordant <- grepl(";", overlap.info$types)
  markers.concordant <- overlap.info[filter.non.concordant,]
  num.concordant <- nrow(markers.concordant)
  if(verbose){message("Found ", num.concordant, " duplicated markers by type.")}
  
  filter.min.overlap <- markers.concordant$overlap.group.rate >= 
    minimum.group.overlap.rate
  markers.overlap <- markers.concordant[filter.min.overlap,]
  num.overlap <- nrow(markers.overlap)
  if(verbose){
    message("Found ", num.overlap, " markers with overlap rate at least ",
            min.overlap.rate,".")}
  
  # return.list
  metadata.list <- list(num.markers.input = num.markers.input,
                        num.concordant = num.concordant,
                        num.overlap = num.overlap,
                        min.overlap.rate = min.overlap.rate)
  return.list <- list(marker.list.final = markers.overlap$marker,
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
#' sce.example[["sample_id"]] <- c(rep("sample1", 10), rep("sample2", 80), rep("sample1", 10))
#' group.markers <- markers_by_group(sce.example)
#' group.overlaps <- marker_overlaps(group.markers, c("gene1", "gene20"))
#' @export
marker_overlaps <- function(group.markers, marker.filter.vector){
  require(dplyr)
  marker.filter.vector <- c("gene1", "gene20")
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
  unique.markers <- overlap.table$marker %>% unique()
  overlap.table <- data.frame(marker = unique.markers)
  overlap.table$groups.overlapping <- overlap.table$types.overlapping <- 
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
#' @param verbose Whether to show verbose status messages.
#' @returns List or table of markers (see \code{return.type} argument), 
#' organized by group IDs.
#' @details Gets the group markers for all individual specified groups, where
#' "group" here means donor, batch, or some other set of cells over which we
#' may be concerned about marker consistency.
#' @examples
#' sce.example <- random_sce(num.cells = 100, num.genes = 100)
#' sce.example[["sample_id"]] <- c(rep("sample1", 10), rep("sample2", 80), rep("sample1", 10))
#' group.markers <- markers_by_group(sce.example)
#' @seealso \code{filter_group_markers}
markers_by_group <- function(sce, 
                             group.variable = "sample_id", 
                             celltype.variable = "celltype", 
                             assay.name = "counts", 
                             markers.per.type = 20, 
                             typemarker.algorithm = "meanratios",
                             return.type = "list",
                             verbose = FALSE){
  require(dplyr)
  unique.group.vector <- sce[[group.variable]] %>% unique() %>% as.character()
  if(verbose){message("get gene markers for each specified batch")}
  group.markers.list <- lapply(unique.group.vector, function(group.id){
    if(verbose){message("getting markers for batch id: ", group.id, "...")}
    filter <- sce[[group.variable]] == group.id
    result.table <- lute(sce = sce[,filter], 
         celltype.variable = celltype.variable, 
         assay.name = assay.name, 
         markers.per.type = markers.per.type,
         typemarker.algorithm = typemarker.algorithm,
         deconvolution.algorithm = NULL,
         return.info = TRUE)$typemarker.results$result.info
    result.table$group.id <- group.id
    return(result.table)
  })
  names(group.markers.list) <- unique.group.vector
  if(return.type == "list"){
    return.data <- group.markers.list
    names(return.data)
  } else{
    return.data <- do.call(rbind, group.markers.list) %>% as.data.frame()
  }
  if(verbose){message("finished marker lists. Returning...")}
  return(return.data)
}
