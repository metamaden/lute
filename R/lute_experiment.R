#!/usr/bin/env R

# Author: Sean Maden
#
# Helper functions to conduct a deconvolution experiment using the lute framework.

#' deconvolution.experiment.info
#'
#' Get terms for a new deconvolution experiment, startin from an object of type 
#' \linkS4class{SingleCellExperiment}.
#'
#' @param y Bulk/convoluted signals matrix.
#' @param sce Object of type \linkS4class{SingleCellExperiment}.
#' @param z Signature matrix. If included, uses this for experiments instead of
#' group-wise signature matrix.
#' @param p.true.proportions Named numeric vector of true cell type proportions. 
#' If NULL, uses the proportions from the sample-wise sce data.
#' @param s Named numeric vector of cell factors. Used to make the pseudobulk, 
#' or used by deconvolution.results.table to transform the Z signature matrix.
#' @param assay.name Name of expression matrix in \code{sce} assays.
#' @param sample.id.variable Name of experiment group variable in \code{sce} 
#' coldata. Used to group separate deconvolution experiments in a series.
#' @param celltype.variable Name of cell type variable in \code{sce} coldata.
#' @param markers.vector Vector of marker IDs, used to filter the \code{sce} and 
#' \code{y} objects prior to deconvolution.
#' @param experiment.labels Character string vector of labels for the experiment.
#' @param verbose Whether to show verbose status messages.
#' 
#' @returns List of length equal to the number of unique sample ids, containing
#' metadata and key variables for a new deconvolution experiment series.
#' 
#' @examples
#' sce.example <- random_sce()
#' sce.example[["sample.id"]] <- "sample1"
#' experiment.info.list <- deconvolution.experiment.info(sce.example)
#' experiment.info.list$sample1$experiment.labels # ";pseudobulk;"
#' 
#' @seealso \code{deconvolution.results.table}, \code{deconvolution.experiment}
#' 
#' @export
deconvolution.experiment.info <- function(sce, s = NULL, z = NULL, 
                                          p.true.proportions = NULL,
                                          assay.name = "counts",
                                          sample.id.variable = "sample.id",
                                          celltype.variable = "celltype",
                                          markers.vector = NULL, y = NULL,
                                          experiment.labels = "", verbose = FALSE){
  require(dplyr)
  samples.id.vector <- sce[[sample.id.variable]]
  unique.samples <- unique(samples.id.vector)
  unique.cell.types <- unique(sce[[celltype.variable]])
  if(!is(markers.vector, "NULL")){sce <- sce[markers.vector,]}
  if(nrow(sce)==0){stop("Error, no markers remaining in provided sce.")}
  if(verbose){message("getting experiment info...")}
  experiment.list <- lapply(unique.samples, function(sample.id){
    if(verbose){message("making experiment info for group id: ", sample.id)}
    sce.filter <- sce[[sample.id.variable]]==sample.id
    sce.sample <- sce[,sce.filter]
    if(verbose){message("parsing z signature matrix...")}
    z.sample <- z
    if(is(z, "NULL")){
      z.sample <- signature_matrix_from_sce(sce.sample, 
                                         celltype.variable = celltype.variable, 
                                         summary.method = "mean", 
                                         assay.name = assay.name)
    }
    if(verbose){message("parsing y bulk convoluted signals matrix...")}
    if(is(y, "NULL")){
      y.sample <- ypb_from_sce(sce = sce.sample, assay.name = assay.name,
                               celltype.variable = celltype.variable, S = s) %>% 
        as.matrix()
      experiment.labels <- paste0(experiment.labels, ";pseudobulk;")
    } else{
      y.sample <- y %>% as.matrix() # y[,sample.id,drop = FALSE] %>% as.matrix()
    }
    if(verbose){message("parsing p.true.proportions...")}
    if(is(p.true.proportions, "NULL")){
      p.true.counts <- sce.sample[[celltype.variable]] %>% table()
      p.true.proportions <- p.true.counts %>% prop.table()
    } else{
      p.true.proportions <- p.true.proportions
    }
    return(list(sample.id = sample.id,
                sce.sample = sce.sample,
                z.sample = z.sample,
                y.sample = y.sample,
                s.experiment = s,
                p.true.counts = p.true.counts,
                p.true.proportions = p.true.proportions,
                experiment.labels = experiment.labels))
  })
  names(experiment.list) <- unique.samples
  return(experiment.list)
}

#' deconvolution.results.table
#'
#' Get the results table for a deconvolution experiment.
#' 
#' @param experiment.list An experiment info list output from \code{experiment_info}
#' @param deconvolution.algorithm Valid deconvolution algorithm supported by 
#' the \code{lute} framework (see \code{?lute} for details).
#' @param s Cell factors for Z signature matrix adjustment. If NULL, use the term 
#' s.experiment from experiment list sample data.
#' @param typemarker.algorithm Valid algorithm for marker selection supported 
#' by the \code{lute} framework, or NULL to skip this.
#' @param celltype.variable Name of cell type variable in \code{sce} and 
#' \code{z} coldata.
#' @param verbose Whether to show verbose status messages.
#' 
#' @returns Results table, summarizing results and metadata for a deconvolution 
#' experiment.
#' 
#' @details Performs a deconvolution experiment, aggregating the results into
#' a results table. Columns in the results table have controlled definitions
#' and are recognized by downstream functions. The column definitions are:
#' 
#' @examples 
#' sce.example <- random_sce()
#' sce.example[["sample.id"]] <- "sample1"
#' experiment.info.list <- deconvolution.experiment.info(sce.example)
#' results.table <- deconvolution.results.table(experiment.info.list)
#' 
#' @seealso \code{deconvolution.experiment.info}, \code{deconvolution.experiment}
#' 
#' @export
deconvolution.results.table <- function(experiment.list, 
                                        celltype.variable = "celltype",
                                        typemarker.algorithm = NULL,
                                        deconvolution.algorithm = "nnls", 
                                        s = NULL, verbose = FALSE){
  require(dplyr)
  if(verbose){message("getting deconvolution results...")}
  results.list <- lapply(experiment.list, function(sample.data){
    if(verbose){
      message("getting results for experiment group id: ", sample.data$sample.id)}
    if(is(s, "NULL")){s <- sample.data$s.experiment}
    list.sizes <- list(null = NULL, adjustment = s)
    results.list.sample <- lapply(list.sizes, function(size.data){
      results <- lute(y = sample.data$y.sample, 
                      sce = sample.data$sce.sample,
                      z = sample.data$z.sample, 
                      celltype.variable = celltype.variable,
                      s = size.data, return.info = FALSE, 
                      typemarker.algorithm = typemarker.algorithm,
                      deconvolution.algorithm = deconvolution.algorithm) %>% 
        unlist()
      if(sum(results) > 1){results <- results/sum(results)}
      return(results)
    })
    results.table.sample <- do.call(rbind, results.list.sample) %>% 
      as.data.frame()
    # format colnames
    colnames(results.table.sample) <- gsub(".*\\.", "", 
                                           colnames(results.table.sample))
    colnames(results.table.sample) <- paste0(colnames(results.table.sample), 
                                             ".predicted.proportion")
    # append size data names
    results.table.sample$cell.size.adjustment.type <- names(list.sizes)
    # append experiment metadata
    results.table.sample$sample.id <- sample.data$sample.id
    results.table.sample$deconvolution.algorithm <- deconvolution.algorithm
    # parse true cell amounts
    unique.cell.types <- names(sample.data$p.true.proportions)
    for(type in unique.cell.types){
      true.proportion <- sample.data$p.true.proportions[[type]]
      true.count <- sample.data$p.true.count[[type]]
      error <- results.table.sample[,paste0(type,".predicted.proportion")] -
        true.proportion
      abs.error <- error %>% abs()
      # append
      results.table.sample[,paste0(type,".true.proportion")] <- true.proportion
      results.table.sample[,paste0(type,".true.count")] <- true.count
      results.table.sample[,paste0(type,".error")] <- error
      results.table.sample[,paste0(type,".absolute.error")] <- abs.error
    }
    results.table.sample
  })
  results.table <- do.call(rbind, results.list)
  return(results.table)
}

#' deconvolution.experiment
#'
#' Performs a deconvolution experiment.
#' 
#' @param sce Object of type \linkS4class{SingleCellExperiment}.
#' @param y Bulk/convoluted signals matrix.
#' @param s Named numeric vector of cell factors. Used to make the pseudobulk, 
#' or used by deconvolution.results.table to transform the Z signature matrix.
#' @param p.true.proportions Named numeric vector of true cell type proportions. 
#' If NULL, uses the proportions from the sample-wise sce data.
#' @param z Signature matrix. If included, uses this for experiments instead of
#' group-wise signature matrix.
#' @param assay.name Name of expression matrix in \code{sce} assays.
#' @param sample.id.variable Name of experiment group variable in \code{sce} 
#' coldata. Used to group separate deconvolution experiments in a series.
#' @param celltype.variable Name of cell type variable in \code{sce} coldata.
#' @param markers.vector Vector of marker IDs, used to filter the \code{sce} and 
#' \code{y} objects prior to deconvolution.
#' @param experiment.labels Character string vector of labels for the experiment.
#' @param deconvolution.algorithm Valid deconvolution algorithm supported by 
#' the \code{lute} framework (see \code{?lute} for details).
#' @param typemarker.algorithm Valid algorithm for marker selection supported 
#' by the \code{lute} framework, or NULL to skip this.
#' @param verbose Whether to show verbose status messages.
#'
#' @details Perform a full deconvolution experiment.
#' 
#' @returns List of experiment groups (smaple IDs) containing variables and 
#' metadata for a deconvolution experiment series.
#'
#' @examples 
#' sce.example <- random_sce()
#' sce.example[["sample.id"]] <- "sample1"
#' full.experiment <- deconvolution.experiment(sce = sce.example, s = c(1, 10))
#' 
#' @seealso \code{deconvolution.results.summaries}, 
#' \code{deconvolution.experiment}, \code{plot.deconvolution.results}
#'
#' @export
deconvolution.experiment <- function(sce, y = NULL, s = NULL, z = NULL,
                                     p.true.proportions = NULL,
                                     assay.name = "counts",
                                     sample.id.variable = "sample.id", 
                                     celltype.variable = "celltype",
                                     deconvolution.algorithm = "nnls",
                                     typemarker.algorithm = NULL,
                                     markers.vector = NULL, 
                                     experiment.labels = "",
                                     verbose = FALSE){
  if(verbose){message("beginning deconvolution experiment...")}
  experiment.list <- deconvolution.experiment.info(
                      y = y, sce = sce, s = s, 
                      assay.name = assay.name, z = z,
                      sample.id.variable = sample.id.variable,
                      celltype.variable = celltype.variable,
                      markers.vector = markers.vector,
                      experiment.labels = experiment.labels,
                      verbose = verbose)
  results.table <- deconvolution.results.table(
    experiment.list = experiment.list, s = s,
    celltype.variable = celltype.variable,
    deconvolution.algorithm = deconvolution.algorithm,
    typemarker.algorithm = typemarker.algorithm)
  # results.summary <- results_summaries(results.table = results.table)
  plots.list <- deconvolution.results.plots(results.table = results.table,
                                            verbose = verbose)
  if(verbose){message("finished deconvolution experiment.")}
  return(list(experiment.info = experiment.list, results.table = results.table,
              plots.list = plots.list))
}

#' deconvolution.results.plots
#' 
#' Get ggplot2 plots summarizing the results of a deconvolution experiment.
#' 
#' 
#' @param results.table Experiment results table returned by 
#' \code{deconvolution.results.table()} (see 
#' \code{?deconvolution.results.table} for details).
#' @param verbose Whether to show verbose status messages.
#' 
#' @details Makes the following results summary plots:
#' 
#' * proportions.scatterplot : Scatter plot of the (x-axis) true and (y-axis) 
#' predicted cell type proportions.
#' 
#' * abs.error.barplot : Bar plots of (y-axis) absolute errors organized by 
#' (x-axis) experiment groups. 
#' 
#' * abs.error.jitterbox : Composite jittered points and boxplot overlay of
#' the (y-axis) absolute errors across experiment groups, organized by (x-axis) 
#' cell factor adjustment type.
#' 
#' @returns List of ggplot2 plots organized by cell types.
#' 
#' @examples 
#' sce.example <- random_sce()
#' sce.example[["sample.id"]] <- "sample1"
#' experiment.info.list <- deconvolution.experiment.info(sce.example)
#' results.table <- deconvolution.results.table(experiment.info.list)
#' results.plots <- deconvolution.results.plots(results.table)
#'
#' @seealso \code{deconvolution.results.table}, \code{deconvolution.experiment}
#'
#' @export
deconvolution.results.plots <- function(results.table, verbose = FALSE){
  require(ggplot2)
  if(verbose){message("getting results table plots...")}
  plots.list <- list()
  results.colnames <- colnames(results.table)
  colnames.filter <- grepl(".*\\.predicted\\.proportion$", results.colnames)
  unique.cell.types <- results.colnames[colnames.filter]
  unique.cell.types <- gsub("\\..*", "", unique.cell.types)
  plots.list <- lapply(unique.cell.types, function(type){
    results.type.filter <- grepl(type, results.colnames)
    results.type <- results.table[,results.type.filter]
    colnames(results.type) <- gsub(paste0(type, "\\."), "", colnames(results.type))
    results.type$adjustment.type <- results.table$cell.size.adjustment.type
    results.type$sample.id <- results.table$sample.id
    # proportions scatterplot
    proportions.scatterplot <- ggplot(results.type, 
                                      aes(x = true.proportion, 
                                          y = predicted.proportion)) + 
      geom_point() + 
      geom_abline(slope = 1, intercept = 0) + theme_bw() +
      ggtitle(type) + facet_wrap(~adjustment.type)
    # abs.error barplot
    abs.error.barplot <- ggplot(results.type, 
                                aes(x = sample.id, 
                                    y = absolute.error)) + 
      geom_bar(stat = "identity") + theme_bw() + ggtitle(type) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      facet_wrap(~adjustment.type)
    # abs.error jitterbox
    abs.error.jitterbox <- ggplot(results.type, 
                                  aes(x = adjustment.type, 
                                      y = absolute.error)) + 
      geom_jitter() + geom_boxplot(alpha = 0, color = 'cyan') + theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle(type)
    return(list(proportions.scatterplot = proportions.scatterplot,
                abs.error.barplot = abs.error.barplot,
                abs.error.jitterbox = abs.error.jitterbox))
  })
  names(plots.list) <- unique.cell.types
  return(plots.list)
}

#' deconvolution.experiment.permute.groups
#'
#' Run a series of permutation experiments.
#'
#' @param sce Object of type \linkS4class{SingleCellExperiment}.
#' @param s Named numeric vector of cell factors. Used to make the pseudobulk, 
#' or used by deconvolution.results.table to transform the Z signature matrix.
#' @param assay.name Name of expression matrix in \code{sce} assays.
#' @param sample.id.variable Name of experiment group variable in \code{sce} 
#' coldata. Used to group separate deconvolution experiments in a series.
#' @param celltype.variable Name of cell type variable in \code{sce} coldata.
#' @param summary.method Statistical method to make the signature matrix from 
#' an \code{sce} object.
#' @param experiment.labels Character string vector of labels for the experiment.
#' @param deconvolution.algorithm Valid deconvolution algorithm supported by 
#' the \code{lute} framework (see \code{?lute} for details).
#' @param verbose Whether to show verbose status messages.
#' @param ... Additional arguments passed to \code{deconvolution.experiment()} 
#' (see \code{?deconvolution.experiment} for details).
#'
#' @details Run a series of deconvolution experiments, permuting the Z signature 
#' matrix and Y pseudobulk signals according to the specified experiment group 
#' id.
#' 
#' The returned experiment results table is identical to the results table from
#' \code{deconvolution.experiment}, except it adds a new column 
#' "group.id.signature" for the group id corresponding to the signature matrix.
#' The column "sample.id" correspond to the group id of the pseudobulk data.
#'
#' @returns List containing the experiment results table and results plot list.
#'
#' @examples 
#' sce.example <- random_sce(num.genes = 100, num.cells = 100)
#' s <- c("type1" = 3, "type2" = 10)
#' sce.example[["sample.id"]] <- c(rep("sample1", 10), rep("sample2", 50), rep("sample1", 40))
#' experiment <- deconvolution.experiment.permute.groups(sce = sce.example, s = s)
#' 
#' @seealso \code{deconvolution.experiment}, 
#' \code{deconvolution.results.plots.permutations}
#' 
#' @export
deconvolution.experiment.permute.groups <- function(sce, s, 
                  sample.id.variable = "sample.id", 
                  celltype.variable = "celltype", summary.method = "mean", 
                  assay.name = "counts", deconvolution.algorithm = "nnls", 
                  experiment.labels = "permutation",
                  verbose = FALSE, ...){
  require(dplyr)
  if(verbose){message("running permutation experiments...")}
  group.id.vector <- sce[[sample.id.variable]]
  unique.group.id.vector <- group.id.vector %>% unique() %>% as.character()
  
  if(verbose){message("making pseudobulk table from unique group ids...")}
  ypb.list <- lapply(unique.group.id.vector, function(group.id){
    filter.group <- sce[[sample.id.variable]]==group.id
    ypb_from_sce(sce = sce[,filter.group], 
                 assay.name = assay.name, 
                 celltype.variable = celltype.variable, 
                 S = s)
  })
  ypb.table <- do.call(cbind, ypb.list) %>% as.matrix()
  colnames(ypb.table) <- unique.group.id.vector
  
  results.table.list <- lapply(unique.group.id.vector, 
                               function(group.id.z){
    if(verbose){message("working on group id ", group.id.z, "...")}
    if(verbose){message("getting main signature matrix...")}
    filter.sce <- sce[[sample.id.variable]]==group.id.z
    sce.group.z <- sce[,filter.sce]
    z <- signature_matrix_from_sce(sce = sce.group.z,
                                   celltype.variable = celltype.variable,
                                   summary.method = summary.method,
                                   assay.name = assay.name)
    
    if(verbose){message("getting pseudobulks...")}
    pb.filter <- !unique.group.id.vector==group.id.z
    unique.group.id.pb <- unique.group.id.vector[pb.filter] %>% unique()
    ypb.table.iteration <- ypb.table[, unique.group.id.pb, drop = F]
    
    if(verbose){message("getting experiment series...")}
    experiment <- deconvolution.experiment(sce = sce, 
                                           y = ypb.table.iteration, 
                                           s = s, z = z, 
                                           assay.name = assay.name, 
                                           sample.id.variable = 
                                             sample.id.variable, 
                                           experiment.labels = 
                                             paste0(experiment.labels, group.id.z),
                                           deconvolution.algorithm = 
                                             deconvolution.algorithm,
                                           celltype.variable = 
                                             celltype.variable,
                                           ...)
    if(verbose){message("finished with group id, returning results table...")}
    results.table.iteration <- experiment$results.table
    results.table.iteration$group.id.signature <- group.id.z
    results.table.iteration
  })
  results.table <- do.call(rbind, results.table.list) %>% as.data.frame()
  # make plots
  plots.list <- deconvolution.results.plots.permutations(results.table)
  return(
    list(
      results.table = results.table, plots.list = plots.list))
}

#' deconvolution.results.plots.permutations
#' 
#' Permutation experiment plots.
#' 
#' @param results.table Table of permutation experiment results.
#' @param verbose Whether to show verbose status messages.
#' 
#' @returns List containing experiment results plots.
#' 
#' @details Plot functions for permutation experiments. The \code{results.table}
#' should contain columns for "group.id.signature", for the signature matrix id, 
#' and "sample.id", for the pseudobulk sample id.
#' 
#' @examples 
#' sce.example <- random_sce(num.genes = 100, num.cells = 100)
#' s <- c("type1" = 3, "type2" = 10)
#' sce.example[["sample.id"]] <- c(rep("sample1", 10), rep("sample2", 50), rep("sample1", 40))
#' experiment <- deconvolution.experiment.permute.groups(sce = sce.example, s = s)
#' experiment$plots.list$`miss-matched`$type1$proportions.scatterplot
#' 
#' @seealso \code{deconvolution.experiment}, \code{plot.deconvolution.results}
#' 
#' @export
deconvolution.results.plots.permutations <- function(results.table, verbose = FALSE){
  require(ggplot2)
  if(verbose){"Making results plots for permutation experiments..."}
  filter.results <- results.table$sample.id==results.table$group.id.signature
  filter.types.vector <- c("matched", "miss-matched", "all")
  lgg.permute <- lapply(filter.types.vector, function(filter.type){
    if(filter.type=="matched"){
      filter.final <- filter.results
    } else if(filter.type=="miss-matched"){
      filter.final <- !filter.results
    } else{
      filter.final <- seq(nrow(results.table))
    }
    results.filtered <- results.table[filter.final,]
    results.filtered$sample.id <- results.filtered$group.id.signature
    lgg.filter <- deconvolution.results.plots(results.filtered)
    return(deconvolution.results.plots(results.filtered))
  })
  names(lgg.permute) <- filter.types.vector
  return(lgg.permute)
}
