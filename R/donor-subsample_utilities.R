#!/usr/bin/env R

# Author: Sean Maden
#
# Utilities to perform subsampling experiments on donors/samples/batches.
#


#' prepare_subsample_experiment
#'
#' @param sce SingleCellExperiment
#' @param iterations Total iterations to perform, per valid method specified in 
#' the methods argument.
#' @param groups.per.iteration Number of groups to sample per iteration. Total 
#' groups should exceed this value, otherwise data is treated as a single group.
#' @param method.vector Vector of valid deconvolution methods to test.
#' @param celltype.variable Variable containing cell type labels.
#' @param assay.name Name of assay in sce to use for pseudobulk.
#' @param group.variable Variable containing group labels for donor/sample/batch.
#' @param fraction.cells Fraction of cells, per group and type, to subsample in each iteration.
#' @param scale.factor Vector of scale factors corresponding to alphebetical order of
#' unique type labels in celltype.variable.
#' @param seed.num Random seed for computational reproducibility.
#' @param save.fnstem Character string stem to append to new saved filenames.
#' @param base.path Path of directory to store data for workflow run.
#' @param which.save Which items to save when running this function. Either "sce" for 
#' SingleCellExperiment, "tp" for true proportions, "ypb" for pseudobulk, or "wt"
#' for the `r-nf` compatible workflow table. These objects correspond to filepaths
#' in save.paths (see details).
#' @param save.names List of new file names to save.
#' @returns List of iterations data, saving new files as side-effect.
#' @details Prepares files for a subsampling experiment to evaluate the impact 
#' of a group-wise bias effect on deconvolution outcomes. Either real or simulated
#' data can be provided to the sce argument. This is then randomly subsampled 
#' where remaining arguments specify details for random iterations. 
#' 
#' The fraction of cells per type and group is specified by the fraction.cells 
#' argument. To ensure the same number of cells are sampled across iterations,
#' the final total of cells by type to sample is taken from the group-wise 
#' minima.
#' 
#' Since it is assumed you will use the same bulk data and true proportions to
#' evaluate iterations, these data are calculated once using the full sce object.
#' The pseudobulk dataset is calculated as the produce of the mean signature 
#' gene expression and the scale.factor cell size scale factor object. The true
#' cell type proportions correspond to this pseudobulk sample.
#'
#' The final output files are saved at the base.path specified in save.paths list.
#' These are suitable for use in an r-nf workflow, or using some other means of
#' parallelization of the iteration runs.
#'
#' @examples 
#' sce <- random_sce()
#' sce[["donor"]] <- "donor1"
#' iterations.list <- prepare_subsample_experiment(sce, 
#' celltype.variable = "celltype", group.variable = "donor")
#' names(iterations.list)
#' 
#' @export
prepare_subsample_experiment <- function(sce, scale.factor, iterations = 10, 
                                         groups.per.iteration = 3,
                                         method.vector = c("nnls", "music"),
                                         celltype.variable = "k2", 
                                         group.variable = "Sample",
                                         assay.name = "counts_adj",
                                         fraction.cells = 0.25,
                                         seed.num = 0,
                                         save.fnstem = "",
                                         which.save = c("sce", "tp", "ypb", "wt"),
                                         base.path = "data",
                                         save.names = list(sce.name = "sce.rda",
                                                           wt.name = "workflow-table.csv",
                                                           tp.name = "true-proportions.rda",
                                                           ypb.name = "ypb.rda",
                                                           li.name = "lindex.rda"),
                                         verbose = TRUE){
  set.seed(seed.num)
  
  # get metadata vectors
  cd <- colData(sce)
  groups.vector <- cd[,group.variable]
  celltype.vector <- cd[,celltype.variable]
  # get unique groups
  unique.groups <- unique(groups.vector)
  # get alphabetized celltype labels
  unique.types <- unique(celltype.vector)
  unique.types <- unique.types[order(unique.types)]
  message("Found ",length(unique.types), " cell type labels: ", 
          paste0(unique.types, collapse = ";"))
  # parse scale factors
  S <- scale.factor
  S <- S[names(S) %in% unique.types]
  S <- S[order(match(names(S), names(unique.types)))]
  message("Found ",length(S), " cell types with scale factors")
  
  # prepare iterations parameters
  if(groups.per.iteration > length(unique.groups)){
    groups.per.iteration <- "NULL"}
  # get the exact numbers of cells of each type to sample
  dft <- as.data.frame(table(celltype.vector, groups.vector))
  # get cells by group as minima
  if(fraction.cells > 1){fraction.cells = fraction.cells/100}
  message("Using cell fraction: ", fraction.cells)
  num.cells.vector <- sapply(unique.types, function(ti){
    type.filter <- dft[,1]==ti; dff <- dft[type.filter,]
    round(min(dff[,3])*fraction.cells, 0)
  })
  
  if(verbose){message("Getting random cell indices for iterations...")}
  lindex <- lapply(seq(iterations), function(ii){
    if(!is(groups.per.iteration, "NULL")){
      random.groups <-  sample(unique.groups, groups.per.iteration)
      cdf <- cd[cd[,group.variable] %in% random.groups,]
    } else{
      cdf <- cd
      random.groups = unique.groups
    }
    cell.index.vector <- unlist(lapply(unique.types, function(typei){
      type.filter <- cdf[,celltype.variable]==typei
      cell.id.vector <- rownames(cdf[type.filter,])
      message("From ", length(cell.id.vector), " getting ", num.cells.vector[typei], " cells.")
      sample(cell.id.vector, num.cells.vector[typei])
    }))
    vindex <- which(colnames(sce) %in% cell.index.vector)
    list(vindex = vindex, samples = random.groups) # return
  })
  
  if(verbose){message("Getting pseudobulk and true proportions...")}
  tp <- as.data.frame(table(sce[[celltype.variable]]))
  tp.prop <- tp[,2]; names(tp.prop) <- tp[,1]
  P <- tp.prop/sum(tp.prop)
  Z <- do.call(cbind, lapply(unique.types, function(typei){
    type.filter <- sce[[celltype.variable]]==typei
    scef <- sce[,type.filter]
    rowMeans(assays(scef)[[assay.name]])
  }))
  ZS <- sweep(Z, 2, S, "*")
  ypb <- t(t(P) %*% t(ZS))
  
  if(verbose){message("Parsing save filenames...")}
  if(!save.fnstem==""){
    save.names <- lapply(save.names, function(namei){
      name.str <- unlist(strsplit(namei, "\\."))
      new.name <- paste(name.str[1:length(name.str)-1], save.fnstem, sep = "_")
      new.name <- paste(new.name, name.str[length(name.str)], sep = ".")
      return(new.name)
    })
  }
  
  if(verbose){message("Getting new workflow table...")}
  wti <- data.frame(iterations_index = seq(iterations))
  wti$method <- methodi
  wti$sample_id <- unlist(lapply(lindex, function(li){
    paste0(li$samples, collapse = ";")}))
  wti$celltype_variable <- celltype.variable
  wti$assay_name <- assay.name
  # manage filepaths
  cnamev <- c("sce_filepath", "bulk_filepath", "list_index_filepath",
              "true_proportions_filepath")
  fpathv <- c(file.path(base.path, save.names[["sce.name"]]), 
              file.path(base.path, save.names[["ypb.name"]]),
              file.path(base.path, save.names[["li.name"]]), 
              file.path(base.path, save.names[["tp.name"]]))
  for(ii in seq(length(cnamev))){
    wti[,cnamev[ii]] <- paste0('"$launchDir/', fpathv[ii], '"')}
  wt <- do.call(rbind, lapply(method.vector, function(methodi){
    wti$method <- methodi; wti
  }))
  
  if(verbose){message("Saving new data...")}
  if("ypb" %in% which.save){
    ypb.fpath <- file.path(base.path, save.names[["ypb.name"]])
    save(ypb, file = ypb.fpath)
  }
  if("tp" %in% which.save){
    tp.fpath <- file.path(base.path, save.names[["tp.name"]])
    save(P, file = tp.fpath)
  }
  if("sce" %in% which.save){
    sce.fpath <- file.path(base.path, save.names[["sce.name"]])
    save(sce, file = sce.fpath)
  }
  if("lindex" %in% which.save){
    li.fpath <- file.path(base.path, save.names[["li.name"]])
    save(lindex, file = li.fpath)
  }
  if("wt" %in% which.save){
    wt.fpath <- file.path(base.path, save.names[["wt.name"]])
    write.csv(wt, file = wt.fpath, row.names = F)
  }
  message("Finished run prep. Returning iterations index list.")
  return(lindex)
}

#'
#'
#'
#'
analyze_subsample_results <- function(){
  
}

#'
#'
#'
#'
#'
subsample_summary <- function(){
  
}

#'
#'
#'
#'
#'
subsample_plots <- function(){
  
}

