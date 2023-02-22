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
#' @param methods Vector of valid deconvolution methods to test.
#' @param celltype.variable Variable containing cell type labels.
#' @param group.variable Variable containing group labels for donor/sample/batch.
#' @param fraction.cells Fraction of cells, per group and type, to subsample in each iteration.
#' @param number.of.groups Total groups to sample for each iteration.
#' @param scale.factor Vector of scale factors corresponding to alphebetical order of
#' unique type labels in celltype.variable.
#' @param seed.num Random seed for computational reproducibility.
#' @param which.save Which items to save when running this function. Either "sce" for 
#' SingleCellExperiment, "tp" for true proportions, "ypb" for pseudobulk, or "wt"
#' for the `r-nf` compatible workflow table. These objects correspond to filepaths
#' in save.paths (see details).
#' @param save.paths List of save path details for new files.
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
prepare_subsample_experiment <- function(sce, iterations = 10, 
                                         groups.per.iteration = 3,
                                         methods = c("nnls", "music"),
                                         celltype.variable = "k2", 
                                         group.variable = "Sample",
                                         fraction.cells = 25,
                                         number.of.groups = 3,
                                         scale.factor = c(10, 2),
                                         seed.num = 0,
                                         which.save = c("sce", "tp", "ypb", "wt"),
                                         save.paths = list(base.path = "./data/", 
                                                           sce.name = "sce.rda",
                                                           wt.name = "workflow-table.rda",
                                                           tp.name = "true-proportions.rda",
                                                           pb.name = "ypb.rda"),
                                         verbose = TRUE){
  set.seed(seed.num)
  
  # get metadata vectors
  cd <- colData(sce)
  groups.vector <- cd[,group.variable]
  celltype.vector <- cd[,celltype.variable]
  # get unique groups
  unique.groups <- unique(groups.vector)
  # get alphabetized celltype labels
  unique.type <- unique(celltype.vector)
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
  num.cells.vector <- sapply(unique.types, function(ti){
    round(min(dft[dft[,1]==ti,3])*fraction.cells, 0)})
  
  if(verbose){message("Getting random cell indices for iterations...")}
  lindex <- lapply(seq(num.iter), function(ii){
    if(!is(groups.per.iteration, "NULL")){
      random.groups <-  sample(unique.groups, groups.per.iteration)
      cdf <- cd[cd[,group.variable] %in% random.groups,]
    } else{
      cdf <- cd
      random.groups = unique.groups
    }
    cell.index.vector <- unlist(lapply(unique.types, function(typei){
      cell.id.vector <- rownames(cdf[cdf[,celltype.variable]==typei,])
      sample(cell.id.vector, num.cells.vector[typei])
    })))
    vindex <- which(colnames(sce) %in% cell.index.vector)
    list(vindex = vindex, samples = random.groups) # return
  })
  
  if(verbose){message("Getting pseudobulk and true proportions...")}
  tp <- as.data.frame(table(sce[[celltype.variable]]))
  tp.prop <- tp[,2]; names(tp.prop) <- tp[,1]
  P <- tp.prop/sum(tp.prop)
  Z <- do.call(cbind, lapply(unique.types, function(typei){
    rowMeans(assays(sce)[[assay.name]])
  }))
  ZS <- sweep(Z, 2, S, "*")
  ypb <- t(t(P) %*% t(ZS))
  
  if(verbose){message("Writing new workflow table...")}
  wt <- do.call(rbind, lapply(methods, function(methodi){
    wti <- data.frame(iterations_index = seq(num.iter))
    wti$method <- methodi
    wti$sample_id <- unlist(lapply(lindex, function(li){
      paste0(li$samples, collapse = ";")}))
    wti$celltype_variable <- celltype.variable
    wti$assay_name <- assay.name
    # manage filepaths
    cnamev <- c("sce_filepath", "bulk_filepath", "list_index_filepath",
                "true_proportions_filepath")
    fpathv <- c(file.path("data", save.paths[["sce.name"]]), 
                file.path("data", save.paths[["ypb.name"]]),
                file.path("data", save.paths[["li.name"]]), 
                file.path("data", save.paths[["tp.name"]]))
    for(ii in seq(length(cnamev))){
      wti[,cnamev[ii]] <- paste0('"$launchDir/', fpathv[ii], '"')}
    wti
  }))
  
  if(verbose){message("Saving new data...")}
  base.path <- sve.paths[["base.path"]]
  if("ypb" %in% which.save){
    ypb.fpath <- file.path(base.path, save.paths[["ypb.name"]])
    save(ypb, file = ypb.fpath)
  }
  if("tp" %in% which.save){
    tp.fpath <- file.path(base.path, save.paths[["tp.name"]])
    save(P, file = tp.fpath)
  }
  if("sce" %in% which.save){
    sce.fpath <- file.path(base.path, save.paths[["sce.name"]])
    save(sce, file = sce.fpath)
  }
  if("lindex" %in% which.save){
    li.fpath <- file.path(base.path, save.paths[["li.name"]])
    save(lindex, file = li.fpath)
  }
  if("wt" %in% which.save){
    wt.fpath <- file.path(base.path, save.paths[["wt.name"]])
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

